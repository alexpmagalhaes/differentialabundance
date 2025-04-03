/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { GUNZIP as GUNZIP_GTF                              } from '../modules/nf-core/gunzip/main'
include { UNTAR                                             } from '../modules/nf-core/untar/main.nf'
include { SHINYNGS_APP                                      } from '../modules/nf-core/shinyngs/app/main'
include { SHINYNGS_STATICEXPLORATORY as PLOT_EXPLORATORY    } from '../modules/nf-core/shinyngs/staticexploratory/main'
include { SHINYNGS_STATICDIFFERENTIAL as PLOT_DIFFERENTIAL  } from '../modules/nf-core/shinyngs/staticdifferential/main'
include { SHINYNGS_VALIDATEFOMCOMPONENTS as VALIDATOR       } from '../modules/nf-core/shinyngs/validatefomcomponents/main'
include { CUSTOM_MATRIXFILTER                               } from '../modules/nf-core/custom/matrixfilter/main'
include { ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION as GTF_TO_TABLE } from '../modules/nf-core/atlasgeneannotationmanipulation/gtf2featureannotation/main'
include { RMARKDOWNNOTEBOOK                                 } from '../modules/nf-core/rmarkdownnotebook/main'
include { AFFY_JUSTRMA as AFFY_JUSTRMA_RAW                  } from '../modules/nf-core/affy/justrma/main'
include { AFFY_JUSTRMA as AFFY_JUSTRMA_NORM                 } from '../modules/nf-core/affy/justrma/main'
include { PROTEUS_READPROTEINGROUPS as PROTEUS              } from '../modules/nf-core/proteus/readproteingroups/main'
include { GEOQUERY_GETGEO                                   } from '../modules/nf-core/geoquery/getgeo/main'
include { ZIP as MAKE_REPORT_BUNDLE                         } from '../modules/nf-core/zip/main'
include { IMMUNEDECONV                                      } from '../modules/nf-core/immunedeconv/main'
include { softwareVersionsToYAML                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'

//
// SUBWORKFLOW: Installed directly from nf-core/modules
//
include { ABUNDANCE_DIFFERENTIAL_FILTER                     } from '../subworkflows/nf-core/abundance_differential_filter/main'
include { DIFFERENTIAL_FUNCTIONAL_ENRICHMENT                } from '../subworkflows/nf-core/differential_functional_enrichment/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIFFERENTIALABUNDANCE {

    take:
    ch_tools

    main:

    ch_versions = Channel.empty()

    // ========================================================================
    // Handle input parameters
    // ========================================================================

    // Add id to meta once at the beginning
    ch_tools_with_id = ch_tools
        .map { tools -> tools + [id: tools.study_name] }

    // Get the sample sheets
    ch_samplesheet = ch_tools_with_id
        .map { tools ->
            [ tools, file(tools.input, checkIfExists: true) ]
        }

    // Create input channels based on study type
    ch_input = ch_samplesheet
        .branch { meta, input ->
            affy_array: meta.study_type == 'affy_array'
            maxquant: meta.study_type == 'maxquant'
            geo_soft_file: meta.study_type == 'geo_soft_file'
            rnaseq: meta.study_type == 'rnaseq'
        }

    // Handle Affy array inputs
    ch_celfiles = ch_input.affy_array
        .map { meta, input ->
            [ meta, file(meta.affy_cel_files_archive, checkIfExists: true) ]
        }

    // Handle Maxquant inputs
    ch_maxquant_inputs = ch_input.maxquant
        .map { meta, input ->
            matrix_file = file(meta.matrix, checkIfExists: true)
            [ meta, file(meta.input), matrix_file ]
        }

    // Handle GEO soft file inputs
    ch_querygse = ch_input.geo_soft_file
        .map { meta, input ->
            [ meta, meta.querygse ]
        }

    // Create optional parameter channels based on ch_tools
    ch_transcript_lengths = ch_tools_with_id
        .map { tools -> [ tools, tools.transcript_length_matrix ? file(tools.transcript_length_matrix, checkIfExists: true) : [] ] }

    ch_control_features = ch_tools_with_id
        .map { tools -> [ tools, tools.control_features ? file(tools.control_features, checkIfExists: true) : [] ] }

    ch_gene_sets = ch_tools_with_id
        .map { tools -> [ tools, tools.gene_sets_files ? tools.gene_sets_files.split(",").collect { file(it, checkIfExists: true) } : [] ] }

    // ========================================================================
    // Handle contrasts
    // ========================================================================

    // Create contrasts channels
    ch_contrasts_file = ch_tools_with_id
        .map { tools ->
            [ tools, file(tools.contrasts_yml ?: tools.contrasts, checkIfExists: true) ]
        }

    ch_contrasts_file_with_extension = ch_contrasts_file
        .map {
            tools, file -> [tools, file, file.extension]
        }

    ch_contrast_variables_input = ch_contrasts_file_with_extension
        .branch{ tools, file, extension ->
            yml: extension == 'yml' || extension == 'yaml'
            csv: extension == 'csv'
            tsv: extension == 'tsv'
        }

    ch_contrasts_variables_from_yml = ch_contrast_variables_input.yml
        .map { tools, yaml_file, ext ->
            def yaml_data = new groovy.yaml.YamlSlurper().parse(yaml_file)
            yaml_data.contrasts.collect { contrast ->
                tuple('id': contrast.comparison[0], 'study_meta': tools)
            }
        }
        .flatten()
        .unique() // Uniquify to keep each contrast variable only once (in case it exists in multiple lines for blocking etc.)

    ch_contrasts_variables_from_other = ch_contrast_variables_input.csv.splitCsv(header:true)
        .mix(ch_contrast_variables_input.tsv.splitCsv(header:true, sep:'\t'))
        .map { tools, row, ext ->
            ['id': row.variable, study_meta: tools]
        }
        .unique()

    ch_contrast_variables = ch_contrasts_variables_from_yml
        .mix(ch_contrasts_variables_from_other)

    // ========================================================================
    // Data preprocessing
    // ========================================================================

    // If we have affy array data in the form of CEL files we'll be deriving
    // matrix and annotation from them

    // Uncompress the CEL files archive
    UNTAR ( ch_celfiles )

    ch_affy_input = ch_input.affy_array
        .join(UNTAR.out.untar)

    // Run affy to derive the matrix. Reset the meta so it can be used to
    // define a prefix for different matrix flavours

    AFFY_JUSTRMA_RAW (
        ch_affy_input,
        [[],[]]
    )
    AFFY_JUSTRMA_NORM (
        ch_affy_input,
        [[],[]]
    )

    ch_affy_platform_features = AFFY_JUSTRMA_RAW.out.annotation

    ch_versions = ch_versions
        .mix(AFFY_JUSTRMA_RAW.out.versions)

    // We'll be running Proteus once per unique contrast variable to generate plots

    // Run proteus to import protein abundances
    PROTEUS(
        ch_contrast_variables.combine(ch_maxquant_inputs)
    )

    ch_versions = ch_versions.mix(PROTEUS.out.versions)

    // Run GEO query to get the annotation

    GEOQUERY_GETGEO(ch_querygse)
    ch_soft_features = GEOQUERY_GETGEO.out.annotation

    ch_versions = ch_versions
        .mix(GEOQUERY_GETGEO.out.versions)

    // Raw inputs
    // Note that we re-map the proteus output tables to the study ID as the
    // tables are the same across contrasts, only one norm table will be necessary

    ch_in_raw = ch_input.rnaseq
        .map{meta, file -> [meta, meta.matrix]}
        .mix(AFFY_JUSTRMA_RAW.out.expression)
        .mix(PROTEUS.out.raw_tab.first().map{ meta, matrix -> tuple(meta.study_meta, matrix) })

    // Normalised inputs

    ch_in_norm = AFFY_JUSTRMA_NORM.out.expression
        .mix(PROTEUS.out.norm_tab.first().map{ meta, matrix -> tuple(meta.study_meta, matrix) })
        .mix(GEOQUERY_GETGEO.out.expression)

    //// Fetch or derive a feature annotation table

    // Branch ch_tools_with_id based on feature source
    ch_feature_sources = ch_tools_with_id
        .branch {
            user_features: it.features
            affy_features: it.study_type == 'affy_array'
            geo_features: it.study_type == 'geo_soft_file'
            gtf_features: it.gtf
            matrix_features: true
        }

    // Handle user-provided feature annotations
    ch_user_features = ch_feature_sources.user_features
        .map { tools -> [ tools, file(tools.features, checkIfExists: true) ] }

    // Handle Affy array platform features
    ch_affy_features = ch_feature_sources.affy_features
        .map { tools -> ch_affy_platform_features }

    // Handle GEO soft file features
    ch_geo_features = ch_feature_sources.geo_features
        .map { tools -> ch_soft_features }

    // Handle GTF-based feature annotations
    ch_gtf_files = ch_feature_sources.gtf_features
        .map { tools -> [ tools, file(tools.gtf, checkIfExists: true) ] }

    // Process GTF files if needed
    ch_gtf_for_processing = ch_gtf_files
        .branch {
            compressed: it[1].name.endsWith('.gz')
            uncompressed: !it[1].name.endsWith('.gz')
        }

    // Decompress GTF files if needed
    GUNZIP_GTF(ch_gtf_for_processing.compressed)
    ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)

    // Combine compressed and uncompressed GTF files
    ch_gtf_processed = ch_gtf_for_processing.uncompressed
        .mix(GUNZIP_GTF.out.gunzip)

    // Convert GTF to feature annotation table
    GTF_TO_TABLE(ch_gtf_processed, [tuple('id':""), []])
    ch_gtf_features = GTF_TO_TABLE.out.feature_annotation
    ch_versions = ch_versions.mix(GTF_TO_TABLE.out.versions)

    ch_pre_matrix_features = ch_feature_sources.matrix_features.branch{ tools ->
        maxquant: tools.study_type == 'maxquant'
        other: true
    }
    ch_matrix_features = ch_pre_matrix_features.maxquant.join(ch_in_norm)
        .mix(ch_pre_matrix_features.other.join(ch_in_raw))
        .map{ meta, matrix ->
            matrix_as_anno_filename = "${workflow.workDir}/${matrix.getBaseName()}_as_anno.${matrix.getExtension()}"
            matrix_copy = file(matrix_as_anno_filename)
            matrix_copy.exists() && matrix.getText().md5().equals(matrix_copy.getText().md5()) ?: matrix.copyTo(matrix_as_anno_filename)
            [ meta, file(matrix_as_anno_filename) ]
        }

    // Combine all feature annotation channels
    ch_features = ch_user_features
        .mix(ch_affy_features)
        .mix(ch_geo_features)
        .mix(ch_gtf_features)
        .mix(ch_matrix_features)


    // Check compatibility of FOM elements and contrasts

    ch_matrix_sources = ch_tools_with_id.branch{
        affy_or_maxquant: params.study_type == 'affy_array' || params.study_type == 'maxquant'
        geo_soft_file: params.study_type == 'geo_soft_file'
        other: true
    }

    // Create a channel with the matrices for validation. For Affy and MaxQuant
    // we have both raw and norm matrices, for GEO soft file we only have norm,
    // for other we only have raw

    ch_matrices_for_validation = ch_matrix_sources.affy_or_maxquant
        .join(ch_in_raw)
        .join(ch_in_norm)
        .map{tuple(it[0], [it[1], it[2]])}
        .mix(ch_matrix_sources.geo_soft_file.join(ch_in_norm)
        .mix(ch_matrix_sources.other.join(ch_in_raw)))

    // Validate the components agains one another. We try not to assume too
    // much about chanel order, so we deploy the multiMap workaround.

    validator_input = ch_samplesheet
        .join(ch_matrices_for_validation)
        .join(ch_features)
        .join(ch_contrasts_file)
        .multiMap{meta, samplesheet, matrices, features_file, contrasts_file ->
            samplesheet_matrices: [meta, samplesheet, matrices]
            features_file: [meta, features_file]
            contrasts_file: [meta, contrasts_file]
        }

    VALIDATOR(
        validator_input.samplesheet_matrices,
        validator_input.features_file,
        validator_input.contrasts_file
    )

    // For Affy, we've validated multiple input matrices for raw and norm,
    // we'll separate them out again here

    ch_multi_validated_assays = VALIDATOR.out.assays
        .filter{meta, assays -> meta.study_type == 'affy_array' || meta.study_type == 'maxquant'}
        .transpose()
        .branch { meta, assay ->
            raw: assay.name.contains('raw')
            normalised: assay.name =~ /normali[sz]ed/
        }

    // Get raw matrices from the validation

    ch_raw = ch_multi_validated_assays.raw
        .mix(VALIDATOR.out.assays.filter{meta, assay -> meta.study_type== 'rnaseq'})

    // For RNASeq and GEO soft file we've validated a single matrix, raw in the
    // case of RNASeq and norm in the case of GEO soft file, and these are the
    // matrices we'll use for differential abundance testing. For affy or maxquant
    // we'll use the normalised matrices.

    ch_matrix_for_differential = ch_multi_validated_assays.normalised
        .mix(VALIDATOR.out.assays.filter{meta, assay -> meta.study_type == 'rnaseq' || meta.study_type == 'geo_soft_file'})

    // Split the contrasts up so we can run differential analyses and
    // downstream plots separately.
    // Replace NA strings that might have snuck into the blocking column

    ch_contrasts = VALIDATOR.out.contrasts
        .splitCsv ( header:true, sep:'\t' )
        .map{study_meta, contrast ->
            contrast.blocking = contrast.blocking.replaceAll('^NA$', '')
            if (!contrast.id){
                contrast.id = contrast.values().join('_')
            }
            meta = contrast + [study_meta: study_meta]
            tuple(study_meta, meta,contrast.variable, contrast.reference, contrast.target)
        }
        .groupTuple() // [study_meta, [meta], [variable], [reference], [target]]
        .map{ study_meta, meta, variable, reference, target -> [study_meta, [meta, variable, reference, target]]}

    // Firstly Filter the input matrix

    ch_matrixfilter_input = ch_matrix_for_differential
        .join(VALIDATOR.out.sample_meta)
        .multiMap{meta, matrix, samplesheet ->
            matrix: [meta, matrix]
            samplesheet: [meta, samplesheet]
        }

    CUSTOM_MATRIXFILTER(
        ch_matrixfilter_input.matrix,
        ch_matrixfilter_input.samplesheet
    )

    // ========================================================================
    // Differential analysis
    // ========================================================================

    // Prepare inputs for differential processes

    ch_differential_input = CUSTOM_MATRIXFILTER.out.filtered
        .map { meta, matrix ->
            [
                meta,
                matrix,
                meta.differential_method,
                meta.differential_min_fold_change,
                meta.differential_max_qval
            ]
        }
        .join(VALIDATOR.out.sample_meta)
        .join(ch_transcript_lengths)
        .join(ch_control_features)
        .join(ch_contrasts)
        .multiMap{meta, matrix, diff_method, diff_fc_threshold, diff_qval_threshold, samplesheet, transcript_lengths, control_features, contrasts ->
            input: [meta, matrix, diff_method, diff_fc_threshold, diff_qval_threshold]
            samplesheet: [meta, samplesheet]
            transcript_lengths: [meta, transcript_lengths]
            control_features: [meta, control_features]
            contrasts: [meta, contrasts]
        }

    // Run differential analysis

    ABUNDANCE_DIFFERENTIAL_FILTER(
        ch_differential_input.input,
        ch_differential_input.samplesheet,
        ch_differential_input.transcript_lengths,
        ch_differential_input.control_features,
        ch_differential_input.contrasts.map{ study_meta, contrasts -> contrasts }
    )

    // collect differential results

    ch_differential_results = ABUNDANCE_DIFFERENTIAL_FILTER.out.results_genewise//.view()
    ch_differential_results_filtered = ABUNDANCE_DIFFERENTIAL_FILTER.out.results_genewise_filtered
    ch_differential_model = ABUNDANCE_DIFFERENTIAL_FILTER.out.model
    ch_differential_norm = ABUNDANCE_DIFFERENTIAL_FILTER.out.normalised_matrix
    ch_differential_varstab = ABUNDANCE_DIFFERENTIAL_FILTER.out.variance_stabilised_matrix

    ch_versions = ch_versions
        .mix(ABUNDANCE_DIFFERENTIAL_FILTER.out.versions)

    // Derive a channel of normalised matrices
    // - from differential analysis for RNASeq
    // - from normalisedvalidated assays for Affy and MaxQuant
    // - from validated assays for GEO soft file

    ch_norm = ch_differential_norm
        .filter{meta, matrix -> meta.study_type == 'rnaseq'}
        .mix(ch_multi_validated_assays.normalised)
        .mix(VALIDATOR.out.assays.filter{meta, assay -> meta.study_type == 'geo_soft_file'})

    // Prepare channel with normalized matrix, and variance stabilized matrices when available
    ch_processed_matrices = ch_norm.join(ch_differential_varstab, remainder: true)
        .map { meta, norm, vs ->
            def matrices = vs ? [norm] + vs : [norm]
            [meta, matrices]
        }

    // ========================================================================
    // Functional analysis
    // ========================================================================

    // Prepare background file - for the moment it is only needed for gprofiler2

    ch_background = CUSTOM_MATRIXFILTER.out.filtered
        .filter{meta, matrix -> meta.functional_method == 'gprofiler2' && params.gprofiler2_background_file == "auto"}
        .mix(
            ch_tools_with_id
                .filter{tools -> tools.functional_method == 'gprofiler2' && ! tools.gprofiler2_background_file }
                .map{tools -> [tools, file(tools.gprofiler2_background_file, checkIfExists: true)]}
        )
        .mix(
            ch_tools_with_id
                .filter{ tools -> tools.functional_method != 'gprofiler2'}
                .map{tools -> [tools, []]}
        )

    // Prepare input for functional analysis

    ch_functional_analysis_matrices = ch_norm
        .filter{meta, matrix -> meta.functional_method == 'gsea'}
        .mix(
            ch_differential_results_filtered
                .filter{meta, input -> meta.functional_method == 'gprofiler2'}
        )

    ch_functional_input = ch_functional_analysis_matrices
        .join(ch_gene_sets)
        .join(ch_background)
        .join(ch_contrasts)
        .join(VALIDATOR.out.sample_meta)
        .join(VALIDATOR.out.feature_meta)
        .multiMap{meta, input, gene_sets, background, contrasts, samplesheet, features ->
            input: [meta, input, gene_sets, background, meta.functional_method]
            contrasts: [meta, contrasts]
            samplesheet: [meta, samplesheet]
            features: [meta, features, meta.features_id_col, meta.features_name_col]
        }

    // Run functional analysis

    DIFFERENTIAL_FUNCTIONAL_ENRICHMENT(
        ch_functional_input.input,
        ch_functional_input.contrasts.map{study_meta, contrasts -> contrasts},
        ch_functional_input.samplesheet,
        ch_functional_input.features
    )

    // Collect functional analysis results

    ch_gsea_results = DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gsea_report
    gprofiler2_plot_html = DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gprofiler2_plot_html
    gprofiler2_all_enrich = DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gprofiler2_all_enrich
    gprofiler2_sub_enrich = DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gprofiler2_sub_enrich

    ch_functional_results = gprofiler2_plot_html
        .join(gprofiler2_all_enrich, remainder: true)
        .join(gprofiler2_sub_enrich, remainder: true)
        .mix(ch_gsea_results)

    ch_versions = ch_versions
        .mix(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.versions)

    // ========================================================================
    // Other analyses
    // ========================================================================

    // Run IMMUNEDECONV

    ch_immunedeconv_input = ch_raw
        .filter{meta, raw -> meta.immunedeconv_run}
        .multiMap{meta, raw ->
            input: [meta, raw, meta.immunedeconv_method, meta.immunedeconv_function]
            name_col: meta.features_name_col
        }

    IMMUNEDECONV(
        ch_immunedeconv_input.input,
        ch_immunedeconv_input.name_col
    )

    ch_versions = ch_versions
        .mix(IMMUNEDECONV.out.versions)

    // ========================================================================
    // Plot figures
    // ========================================================================

    // For geoquery we've done no matrix processing and been supplied with the
    // normalised matrix, which can be passed through to downstream analysis

    ch_mat = ch_norm.filter{meta, matrix -> meta.study_type == 'geo_soft_file'}
        .map{meta, norm -> [meta, [norm]]}
        .mix(
            ch_raw
                .filter{meta, raw -> meta.study_type != 'geo_soft_file'}
                .join(ch_processed_matrices)
                .map { meta, raw, matrices ->
                    [meta, [raw] + matrices]
                }
        )

    ch_all_matrices = VALIDATOR.out.sample_meta              // meta_exp, samples
        .join(VALIDATOR.out.feature_meta)                    // meta_exp, samples, features
        .join(ch_mat)                                        // meta_mat, list of matrices (raw, norm, variance stabilized)
        .map { meta, samples, features, matrices ->
            return [meta, samples, features, matrices]
        }

    // Exploratory analysis

    ch_exploratory_input = ch_contrast_variables         // [meta]
        .map{meta -> [ meta.study_meta, meta.id] }       // [study_meta, variable]
        .groupTuple()                                    // [study_meta, [variable]]
        .combine(ch_all_matrices.groupTuple(), by: 0)    // [study_meta, [variable], [samples, features, [matrices]]]
        .transpose()                                     // [study_meta, variable, samples, features, [matrices]]
        .map{study_meta, variable, samples, features, matrices ->
            [['id': variable] + ['study_meta': study_meta], samples, features, matrices]
        }

    PLOT_EXPLORATORY(
        ch_exploratory_input
    )

    // Plot differential analysis results

    ch_plot_differential_input = ch_differential_results
        .map{meta, differential_results -> [meta.study_meta, meta, differential_results]}
        .combine(ch_all_matrices, by: 0)
        .multiMap{study_meta, meta, differential_results, samples, features, matrices ->
            differential_results: [meta, differential_results]
            samples_features_matrices: [meta, samples, features, matrices]
        }

    PLOT_DIFFERENTIAL(
        ch_plot_differential_input.differential_results,
        ch_plot_differential_input.samples_features_matrices
    )

    // Gather software versions

    ch_versions = ch_versions
        .mix(VALIDATOR.out.versions)
        .mix(PLOT_EXPLORATORY.out.versions)
        .mix(PLOT_DIFFERENTIAL.out.versions)

    // Gather software versions

    ch_versions = ch_versions
        .mix(VALIDATOR.out.versions)
        .mix(PLOT_EXPLORATORY.out.versions)
        .mix(PLOT_DIFFERENTIAL.out.versions)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
