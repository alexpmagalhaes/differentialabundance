/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { prepareModuleInput     } from '../subworkflows/local/utils_nfcore_differentialabundance_pipeline/main'
include { prepareModuleOutput    } from '../subworkflows/local/utils_nfcore_differentialabundance_pipeline/main'
include { getMetaWithoutContrast } from '../subworkflows/local/utils_nfcore_differentialabundance_pipeline/main'

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
    ch_paramsets

    main:

    ch_versions = Channel.empty()
    ch_paramsets_by_name = ch_paramsets.map{paramset -> [paramset.analysis_name, paramset]}

    // ========================================================================
    // Handle input
    // ========================================================================

    // Get the sample sheets
    ch_samplesheet = ch_paramsets
        .map { paramset ->
            [ paramset, file(paramset.input, checkIfExists: true) ]
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
            [ meta, input, file(meta.matrix, checkIfExists: true) ]
        }

    // Handle GEO soft file inputs
    ch_querygse = ch_input.geo_soft_file
        .map { meta, input ->
            [ meta, meta.querygse ]
        }

    // Create optional parameter channels based on ch_paramsets
    ch_transcript_lengths = ch_paramsets
        .map { paramset -> [ paramset, paramset.transcript_length_matrix ? file(paramset.transcript_length_matrix, checkIfExists: true) : [] ] }

    ch_control_features = ch_paramsets
        .map { paramset -> [ paramset, paramset.control_features ? file(paramset.control_features, checkIfExists: true) : [] ] }

    ch_gene_sets = ch_paramsets
        .map { paramset -> [ paramset, paramset.gene_sets_files ? paramset.gene_sets_files.split(",").collect { file(it, checkIfExists: true) } : [] ] }

    // ========================================================================
    // Handle contrasts
    // ========================================================================

    // Create contrasts channels
    ch_contrasts_file = ch_paramsets
        .map { paramset ->
            [ paramset, file(paramset.contrasts_yml ?: paramset.contrasts, checkIfExists: true) ]
        }

    ch_contrasts_file_with_extension = ch_contrasts_file
        .map {
            paramset, file -> [paramset, file, file.extension]
        }

    ch_contrast_variables_input = ch_contrasts_file_with_extension
        .branch{ paramset, file, extension ->
            yml: extension == 'yml' || extension == 'yaml'
            csv: extension == 'csv'
            tsv: extension == 'tsv'
        }

    ch_contrasts_variables_from_yml = ch_contrast_variables_input.yml
        .map { paramset, yaml_file, ext ->
            def yaml_data = new groovy.yaml.YamlSlurper().parse(yaml_file)
            yaml_data.contrasts.collect { contrast ->
                tuple(paramset, contrast.comparison[0])
            }
        }
        .map{it[0]}
        .unique() // Uniquify to keep each contrast variable only once (in case it exists in multiple lines for blocking etc.)

    ch_contrasts_variables_from_other = ch_contrast_variables_input.csv.splitCsv(header:true)
        .mix(ch_contrast_variables_input.tsv.splitCsv(header:true, sep:'\t'))
        .map { paramset, row, ext ->
            [paramset, row.variable]
        }
        .unique()

    ch_contrast_variables = ch_contrasts_variables_from_yml
        .mix(ch_contrasts_variables_from_other)

    // ========================================================================
    // Data type specific preprocessing
    // ========================================================================

    //
    // 1. deal with affy array data
    // If we have affy array data in the form of CEL files we'll be deriving
    // matrix and annotation from them
    //

    // Uncompress the CEL files archive

    ch_untar_input = prepareModuleInput(ch_celfiles)

    UNTAR ( ch_untar_input )

    ch_untar_out = prepareModuleOutput(UNTAR.out.untar, ch_paramsets_by_name)

    // run affy

    ch_affy_input = ch_input.affy_array
        .join(ch_untar_out)
    ch_affy_input = prepareModuleInput(ch_affy_input, 'affy', size=3)

    AFFY_JUSTRMA_RAW (
        ch_affy_input,
        [[],[]]
    )
    AFFY_JUSTRMA_NORM (
        ch_affy_input,
        [[],[]]
    )

    ch_affy_raw = prepareModuleOutput(AFFY_JUSTRMA_RAW.out.expression, ch_paramsets_by_name)
    ch_affy_norm = prepareModuleOutput(AFFY_JUSTRMA_NORM.out.expression, ch_paramsets_by_name)
    ch_affy_platform_features = prepareModuleOutput(AFFY_JUSTRMA_RAW.out.annotation, ch_paramsets_by_name)

    ch_versions = ch_versions
        .mix(AFFY_JUSTRMA_RAW.out.versions)

    //
    // 2. deal with maxquant data
    // We'll be running Proteus once per unique contrast variable to generate plots
    //

    ch_proteus_input = prepareModuleInput(ch_maxquant_inputs, 'proteus', size=3)

    // add contrast variable
    ch_proteus_input = ch_proteus_input
        .join(ch_contrast_variables)
        .map { meta, input, matrix, contrast ->
            def meta_new = meta + [contrast: contrast]
            [meta_new, input, matrix]
        }

    // Run proteus to import protein abundances
    PROTEUS( ch_proteus_input )

    // Note that the tables are the same across contrasts, only one table will be necessary
    // that is why here we take the first one and remove the contrast variable from meta
    ch_proteus_raw = prepareModuleOutput(PROTEUS.out.raw_tab, ch_paramsets_by_name)
        .first()
        .map { meta, matrix ->
            def meta_new = meta.findAll { key, value -> key != 'contrast' }
            [meta_new, matrix]
        }
    ch_proteus_norm = prepareModuleOutput(PROTEUS.out.norm_tab, ch_paramsets_by_name)
        .first()
        .map { meta, matrix ->
            def meta_new = meta.findAll { key, value -> key != 'contrast' }
            [meta_new, matrix]
        }

    ch_versions = ch_versions.mix(PROTEUS.out.versions)

    //
    // 3. deal with GEO soft file data
    // Run GEO query to get the annotation
    //

    ch_geoquery_input = prepareModuleInput(ch_querygse, 'geoquery')

    GEOQUERY_GETGEO(ch_geoquery_input)

    ch_soft_norm = prepareModuleOutput(GEOQUERY_GETGEO.out.expression, ch_paramsets_by_name)
    ch_soft_features = prepareModuleOutput(GEOQUERY_GETGEO.out.annotation, ch_paramsets_by_name)

    ch_versions = ch_versions
        .mix(GEOQUERY_GETGEO.out.versions)

    // ========================================================================
    // Parse channels for input matrices and features
    // ========================================================================

    //
    // Define raw and normalised input expression/abundance data
    //

    // Raw inputs

    ch_in_raw = ch_input.rnaseq
        .map{meta, file -> [meta, meta.matrix]}
        .mix(ch_affy_raw)
        .mix(ch_proteus_raw)

    // Normalised inputs

    ch_in_norm = ch_affy_norm
        .mix(ch_proteus_norm)
        .mix(ch_soft_norm)

    //
    // Fetch or derive a feature annotation table
    //

    // Branch ch_paramsets based on feature source
    ch_feature_sources = ch_paramsets
        .branch {
            user_features: it.features
            affy_features: it.study_type == 'affy_array'
            geo_features: it.study_type == 'geo_soft_file'
            gtf_features: it.gtf
            matrix_features: true
        }

    // Handle user-provided feature annotations
    ch_user_features = ch_feature_sources.user_features
        .map { paramset -> [ paramset, file(paramset.features, checkIfExists: true) ] }

    // Handle Affy array platform features
    ch_affy_features = ch_feature_sources.affy_features
        .map { paramset -> ch_affy_platform_features }

    // Handle GEO soft file features
    ch_geo_features = ch_feature_sources.geo_features
        .map { paramset -> ch_soft_features }

    // Handle GTF-based feature annotations
    ch_gtf_files = ch_feature_sources.gtf_features
        .map { paramset -> [ paramset, file(paramset.gtf, checkIfExists: true) ] }

    // Process GTF files if needed
    ch_gtf_for_processing = ch_gtf_files
        .branch {
            compressed: it[1].name.endsWith('.gz')
            uncompressed: !it[1].name.endsWith('.gz')
        }

    // Decompress GTF files if needed
    ch_gunzip_input = prepareModuleInput(ch_gtf_for_processing.compressed)
    GUNZIP_GTF( ch_gunzip_input )
    ch_gunzip_out = prepareModuleOutput(GUNZIP_GTF.out.gunzip, ch_paramsets_by_name)
    ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)

    // Combine compressed and uncompressed GTF files
    ch_gtf_processed = ch_gtf_for_processing.uncompressed
        .mix(ch_gunzip_out)

    // Convert GTF to feature annotation table
    ch_gtf_input = prepareModuleInput(ch_gtf_processed, 'gtf')
    GTF_TO_TABLE(
        ch_gtf_input,
        [tuple('id':""), []]
    )
    ch_gtf_features = prepareModuleOutput(GTF_TO_TABLE.out.feature_annotation, ch_paramsets_by_name)
    ch_versions = ch_versions.mix(GTF_TO_TABLE.out.versions)

    // extract features from matrix
    // note that in the case of maxquant we use the normalised matrix
    // whereas for other study types we use the raw matrix
    ch_pre_matrix_features = ch_feature_sources.matrix_features.branch{ paramset ->
        maxquant: paramset.study_type == 'maxquant'
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

    // ========================================================================
    // Validate input
    // ========================================================================

    // Check compatibility of FOM elements and contrasts

    ch_matrix_sources = ch_paramsets.branch{
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
        .mix(ch_matrix_sources.geo_soft_file.join(ch_in_norm))
        .mix(ch_matrix_sources.other.join(ch_in_raw))

    // Validate the components against one another. We try not to assume too
    // much about chanel order, so we deploy the multiMap workaround.

    validator_input = ch_samplesheet
        .join(ch_matrices_for_validation)
        .join(ch_features)
        .join(ch_contrasts_file)

    validator_input = prepareModuleInput(validator_input, 'validator', size=5)
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

    ch_validated_assays = prepareModuleOutput(VALIDATOR.out.assays, ch_paramsets_by_name)
    ch_validated_contrast = prepareModuleOutput(VALIDATOR.out.contrasts, ch_paramsets_by_name)
    ch_validated_samplemeta = prepareModuleOutput(VALIDATOR.out.sample_meta, ch_paramsets_by_name)
    ch_validated_featuremeta = prepareModuleOutput(VALIDATOR.out.feature_meta, ch_paramsets_by_name)

    // For Affy, we've validated multiple input matrices for raw and norm,
    // we'll separate them out again here

    ch_multi_validated_assays = ch_validated_assays
        .filter{meta, assays -> meta.study_type == 'affy_array' || meta.study_type == 'maxquant'}
        .transpose()
        .branch { meta, assay ->
            raw: assay.name.contains('raw')
            normalised: assay.name =~ /normali[sz]ed/
        }

    // Get raw matrices from the validation

    ch_raw = ch_multi_validated_assays.raw
        .mix(ch_validated_assays.filter{meta, assay -> meta.study_type == 'rnaseq'})

    // For RNASeq and GEO soft file we've validated a single matrix, raw in the
    // case of RNASeq and norm in the case of GEO soft file, and these are the
    // matrices we'll use for differential abundance testing. For affy or maxquant
    // we'll use the normalised matrices.

    ch_matrix_for_differential = ch_multi_validated_assays.normalised
        .mix(ch_validated_assays.filter{meta, assay -> meta.study_type == 'rnaseq' || meta.study_type == 'geo_soft_file'})

    // Split the contrasts up so we can run differential analyses and
    // downstream plots separately.
    // Replace NA strings that might have snuck into the blocking column

    ch_contrasts = ch_validated_contrast
        .splitCsv ( header:true, sep:'\t' )
        .map{paramset, contrast ->
            contrast.blocking = contrast.blocking.replaceAll('^NA$', '')
            if (!contrast.id){
                contrast.id = contrast.values().join('_')
            }
            return [paramset, contrast, contrast.variable, contrast.reference, contrast.target]
        }
        .groupTuple() // [paramset, [contrast], [variable], [reference], [target]]

    // ========================================================================
    // Filter matrix
    // ========================================================================

    ch_matrixfilter_input = ch_matrix_for_differential
        .join(ch_validated_samplemeta)

    ch_matrixfilter_input = prepareModuleInput(ch_matrixfilter_input, 'matrixfilter', size=3)
        .multiMap{meta, matrix, samplesheet ->
            matrix: [meta, matrix]
            samplesheet: [meta, samplesheet]
        }

    CUSTOM_MATRIXFILTER(
        ch_matrixfilter_input.matrix,
        ch_matrixfilter_input.samplesheet
    )

    ch_filtered_matrix = prepareModuleOutput(CUSTOM_MATRIXFILTER.out.filtered, ch_paramsets_by_name)

    // ========================================================================
    // Differential analysis
    // ========================================================================

    ch_differential_input = ch_filtered_matrix
        .join(ch_validated_samplemeta)
        .join(ch_transcript_lengths)
        .join(ch_control_features)
        .join(ch_contrasts)

    // Use a multiMap to generate synched channels for differential analysis
    ch_differential_input = prepareModuleInput(ch_differential_input, 'differential', size=9)
        .multiMap{ meta, matrix, samplesheet, transcript_lengths, control_features, contrast, variable, reference, target ->
            input: [meta, matrix, meta.differential_method, meta.differential_min_fold_change, meta.differential_max_qval]
            samplesheet: [meta, samplesheet]
            transcript_lengths: [meta, transcript_lengths]
            control_features: [meta, control_features]
            contrast: [meta, contrast, variable, reference, target]
        }

    // Run differential analysis

    ABUNDANCE_DIFFERENTIAL_FILTER(
        ch_differential_input.input,
        ch_differential_input.samplesheet,
        ch_differential_input.transcript_lengths,
        ch_differential_input.control_features,
        ch_differential_input.contrast
    )

    // Collect differential results

    // note that these channels, the meta contain contrast info too
    ch_differential_results = prepareModuleOutput(ABUNDANCE_DIFFERENTIAL_FILTER.out.results_genewise, ch_paramsets_by_name)
    ch_differential_results_filtered = prepareModuleOutput(ABUNDANCE_DIFFERENTIAL_FILTER.out.results_genewise_filtered, ch_paramsets_by_name)
    ch_differential_model = prepareModuleOutput(ABUNDANCE_DIFFERENTIAL_FILTER.out.model, ch_paramsets_by_name)

    // whereas these channels, the meta do not contain contrast info, as they come from the NORM modules instead of DIFFERENTIAL modules
    ch_differential_norm = prepareModuleOutput(ABUNDANCE_DIFFERENTIAL_FILTER.out.normalised_matrix, ch_paramsets_by_name)
    ch_differential_varstab = prepareModuleOutput(ABUNDANCE_DIFFERENTIAL_FILTER.out.variance_stabilised_matrix, ch_paramsets_by_name)

    ch_versions = ch_versions
        .mix(ABUNDANCE_DIFFERENTIAL_FILTER.out.versions)

    // Derive a channel of normalised matrices
    // - from differential analysis for RNASeq
    // - from normalisedvalidated assays for Affy and MaxQuant
    // - from validated assays for GEO soft file

    ch_norm = ch_differential_norm
        .filter{meta, matrix -> meta.study_type == 'rnaseq'}
        .mix(ch_multi_validated_assays.normalised)
        .mix(ch_validated_assays.filter{meta, assay -> meta.study_type == 'geo_soft_file'})

    // Prepare channel with normalized matrix, and variance stabilized matrices when available

    ch_processed_matrices = ch_norm.join(ch_differential_varstab, remainder: true)
        .map { meta, norm, vs ->
            def matrices = vs ? [norm, vs] : [norm]
            [meta, matrices]
        }

    // ========================================================================
    // Functional analysis
    // ========================================================================

    // Prepare background file - for the moment it is only needed for gprofiler2
    // The background file might come from:
    // - the filter matrix, if auto
    // - the gprofiler2 background file, if provided
    // - empty, if not gprofiler2

    ch_background = ch_filtered_matrix
        .filter{meta, matrix -> meta.functional_method == 'gprofiler2' && params.gprofiler2_background_file == "auto"}
        .mix(
            ch_paramsets
                .filter{paramset -> paramset.functional_method == 'gprofiler2' && paramset.gprofiler2_background_file }
                .map{paramset -> [paramset, file(paramset.gprofiler2_background_file, checkIfExists: true)]}
        )
        .mix(
            ch_paramsets
                .filter{ paramset -> paramset.functional_method != 'gprofiler2'}
                .map{paramset -> [paramset, []]}
        )

    // Prepare input for functional analysis
    // - use normalized matrix, if method is gsea
    // - use filtered differential results, if method is gprofiler2

    ch_functional_analysis_matrices = ch_norm
        .filter{meta, matrix -> meta.functional_method == 'gsea'}
        .map { meta, matrix ->
            def paramset = meta  // in this case paramset is the same as meta
            [paramset, meta, matrix]
        }
        .mix(
            ch_differential_results_filtered
                .filter{meta, results -> meta.functional_method == 'gprofiler2'}
                // note that differential analysis results contain the contrast info in the meta
                // we need to remove it to allow proper join in the next step
                .map{ meta, file ->
                    def paramset = getMetaWithoutContrast(meta)
                    [paramset, meta, file]
                }
        )

    ch_functional_input = ch_functional_analysis_matrices
        .join(ch_gene_sets)  // [meta, [gmt files]]
        .join(ch_background)
        .join(ch_contrasts)
        .join(ch_validated_samplemeta)
        .join(ch_validated_featuremeta)
        .map { it.tail() }  // remove paramset, and keep meta

    ch_functional_input = prepareModuleInput(ch_functional_input, 'functional', size=10)
        .multiMap{meta, input, gene_sets, background, contrasts, variable, reference, target, samplesheet, features ->
            input: [meta, input, gene_sets, background, meta.functional_method]
            contrasts: [meta, contrasts, variable, reference, target]
            samplesheet: [meta, samplesheet]
            features: [meta, features, meta.features_id_col, meta.features_name_col]
        }

    // Run functional analysis

    DIFFERENTIAL_FUNCTIONAL_ENRICHMENT(
        ch_functional_input.input,
        ch_functional_input.contrasts,
        ch_functional_input.samplesheet,
        ch_functional_input.features
    )

    // Collect functional analysis results

    // note that gsea results, as they were derived from the normalised matrix, meta does not contain contrast info
    ch_gsea_results = prepareModuleOutput(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gsea_report, ch_paramsets_by_name)

    // note that gprofiler2 results, as they were derived from the differential results, meta contains contrast info
    gprofiler2_plot_html = prepareModuleOutput(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gprofiler2_plot_html, ch_paramsets_by_name)
    gprofiler2_all_enrich = prepareModuleOutput(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gprofiler2_all_enrich, ch_paramsets_by_name)
    gprofiler2_sub_enrich = prepareModuleOutput(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gprofiler2_sub_enrich, ch_paramsets_by_name)

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

    ch_immunedeconv_input = ch_in_raw
        .filter{meta, raw -> meta.immunedeconv_run}

    ch_immunedeconv_input = prepareModuleInput(ch_immunedeconv_input, 'immunedeconv')
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

    ch_all_matrices = ch_validated_samplemeta                // meta, samples
        .join(ch_validated_featuremeta)                      // meta, features
        .join(ch_mat)                                        // meta, list of matrices (raw, norm, variance stabilized)
        .map { meta, samples, features, matrices ->
            return [meta, samples, features, matrices]
        }

    // Exploratory analysis

    ch_exploratory_input = ch_contrast_variables         // [meta, variable]
        .combine(ch_all_matrices, by:0)                  // [meta, variable, samples, features, [matrices]]
        .map { meta, variable, samples, features, matrices ->
            // we need to add variable info as id in meta
            // since PLOT_EXPLORATORY uses meta.id as the contrast variable
            def meta_new = meta + [id: variable]
            [meta_new, samples, features, matrices]
        }

    ch_exploratory_input = prepareModuleInput(ch_exploratory_input, 'exploratory', size=5)

    PLOT_EXPLORATORY(
        ch_exploratory_input
    )

    // Plot differential analysis results

    ch_plot_differential_input = ch_differential_results
        .map { meta, differential_results ->
            // differential_results meta contain contrast info
            // whereas ch_all_matrices meta do not
            // we need to remove the contrast info from meta
            def paramset = getMetaWithoutContrast(meta)
            [paramset, meta, differential_results]
        }
        .combine(ch_all_matrices, by: 0)
        .map{it.tail()}  // remove paramset, and keep meta

    ch_plot_differential_input = prepareModuleInput(ch_plot_differential_input, 'plot_differential', size=5)
        .multiMap{meta, differential_results, samples, features, matrices ->
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

    // ========================================================================
    // Generate report
    // ========================================================================

    // Collate and save software versions

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'differentialabundance_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    ch_collated_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'collated_versions.yml', sort: true, newLine: true)

    // Derive the report file

    ch_report_file = ch_paramsets
        .map{ paramset -> tuple(paramset, paramset.report_file) }

    // Generate a list of files that will be used by the markdown report

    ch_report_files = ch_paramsets
        .map { paramset ->
            [ paramset, [
                file(paramset.report_file, checkIfExists: true),
                file(paramset.logo_file, checkIfExists: true),
                file(paramset.css_file, checkIfExists: true),
                file(paramset.citations_file, checkIfExists: true)
            ]]
        }
        .join(ch_all_matrices)
        .join(ch_validated_contrast)
        .join(
            ch_differential_results
                .map { meta, differential_results ->
                    def paramset = getMetaWithoutContrast(meta)
                    [paramset, differential_results]
                }
                .groupTuple()
        )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
