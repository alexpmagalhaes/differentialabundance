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
        .filter { it.transcript_length_matrix }
        .map { tools -> [ tools, file(tools.transcript_length_matrix, checkIfExists: true) ] }
        .ifEmpty { Channel.of([[], []]) }

    ch_control_features = ch_tools_with_id
        .filter { it.control_features }
        .map { tools -> [ tools, file(tools.control_features, checkIfExists: true) ] }
        .ifEmpty { Channel.of([[], []]) }

    ch_gene_sets = ch_tools_with_id
        .filter { it.gene_sets_files }
        .map { tools -> tools.gene_sets_files.split(",").collect { file(it, checkIfExists: true) } }
        .ifEmpty { Channel.of([[]]) }

    // Report related files
    ch_report_files = ch_tools_with_id
        .map { tools ->
            [
                file(tools.report_file, checkIfExists: true),
                file(tools.logo_file, checkIfExists: true),
                file(tools.css_file, checkIfExists: true),
                file(tools.citations_file, checkIfExists: true)
            ]
        }

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
            yml: extension == 'yml'
            csv: extension == 'csv'
            tsv: extension == 'tsv'
        }

    ch_contrasts_variables_from_yml = ch_contrast_variables_input.yml
        //.view()
        .map { tools, yaml_file, ext ->
            def yaml_data = new groovy.yaml.YamlSlurper().parse(yaml_file)
            yaml_data.contrasts.collect { contrast ->
                tuple('id': contrast.comparison[0], saved_meta = tools)
            }
        }
        .flatten()
        .unique() // Uniquify to keep each contrast variable only once (in case it exists in multiple lines for blocking etc.)

    ch_contrasts_variables_from_other = ch_contrast_variables_input.csv.splitCsv(header:true)
        .mix(ch_contrast_variables_input.tsv.splitCsv(header:true, sep:'\t'))
        .map { tools, row, ext ->
            ['id': row.variable, saved_meta: tools]
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
        .mix(PROTEUS.out.raw_tab.first().map{ meta, matrix -> tuple(meta.saved_meta, matrix) })

    // Normalised inputs

    ch_in_norm = AFFY_JUSTRMA_NORM.out.expression
        .mix(PROTEUS.out.norm_tab.first().map{ meta, matrix -> tuple(meta.saved_meta, matrix) })
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
