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

    // Create input channels
    def exp_meta = [ "id": params.study_name  ]
    if (params.input) { ch_input = Channel.of([ exp_meta, file(params.input, checkIfExists: true) ]) }

    if (params.study_type == 'affy_array') {
        ch_celfiles = Channel.of([ exp_meta, file(params.affy_cel_files_archive, checkIfExists: true) ])
    } else if (params.study_type == 'maxquant') {
        matrix_file = file(params.matrix, checkIfExists: true)
        proteus_in = Channel.of([ file(params.input), matrix_file ])
    } else if (params.study_type == 'geo_soft_file') {
        ch_querygse = Channel.of([exp_meta, params.querygse])
    } else if (params.study_type == "rnaseq") {
        matrix_file = file(params.matrix, checkIfExists: true)
        ch_in_raw = Channel.of([ exp_meta, matrix_file])
    }

    // Create optional parameter channels
    if (params.transcript_length_matrix) { ch_transcript_lengths = Channel.of([ exp_meta, file(params.transcript_length_matrix, checkIfExists: true)]).first() } else { ch_transcript_lengths = Channel.of([[],[]]) }
    if (params.control_features) { ch_control_features = Channel.of([ exp_meta, file(params.control_features, checkIfExists: true)]).first() } else { ch_control_features = Channel.of([[],[]]) }

    if (params.gene_sets_files) {
        gene_sets_files = params.gene_sets_files.split(",")
        ch_gene_sets = Channel.of(gene_sets_files).map { file(it, checkIfExists: true) }
    } else {
        ch_gene_sets = Channel.of([[]])
    }

    // Report related files
    report_file = file(params.report_file, checkIfExists: true)
    logo_file = file(params.logo_file, checkIfExists: true)
    css_file = file(params.css_file, checkIfExists: true)
    citations_file = file(params.citations_file, checkIfExists: true)

    // ========================================================================
    // Handle contrasts
    // ========================================================================

    // Channel for the contrasts file
    if (params.contrasts_yml && params.contrasts) {
        error("Both '--contrasts' and '--contrasts_yml' parameters are set. Please specify only one of these options to define contrasts.")
    }
    if (!(params.contrasts_yml || params.contrasts)) {
        error("Either '--contrasts' and '--contrasts_yml' must be set. Please specify one of these options to define contrasts.")
    }

    // SUPPORT BOTH YAML AND CSV CONTRASTS FILE
    if (params.contrasts_yml) {
        //yaml contrasts file processing
        ch_contrasts_file = Channel.from([[exp_meta, file(params.contrasts_yml)]])
        ch_contrast_variables = ch_contrasts_file
            .map { entry ->
                def yaml_file = entry[1]
                def yaml_data = new groovy.yaml.YamlSlurper().parse(yaml_file)
                yaml_data.contrasts.collect { contrast ->
                    tuple('id': contrast.comparison[0])
                }
            }
            .flatten()
            .unique() // Uniquify to keep each contrast variable only once (in case it exists in multiple lines for blocking etc.)
    } else if (params.contrasts) {
        //csv contrasts file processing
        ch_contrasts_file = Channel.from([[exp_meta, file(params.contrasts)]])
        ch_contrast_variables = ch_contrasts_file
            .splitCsv(header:true, sep:(params.contrasts.endsWith('csv') ? ',' : '\t'))
            .map{ it.tail().first() }
            .map{
                tuple('id': it.variable)
            }
            .unique()
    }

    // ========================================================================
    // Data preprocessing
    // ========================================================================

    // If we have affy array data in the form of CEL files we'll be deriving
    // matrix and annotation from them

    if (params.study_type == 'affy_array'){

        // Uncompress the CEL files archive

        UNTAR ( ch_celfiles )

        ch_affy_input = ch_input
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

        // Fetch affy outputs and reset the meta

        ch_in_raw = AFFY_JUSTRMA_RAW.out.expression
        ch_in_norm = AFFY_JUSTRMA_NORM.out.expression

        ch_affy_platform_features = AFFY_JUSTRMA_RAW.out.annotation

        ch_versions = ch_versions
            .mix(AFFY_JUSTRMA_RAW.out.versions)

    } else if (params.study_type == 'maxquant'){

        // We'll be running Proteus once per unique contrast variable to generate plots
        // TODO: there should probably be a separate plotting module in proteus to simplify this

        // Run proteus to import protein abundances
        PROTEUS(
            ch_contrast_variables.combine(proteus_in)
        )

        // Re-map the proteus output tables to the study ID as the tables are the same across contrasts, only one norm table will be necessary
        ch_in_raw = PROTEUS.out.raw_tab
            .first()
            .map{ meta, matrix -> tuple(exp_meta, matrix) }
        ch_in_norm = PROTEUS.out.norm_tab
            .first()
            .map{ meta, matrix -> tuple(exp_meta, matrix) }

        ch_versions = ch_versions.mix(PROTEUS.out.versions)
    } else if(params.study_type == 'geo_soft_file'){

        GEOQUERY_GETGEO(ch_querygse)
        ch_in_norm = GEOQUERY_GETGEO.out.expression
        ch_soft_features = GEOQUERY_GETGEO.out.annotation

        ch_versions = ch_versions
            .mix(GEOQUERY_GETGEO.out.versions)
    }
    //// Fetch or derive a feature annotation table

    // If user has provided a feature annotation table, use that
    if (params.features){
        ch_features = Channel.of([ exp_meta, file(params.features, checkIfExists: true)])
    } else if (params.study_type == 'affy_array'){
        ch_features = ch_affy_platform_features
    } else if(params.study_type == 'geo_soft_file') {
        ch_features = ch_soft_features
    } else if (params.gtf){
        // Get feature annotations from a GTF file, gunzip if necessary

        file_gtf_in = file(params.gtf)
        file_gtf = [ [ "id": file_gtf_in.simpleName ], file_gtf_in ]

        if ( params.gtf.endsWith('.gz') ){
            GUNZIP_GTF(file_gtf)
            file_gtf = GUNZIP_GTF.out.gunzip
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        }

        // Get a features table from the GTF and combine with the matrix and sample
        // annotation (fom = features/ observations/ matrix)

        GTF_TO_TABLE( file_gtf, [[ "id":""], []])
        ch_features = GTF_TO_TABLE.out.feature_annotation
            .map{
                tuple( exp_meta, it[1])
            }

        // Record the version of the GTF -> table tool

        ch_versions = ch_versions
            .mix(GTF_TO_TABLE.out.versions)
    }
    else{

        // Otherwise we can just use the matrix input; save it to the workdir so that it does not
        // just appear wherever the user runs the pipeline
        matrix_as_anno_filename = "${workflow.workDir}/${matrix_file.getBaseName()}_as_anno.${matrix_file.getExtension()}"
        if (params.study_type == 'maxquant'){
            ch_features_matrix = ch_in_norm
        } else {
            ch_features_matrix = ch_in_raw
        }
        ch_features = ch_features_matrix
            .map{ meta, matrix ->
                matrix_copy = file(matrix_as_anno_filename)
                matrix_copy.exists() && matrix.getText().md5().equals(matrix_copy.getText().md5()) ?: matrix.copyTo(matrix_as_anno_filename)
                [ meta, file(matrix_as_anno_filename) ]
            }
    }

    // Check compatibility of FOM elements and contrasts
    if (params.study_type == 'affy_array' || params.study_type == 'maxquant'){
        ch_matrices_for_validation = ch_in_raw
            .join(ch_in_norm)
            .map{tuple(it[0], [it[1], it[2]])}
    }
    else if (params.study_type == 'geo_soft_file'){
        ch_matrices_for_validation = ch_in_norm
    }
    else{
        ch_matrices_for_validation = ch_in_raw
    }

    VALIDATOR(
        ch_input.join(ch_matrices_for_validation),
        ch_features,
        ch_contrasts_file
    )

    // For Affy, we've validated multiple input matrices for raw and norm,
    // we'll separate them out again here

    if (params.study_type == 'affy_array' || params.study_type == 'maxquant'){
        ch_validated_assays = VALIDATOR.out.assays
            .transpose()
            .branch {
                raw: it[1].name.contains('raw')
                normalised: it[1].name =~ /normali[sz]ed/
            }
        ch_raw = ch_validated_assays.raw
        ch_norm = ch_validated_assays.normalised
    }
    else if (params.study_type == 'geo_soft_file') {
        ch_norm = VALIDATOR.out.assays
    }

    if(params.study_type != 'rnaseq') {
        ch_matrix_for_differential = ch_norm
    }
    else{
        ch_raw = VALIDATOR.out.assays
        ch_matrix_for_differential = ch_raw
    }

    // Split the contrasts up so we can run differential analyses and
    // downstream plots separately.
    // Replace NA strings that might have snuck into the blocking column

    ch_contrasts = VALIDATOR.out.contrasts
        .map{it[1]}
        .splitCsv ( header:true, sep:'\t' )
        .map{
            it.blocking = it.blocking.replaceAll('^NA$', '')
            if (!it.id){
                it.id = it.values().join('_')
            }
            tuple(it, it.variable, it.reference, it.target)
        }

    // Firstly Filter the input matrix

    CUSTOM_MATRIXFILTER(
        ch_matrix_for_differential,
        VALIDATOR.out.sample_meta
    )

    // ========================================================================
    // Differential analysis
    // ========================================================================

    // Prepare inputs for differential processes

    ch_differential_input = CUSTOM_MATRIXFILTER.out.filtered
        .combine(ch_tools)
        .map { meta, matrix, tools_norm, tools_diff, tools_func ->
            [
                // as the differential subworkflow will automatically add
                // method_differential into meta, here we only need to add
                // args_differential into meta. In this way, the differential
                // outputs should always have in the meta the method_differential
                // and args_differential
                meta + [args_differential: tools_diff.args],
                matrix,
                tools_diff.method,
                tools_diff.fc_threshold,
                tools_diff.stat_threshold
            ]
        }

    // Run differential analysis

    ABUNDANCE_DIFFERENTIAL_FILTER(
        ch_differential_input,
        VALIDATOR.out.sample_meta,
        ch_transcript_lengths,
        ch_control_features,
        ch_contrasts
    )

    // collect differential results

    ch_differential_results = ABUNDANCE_DIFFERENTIAL_FILTER.out.results_genewise
    ch_differential_results_filtered = ABUNDANCE_DIFFERENTIAL_FILTER.out.results_genewise_filtered
    ch_differential_model = ABUNDANCE_DIFFERENTIAL_FILTER.out.model
    ch_differential_norm = ABUNDANCE_DIFFERENTIAL_FILTER.out.normalised_matrix
    ch_differential_varstab = ABUNDANCE_DIFFERENTIAL_FILTER.out.variance_stabilised_matrix

    ch_versions = ch_versions
        .mix(ABUNDANCE_DIFFERENTIAL_FILTER.out.versions)

    // when the study type is rnaseq, we use the normalised matrix from the differential subworkflow
    if (params.study_type == 'rnaseq') {
        ch_norm = ch_differential_norm
    } else {
        // when the study type is not rnaseq, the normalised matrix comes directly from VALIDATOR.out
        // we need to update the meta so that it can match with tools_norm
        ch_norm = ch_norm
            .map { meta, norm ->
                def meta_new = meta + [method_differential: 'validator', args_differential: [:]]
                [meta_new, norm]
            }
    }

    // prepare channel with normalized matrix, and variance stabilized matrices when available
    ch_processed_matrices = ch_norm.join(ch_differential_varstab, remainder: true)
        .map { meta, norm, vs ->
            def matrices = vs ? [norm] + vs : [norm]
            [meta, matrices]
        }

    // ========================================================================
    // Functional analysis
    // ========================================================================

    // Prepare background file - for the moment it is only needed for gprofiler2
    // TODO: prepare this channel based on ch_tools

    ch_background = Channel.of([[]])
    if (params.functional_method == 'gprofiler2') {
        if (params.gprofiler2_background_file == "auto") {
            // If auto, use input matrix as background
            ch_background = CUSTOM_MATRIXFILTER.out.filtered.map{it.tail()}.first()
        } else {
            ch_background = Channel.from(file(params.gprofiler2_background_file, checkIfExists: true))
        }
    }

    // Prepare input for functional analysis

    // Some functional analysis methods act directly on the differential analysis results, some on the normalised matrix.
    // By crossing with ch_tools, we pair the correct input file with the correct functional analysis method, taking into
    // account the upstream differential method, args and data type (normalized matrix or filtered differential analysis
    // results).
    ch_functional_input = ch_norm.map { meta, input ->
            [[method: meta.method_differential, args: meta.args_differential, type: 'norm'], meta, input]
        }
        .mix(
            ch_differential_results_filtered.map { meta, input ->
                [[method: meta.method_differential, args: meta.args_differential, type: 'filtered'], meta, input]
            }
        )
        .cross(
            ch_tools.map { tools_norm, tools_diff, tools_func ->
                [[method: tools_norm.method, args: tools_norm.args, type: tools_func.input_type], tools_func]
            }
        )
        .map { input, tools ->
            // the functional analysis subworkflow will automatically add method_functional into meta
            // so here we only need to add args_functional into meta. In this way, the functional outputs
            // will have in the meta: method_differential, args_differential, method_functional and args_functional
            def meta = input[1] + [args_functional: tools[1].args]
            [meta, input[2], tools[1].method]  // meta, input, functional analysis method
        }
        .combine(ch_gene_sets)
        .combine(ch_background)
        .map { meta, input, method, gene_sets, background ->
            [meta, input, gene_sets, background, method]
        }

    // Run functional analysis

    DIFFERENTIAL_FUNCTIONAL_ENRICHMENT(
        ch_functional_input,
        ch_contrasts,
        VALIDATOR.out.sample_meta,
        VALIDATOR.out.feature_meta.combine(Channel.of([params.features_id_col, params.features_name_col]))
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
    // Plot figures
    // ========================================================================

    // The exploratory plots are made by coloring by every unique variable used
    // to define contrasts

    ch_contrast_variables = ch_contrasts
        .map{
            [ "id": it[1] ]
        }
        .unique()

    // For geoquery we've done no matrix processing and been supplied with the
    // normalised matrix, which can be passed through to downstream analysis
    if (params.study_type == "geo_soft_file") {
        ch_mat = ch_norm.map {meta, norm -> [meta, [norm]]}
    } else {
        // otherwise, we take all the matrices: raw, normalized matrix,
        // and variance stabilized matrices if available
        ch_mat = ch_raw.map{ it[1] }
            .combine(ch_processed_matrices)
            .map { raw, meta, matrices ->
                [meta, [raw] + matrices]
            }
    }

    ch_all_matrices = VALIDATOR.out.sample_meta                 // meta_exp, samples
        .join(VALIDATOR.out.feature_meta)                       // meta_exp, samples, features
        .combine(ch_mat)                                        // meta_mat, list of matrices (raw, norm, variance stabilized)
        .map { meta_exp, samples, features, meta_mat, matrices ->
            return [meta_mat, samples, features, matrices]
        }

    // Exploratory analysis

    PLOT_EXPLORATORY(
        ch_contrast_variables
            .combine(ch_all_matrices.map{ it.tail() })
    )

    // Plot differential analysis results

    PLOT_DIFFERENTIAL(
        ch_differential_results,
        ch_all_matrices.first()
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

    // Generate a list of files that will be used by the markdown report

    ch_report_file = Channel.from(report_file)
        .map{ tuple(exp_meta, it) }

    ch_logo_file = Channel.from(logo_file)
    ch_css_file = Channel.from(css_file)
    ch_citations_file = Channel.from(citations_file)

    // create a list of input files for the report
    // This list will be created considering the tools used in the pipeline.
    // So all the files that were produced following a concrete set of tools
    // will go together to the report. This is needed later on for running
    // multiple tools in parallel (future implementation). It is also needed
    // to transfer the args in ch_tools properly to the report module.

    // first we parse the matrices channel with a tool-based key
    ch_matrices_with_key = ch_all_matrices
        .map { meta, samples, features, matrices ->
            [[method: meta.method_differential, args: meta.args_differential], samples, features, matrices]
        }
        .cross (
            ch_tools.map { tools_norm, tools_diff, tools_func -> [
                [method: tools_norm.method, args: tools_norm.args],
                [method_differential: tools_diff.method, args_differential: tools_diff.args, method_functional: tools_func.method, args_functional: tools_func.args]
            ]}
        )
        .map { files, tools ->
            [tools[1], files[1], files[2], files[3]]  // key, samples, features, [ matrices ]
        }

    // we parse the differential analysis outputs channel with a tool-based key
    ch_differential_with_key = ch_differential_results
        .join(ch_differential_model, remainder: true)
        .map { meta, results, model ->
            [[method: meta.method_differential, args: meta.args_differential], results, model]
        }
        .groupTuple()   // group all the files produced for different contrasts together
        .cross (
            ch_tools.map { tools_norm, tools_diff, tools_func -> [
                [method: tools_diff.method, args: tools_diff.args],
                [method_differential: tools_diff.method, args_differential: tools_diff.args, method_functional: tools_func.method, args_functional: tools_func.args]
            ]}
        )
        .map { files, tools ->
            [tools[1], files[1], files[2]]  // key, [ results ], [ models ]
        }

    // we parse the functional analysis outputs channel with a tool-based key
    ch_functional_with_key = ch_functional_results
        .map {[
            [method_differential: it[0].method_differential, args_differential: it[0].args_differential, method_functional: it[0].method_functional, args_functional: it[0].args_functional],
            it.tail()
        ]}
        .groupTuple()   // group all the files produced for different contrasts together
        .cross (
            ch_tools.map { tools_norm, tools_diff, tools_func -> [
                // note that the normalized matrix is used to run gsea, whereas differential results are used to run gprofiler2
                [method_differential: (tools_func.method == 'gsea' ? tools_norm.method : tools_diff.method), args_differential: (tools_func.method == 'gsea' ? tools_norm.args : tools_diff.args), method_functional: tools_func.method, args_functional: tools_func.args],
                [method_differential: tools_diff.method, args_differential: tools_diff.args, method_functional: tools_func.method, args_functional: tools_func.args]
            ]}
        )
        .map { files, tools ->
            [tools[1], files.tail()]    // key, [ files ]
        }

    // we combine the channels defined above with the contrasts, collated versions,
    // logo, css and citations files. We also combine with the outputs from differential
    // and functional analysis using the tool-based key:
    // [method_differential, args_differential, method_functional, args_functional]
    ch_report_input_files = ch_matrices_with_key
        .combine(VALIDATOR.out.contrasts.map{it.tail()})
        .combine(ch_collated_versions)
        .combine(ch_logo_file)
        .combine(ch_css_file)
        .combine(ch_citations_file)
        .join(ch_differential_with_key)
        .join(ch_functional_with_key, remainder: true)
        .map {
            // functional analysis is optional. When it is not run, joining an empty
            // ch_functional_with_key channel with remainder:true will result in a
            // channel with a meta containing null values for method_functional and
            // args_functional, and null files. Since args_functional should be a
            // map, here we update meta to keep args_functional as a map [:]. We
            // also remove the null files using grep().
            def meta = (it[0].method_functional) ? it[0] : it[0] + [args_functional: [:]]
            [meta, it.tail().grep().flatten()]    // meta, files
        }


    // Run IMMUNEDECONV
    if (params.immunedeconv_run){
        matrix_file = file(params.matrix, checkIfExists:true)
        IMMUNEDECONV(
            [ [id:matrix_file.baseName], matrix_file, params.immunedeconv_method, params.immunedeconv_function ],
            params.features_name_col
        )
        ch_versions = ch_versions
            .mix(IMMUNEDECONV.out.versions)
    }

    if (params.shinyngs_build_app){

        // Make (and optionally deploy) the shinyngs app

        // parse matrices and differential results considering the tool-based key
        ch_matrices_with_differential = ch_matrices_with_key
            .join(ch_differential_with_key)
            .multiMap { key, samples, features, matrices, results, models ->
                def meta = exp_meta + key
                matrices:
                [meta, samples, features, matrices]   // meta, samples, features, [  matrices ]
                results:
                [meta, results]                       // meta, [ results ]
            }

        // Make a new contrasts file from the differential metas to guarantee the
        // same order as the differential results
        def contrast_columns = ['variable','reference','target','blocking']
        ch_app_differential = ch_differential_results.first().map{contrast_columns.join(',')}
            .concat(
                ch_differential_results.map{it[0].subMap(contrast_columns).values().join(',')}
            )
            .collectFile(name: 'contrasts.csv', newLine: true, sort: false)
            // combine with differential results
            .combine (ch_matrices_with_differential.results)
            .map { contrasts, meta, results ->
                [meta, contrasts, results]
            }

        SHINYNGS_APP(
            ch_matrices_with_differential.matrices, // meta, samples, features, [  matrices ]
            ch_app_differential, // meta, contrasts, [differential results]
            params.exploratory_assay_names.split(',').findIndexOf { it == params.exploratory_final_assay } + 1
        )
        ch_versions = ch_versions.mix(SHINYNGS_APP.out.versions)
    }

    // Make a params list - starting with the input matrices and the relevant
    // params to use in reporting

    def report_file_names = [ 'observations', 'features' ] +
        params.exploratory_assay_names.split(',').collect { "${it}_matrix".toString() } +
        [ 'contrasts_file', 'versions_file', 'logo', 'css', 'citations' ]

    // Condition params reported on study type
    def params_pattern = "report|gene_sets|study|observations|features|filtering|exploratory|differential"
    if (params.study_type == 'affy_array' || params.study_type == 'geo_soft_file'){
        params_pattern += "|affy"
    } else if (params.study_type == 'maxquant'){
        params_pattern += "|proteus"
    }

    // prepare the final report params
    ch_report_params = ch_report_input_files
        .map{ meta, files ->
            // the params scope don't necessarily include the args parsed from toolsheet.
            // Hence we need to consider both params and the args parsed from toolsheet
            // (args_differential and args_functional) to create the final report params.
            def params_with_tools = params + meta.args_differential + meta.args_functional

            // update the pattern with tools-specific info
            def pattern_with_tools = params_pattern +
                "|${meta.method_differential}" +
                (meta.method_functional ? "|functional|${meta.method_functional}" : "")

            // Extract and return parameters for the report by creating a map where
            // the keys are parameter names matching the pattern and the values are
            // their corresponding parameter values. And add this map to the report
            // files map.
            params_with_tools.findAll{ k,v -> k.matches(~/(${pattern_with_tools}).*/) } +
            [report_file_names, files.collect{ f -> f.name}].transpose().collectEntries()
        }

    // Render the final report

    RMARKDOWNNOTEBOOK(
        ch_report_file,
        ch_report_params,
        ch_report_input_files.map{it.tail().flatten()}
    )

    // Make a report bundle comprising the markdown document and all necessary
    // input files

    MAKE_REPORT_BUNDLE(
        RMARKDOWNNOTEBOOK.out.parameterised_notebook
            .combine(ch_report_input_files.map{it.tail().flatten()})
            .map{[it[0], it[1..-1]]}
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
