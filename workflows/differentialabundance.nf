/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HANDLE INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
def exp_meta = [ "id": params.study_name  ]
if (params.input) { ch_input = Channel.of([ exp_meta, file(params.input, checkIfExists: true) ]) } else { exit 1, 'Input samplesheet not specified!' }

if (params.study_type == 'affy_array') {
    if (params.affy_cel_files_archive) {
        ch_celfiles = Channel.of([ exp_meta, file(params.affy_cel_files_archive, checkIfExists: true) ])
    } else {
        error("CEL files archive not specified!")
    }
} else if (params.study_type == 'maxquant') {

    // Should the user have enabled --gsea_run, throw an error
    if (params.gsea_run) {
        error("Cannot run GSEA for maxquant data; please set --gsea_run to false.")
    }
    if (params.gprofiler2_run){
        error("gprofiler2 pathway analysis is not yet possible with maxquant input data; please set --gprofiler2_run false and rerun pipeline!")
    }
    if (!params.matrix) {
        error("Input matrix not specified!")
    }
    matrix_file = file(params.matrix, checkIfExists: true)

    // Make channel for proteus
    proteus_in = Channel.of([ file(params.input), matrix_file ])
} else if (params.study_type == 'geo_soft_file') {

    // To pull SOFT files from a GEO a GSE study identifer must be provided

    if (params.querygse && params.features_metadata_cols) {
        ch_querygse = Channel.of([exp_meta, params.querygse])
    } else {
        error("Query GSE not specified or features metadata columns not specified")
    }
} else {
    // If this is not microarray data or maxquant output, and this an RNA-seq dataset,
    // then assume we're reading from a matrix

    if (params.study_type == "rnaseq" && params.matrix) {
        matrix_file = file(params.matrix, checkIfExists: true)
        ch_in_raw = Channel.of([ exp_meta, matrix_file])
    } else {
        error("Input matrix not specified!")
    }

}

// Check optional parameters
if (params.transcript_length_matrix) { ch_transcript_lengths = Channel.of([ exp_meta, file(params.transcript_length_matrix, checkIfExists: true)]).first() } else { ch_transcript_lengths = Channel.of([[],[]]) }
if (params.control_features) { ch_control_features = Channel.of([ exp_meta, file(params.control_features, checkIfExists: true)]).first() } else { ch_control_features = Channel.of([[],[]]) }

def run_gene_set_analysis = params.gsea_run || params.gprofiler2_run

ch_gene_sets = Channel.of([[]])
if (run_gene_set_analysis) {
    if (params.gene_sets_files) {
        gene_sets_files = params.gene_sets_files.split(",")
        ch_gene_sets = Channel.of(gene_sets_files).map { file(it, checkIfExists: true) }
        if (params.gprofiler2_run && (!params.gprofiler2_token && !params.gprofiler2_organism) && gene_sets_files.size() > 1) {
            error("gprofiler2 can currently only work with a single gene set file")
        }
    } else if (params.gsea_run) {
        error("GSEA activated but gene set file not specified!")
    } else if (params.gprofiler2_run) {
        if (!params.gprofiler2_token && !params.gprofiler2_organism) {
            error("To run gprofiler2, please provide a run token, GMT file or organism!")
        }
    }
}

// Define tool channel
// For the moment, it allows a single combination of tools, but this will be changed with
// the toolsheet and introduction of multi-tool logic.
if (params.study_type == 'affy_array' ||
    params.study_type == 'geo_soft_file' ||
    params.study_type == 'maxquant' ||
    (params.study_type == 'rnaseq' && params.differential_use_limma)
    ) {
        tools_differential = ['method': 'limma']
    } else {
        tools_differential = ['method': 'deseq2']
    }
tools_differential += [
    'fc_threshold': params.differential_min_fold_change,
    'stat_threshold': params.differential_max_qval
]
tools_functional = []
if (params.gsea_run) {
    tools_functional = ['method': 'gsea', 'input_type': 'norm']
}
if (params.gprofiler2_run) {
    tools_functional = ['method': 'gprofiler2', 'input_type': 'filtered']
}
ch_tools = Channel.of([tools_differential, tools_functional])

// report related files
report_file = file(params.report_file, checkIfExists: true)
logo_file = file(params.logo_file, checkIfExists: true)
css_file = file(params.css_file, checkIfExists: true)
citations_file = file(params.citations_file, checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


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
    DEFINE FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def replace_meta_method_id (channel, new_id) {
    return channel.map { it ->
        def meta = it[0] - ['method': it[0].method]
        meta[new_id] = it[0].method
        [ meta, it[1..it.size()-1] ].flatten()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIFFERENTIALABUNDANCE {

    main:

    ch_versions = Channel.empty()
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
        .map { meta, matrix, tools_differential, tools_functional ->
            [
                meta,
                matrix,
                tools_differential.method,
                tools_differential.fc_threshold,
                tools_differential.stat_threshold
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
    // meta.method would be replaced by meta.method_differential

    ch_differential_results = replace_meta_method_id(ABUNDANCE_DIFFERENTIAL_FILTER.out.results_genewise, 'method_differential')
    ch_differential_results_filtered = replace_meta_method_id(ABUNDANCE_DIFFERENTIAL_FILTER.out.results_genewise_filtered, 'method_differential')
    ch_differential_norm = replace_meta_method_id(ABUNDANCE_DIFFERENTIAL_FILTER.out.normalised_matrix, 'method_differential')
    ch_differential_model = replace_meta_method_id(ABUNDANCE_DIFFERENTIAL_FILTER.out.model, 'method_differential')
    rlog_counts = replace_meta_method_id(ABUNDANCE_DIFFERENTIAL_FILTER.out.rlog_counts.filter{ it != null}, 'method_differential')
    vst_counts = replace_meta_method_id(ABUNDANCE_DIFFERENTIAL_FILTER.out.vst_counts.filter{ it != null}, 'method_differential')

    ch_versions = ch_versions
        .mix(ABUNDANCE_DIFFERENTIAL_FILTER.out.versions)

    if (params.study_type == 'rnaseq') ch_norm = ch_differential_norm

    // ========================================================================
    // Functional analysis
    // ========================================================================

    // Prepare background file - for the moment it is only needed for gprofiler2

    ch_background = Channel.of([[]])
    if (params.gprofiler2_run) {
        if (params.gprofiler2_background_file == "auto") {
            // If auto, use input matrix as background
            ch_background = CUSTOM_MATRIXFILTER.out.filtered.map{it.tail()}.first()
        } else {
            ch_background = Channel.from(file(params.gprofiler2_background_file, checkIfExists: true))
        }
    }

    // Prepare input for functional analysis

    ch_functional_input = ch_differential_results_filtered.combine(ch_tools.filter{it[1].input_type == 'filtered'})
        .mix(ch_norm.combine(ch_tools.filter{it[1].input_type == 'norm'}))
        .combine(ch_gene_sets)
        .combine(ch_background)
        .map { meta, input, tools_differential, tools_functional, gene_sets, background ->
            if (!('method_differential' in meta) || (meta.method_differential == tools_differential.method)) {
                return [meta, input, gene_sets, background, tools_functional.method]
            }
        }

    // Run functional analysis

    DIFFERENTIAL_FUNCTIONAL_ENRICHMENT(
        ch_functional_input,
        ch_contrasts,
        VALIDATOR.out.sample_meta,
        VALIDATOR.out.feature_meta.combine(Channel.of([params.features_id_col, params.features_name_col]))
    )

    // Prepare the results for downstream analysis

    ch_gsea_results = replace_meta_method_id(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gsea_report, 'method_functional')
    gprofiler2_plot_html = replace_meta_method_id(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gprofiler2_plot_html, 'method_functional')
    gprofiler2_all_enrich = replace_meta_method_id(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gprofiler2_all_enrich, 'method_functional')
    gprofiler2_sub_enrich = replace_meta_method_id(DIFFERENTIAL_FUNCTIONAL_ENRICHMENT.out.gprofiler2_sub_enrich, 'method_functional')

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
        ch_mat = ch_norm.combine( Channel.of([[],[],[]]) )
    } else {
        ch_processed_matrices = ch_norm
            .combine( rlog_counts.ifEmpty([[],[]]) )
            .combine( vst_counts.ifEmpty([[],[]]) )
            .map { meta, norm, meta_rlog, rlog, meta_vst, vst ->
                def matrices = [meta, norm]
                if (meta_rlog == [] || meta == meta_rlog) matrices += [rlog]
                if (meta_vst == [] || meta == meta_vst) matrices += [vst]
                return matrices   // meta, norm, rlog, vst
            }
        ch_mat = ch_raw.combine(ch_processed_matrices)
            .map { meta_exp, raw, meta_norm, norm, rlog, vst ->
                if (meta_norm.subMap(meta_exp.keySet()) == meta_exp) {
                    return [meta_norm, raw, norm, rlog, vst]
                }
            }
    }

    ch_all_matrices = VALIDATOR.out.sample_meta                 // meta, samples
        .join(VALIDATOR.out.feature_meta)                       // meta, samples, features
        .combine(ch_mat)
        .map { meta_exp, samples, features, meta_mat, raw, norm, rlog, vst ->
            if (meta_mat.subMap(meta_exp.keySet()) == meta_exp) {
                def processed_matrices = [raw]
                if (norm != []) processed_matrices += [norm]
                if (rlog != []) processed_matrices += [rlog]
                if (vst != []) processed_matrices += [vst]
                return [meta_mat, samples, features, processed_matrices]
            }
        }

    // Exploratory analysis

    PLOT_EXPLORATORY(
        ch_contrast_variables
            .combine(ch_all_matrices.map{ it.tail() })
    )

    // Differential analysis

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

    ch_collated_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'collated_versions.yml', sort: true, newLine: true)

    // Generate a list of files that will be used by the markdown report

    ch_report_file = Channel.from(report_file)
        .map{ tuple(exp_meta, it) }

    ch_logo_file = Channel.from(logo_file)
    ch_css_file = Channel.from(css_file)
    ch_citations_file = Channel.from(citations_file)

    ch_report_input_files = ch_all_matrices
        .map{ it.tail() }
        .map{it.flatten()}
        .combine(VALIDATOR.out.contrasts.map{it.tail()})
        .combine(ch_collated_versions)
        .combine(ch_logo_file)
        .combine(ch_css_file)
        .combine(ch_citations_file)
        .combine(ch_differential_results.map{it[1]}.toList())
        .combine(ch_differential_model.map{it[1]}.toList())

    if (params.gsea_run){
        ch_report_input_files = ch_report_input_files
            .combine(ch_gsea_results
                .map{it.tail()}.flatMap().toList()
            )
    }

    if (params.gprofiler2_run){
        ch_report_input_files = ch_report_input_files
            .combine(gprofiler2_plot_html.map{it[1]}.flatMap().toList())
            .combine(gprofiler2_all_enrich.map{it[1]}.flatMap().toList())
            .combine(gprofiler2_sub_enrich.map{it[1]}.flatMap().toList())
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

        // Make a new contrasts file from the differential metas to guarantee the
        // same order as the differential results

        ch_app_differential = ch_differential_results.first().map{it[0].keySet().tail().join(',')}
            .concat(
                ch_differential_results.map{it[0].values().tail().join(',')}
            )
            .collectFile(name: 'contrasts.csv', newLine: true, sort: false)
            .map{
                tuple(exp_meta, it)
            }
            .combine(ch_differential_results.map{it[1]}.collect().map{[it]})

        SHINYNGS_APP(
            ch_all_matrices,     // meta, samples, features, [  matrices ]
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
    if (params.study_type == 'rnaseq'){
        if (params.differential_use_limma){
            params_pattern += "|limma"
        } else {
            params_pattern += "|deseq2"
        }
    }
    if (params.study_type == 'affy_array' || params.study_type == 'geo_soft_file'){
        params_pattern += "|affy|limma"
    }
    if (params.study_type == 'maxquant'){
        params_pattern += "|proteus|limma"
    }
    if (params.gprofiler2_run){
        params_pattern += "|gprofiler2"
    }
    if (params.gsea_run){
        params_pattern += "|gsea"
    }
    params_pattern = ~/(${params_pattern}).*/

    ch_report_params = ch_report_input_files
        .map{
            params.findAll{ k,v -> k.matches(params_pattern) } +
            [report_file_names, it.collect{ f -> f.name}].transpose().collectEntries()
        }

    // Render the final report
    RMARKDOWNNOTEBOOK(
        ch_report_file,
        ch_report_params,
        ch_report_input_files
    )

    // Make a report bundle comprising the markdown document and all necessary
    // input files

    MAKE_REPORT_BUNDLE(
        RMARKDOWNNOTEBOOK.out.parameterised_notebook
            .combine(ch_report_input_files)
            .map{[it[0], it[1..-1]]}
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
