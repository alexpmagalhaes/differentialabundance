//
// Subworkflow with functionality specific to the nf-core/differentialabundance pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    // ===========================================================================
    // Handle toolsheet
    // ===========================================================================

    // Define tool settings
    // Use the toolsheet information if an analysis name is provided
    // otherwise ensamble the tools channel from the command line parameters

    // TODO: for the moment we only run one analysis at a time, but in the future
    // we would enable benchmark mode to run multiple analyses.

    // TODO: define proper checks
    //   - check if analysis_name is in toolsheet
    //   - replace the checks depending on params.differential_method, etc.
    // TODO: change report rmd to fit the new params (eg. functional_method)
    // TODO: update usage.md
    // TODO: remove method-specific fixed params (eg. differential_file_suffix,
    // differential_fc_column, etc.) from user params scope. Instead, define
    // them here based on the chosen tool.

    if (params.analysis_name) {

        // use the corresponding toolsheet given the study type
        if (params.toolsheet_custom) {
            ch_toolsheet = Channel.fromList(samplesheetToList(params.toolsheet_custom, './assets/schema_tools.json'))
        } else if (params.study_type == 'rnaseq') {
            ch_toolsheet = Channel.fromList(samplesheetToList(params.toolsheet_rnaseq, './assets/schema_tools.json'))
        } else if (params.study_type == 'affy_array') {
            ch_toolsheet = Channel.fromList(samplesheetToList(params.toolsheet_affy, './assets/schema_tools.json'))
        } else if (params.study_type == 'geo_soft_file') {
            ch_toolsheet = Channel.fromList(samplesheetToList(params.toolsheet_soft, './assets/schema_tools.json'))
        } else if (params.study_type == 'maxquant') {
            ch_toolsheet = Channel.fromList(samplesheetToList(params.toolsheet_maxquant, './assets/schema_tools.json'))
        } else {
            error("Please make sure to mention the correct study_type")
        }

        // create a channel with the toolsheet meta information
        // Note that here we create a meta with the tool method info, and also
        // tool-specific arguments. We use parseArgs to parse the args from
        // toolsheet, and getParams to get the original args from the pipeline
        // params scope. In this way, when a parameter is defined in toolsheet,
        // it will have the highest priority, and overwrite the original params.
        // These args can be accessed by the module through modules.config.
        ch_tools_meta = ch_toolsheet
            .filter{ it[0].analysis_name == params.analysis_name }
            .map { it ->
                def meta = [
                    analysis_name: it[0].analysis_name,
                    diff_method  : it[0].diff_method,
                    diff_args    : (it[0].diff_args) ?
                        getParams('differential', it[0].diff_method) + parseArgs(it[0].diff_args) + [differential_method: it[0].diff_method] :
                        getParams('differential', it[0].diff_method) + [differential_method: it[0].diff_method],
                    func_method  : it[0].func_method,
                    func_args    : (it[0].func_args) ?
                        getParams('functional', it[0].func_method) + parseArgs(it[0].func_args) + [functional_method: it[0].func_method] :
                        getParams('functional', it[0].func_method) + [functional_method: it[0].func_method]
                ]
                return [meta]
            }
    } else {
        // create a channel with a meta tool-specific info based on the pipeline
        // params scope.
        ch_tools_meta = Channel.of([[
            diff_method: params.differential_method,
            diff_args  : getParams('differential', params.differential_method),
            func_method: params.functional_method,
            func_args  : getParams('functional', params.functional_method)
        ]])
    }

    // parse channel tools into a proper structure for downstream processes
    ch_tools = ch_tools_meta
        .map{ it ->
            // Normalization tool:
            // in rnaseq, we use the data normalized by the differential tool
            // in non-rnaseq studies, like for example array, we use the data
            // normalized from the VALIDATOR module
            def tools_normalization = (params.study_type == 'rnaseq') ?
                [method: it[0].diff_method, args: it[0].diff_args] :
                [method: 'validator', args: [:]]
            // Differential analysis tool:
            // Also set fold change and q-value thresholds, required as input for
            // the differential abundance analysis subworkflow
            def tools_differential = [
                method        : it[0].diff_method,
                args          : it[0].diff_args,
                fc_threshold  : ('differential_min_fold_change' in it[0].diff_args) ?
                    it[0].diff_args.differential_min_fold_change :
                    params.differential_min_fold_change,
                stat_threshold: ('differential_max_qval' in it[0].diff_args) ?
                    it[0].diff_args.differential_max_qval :
                    params.differential_max_qval
            ]
            // Functional analysis tool:
            // Use gprofiler2 (filtered input) if enabled, else gsea (normalized input) if enabled.
            def tools_functional = [
                method    : it[0].func_method,
                args      : it[0].func_args,
                input_type: it[0].func_method == 'gprofiler2' ? 'filtered' : (it[0].func_method == 'gsea' ? 'norm' : null)
            ]
            return [ tools_normalization, tools_differential, tools_functional ]
        }

        ch_tools.view{"ch_tools: $it"}

    emit:
    tools       = ch_tools
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications


    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",

            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [

        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

/**
* Parse a string of arguments into a map.
* @param argsStr The string of arguments to parse.
* @return A map of parameters.
* @example
* parseArgs("--arg1 aa --arg2 bb --arg4 cc") => [arg1: aa, arg2: bb, arg4: cc]
*/
def parseArgs(String argsStr) {
    println("argsStr: $argsStr")
    if (!argsStr) return [:]
    def tokens = argsStr.split().findAll { it }  // Split and remove empty strings
    def pairs = tokens.collate(2)                // Group into pairs
    return pairs.collectEntries {
        [(it[0].replaceAll('^-+', '')): it[1]]   // Remove leading dashes
    }
}

/**
* Get params from the pipeline params scope.
* @param paramsType The base pattern to match the params against.
* @param method The method to get the params for.
* @return A map of params.
*/
def getParams(String basePattern, String method) {
    print("basePattern: $basePattern, method: $method")
    if (!method) return [:]
    pattern = "$basePattern|$method"
    return params.findAll { k, v -> k.matches(~/(${pattern}).*/) }
}
