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

    //
    // Define tool settings based on analysis name or default parameters
    //
    ch_tools = Channel.of(getToolConfigurations())

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
* parseStrArgsToMap("--arg1 aa --arg2 bb --arg4 cc") => [arg1: aa, arg2: bb, arg4: cc]
*/
def parseStrArgsToMap(String argsStr) {
    if (!argsStr) return [:]                     // if null return empty
    def tokens = argsStr.split().findAll { it }  // Split and remove empty strings
    def pairs = tokens.collate(2)                // Group into pairs
    return pairs.collectEntries {
        [(it[0].replaceAll('^-+', '')): it[1]]   // Remove leading dashes
    }
}

/**
* Parse the parameters for a specific analysis type and method.
* @param analysisType The analysis type. Eg. 'differential', 'functional'.
* @param method The method to use for the analysis.
* @param argsStr The string of arguments to parse.
* @return A map of parameters.
*/
def parseParams(String analysisType, String method, String argsStr){
    def parsed_params = [:]

    // parse the params specific to the analysis type from the pipeline params
    parsed_params += params.findAll { k, v -> k.matches(~/(${analysisType}).*/) }

    // parse the params specific to the method
    // First it gets the default values as specified by the pipeline params scope,
    // and then overwrite with the values specified in the argsStr (highest priority),
    // when provided
    if (method) {
        parsed_params += params.findAll { k, v -> k.matches(~/(${method}).*/) }
        parsed_params["${analysisType}_method"] = method
    }
    if (argsStr) {
        parsed_params += parseStrArgsToMap(argsStr)
    }

    return parsed_params
}

/**
* Get tool configurations based on whether analysis_name is provided
* @return A channel with tool configurations.
*/
def getToolConfigurations() {
    if (params.analysis_name) {
        return getToolsheetConfigurations()
    } else {
        return getDefaultConfigurations()
    }
}

/**
* Get configurations from toolsheet
* @return A list of tool configurations by parsing toolsheet.
*/
def getToolsheetConfigurations() {
    // Load appropriate toolsheet based on study type
    def toolsheet_path = getToolsheetPath()

    // Load and filter toolsheet for specific analysis
    def toolsheet = samplesheetToList(toolsheet_path, "${projectDir}/assets/schema_tools.json")
        .find { it[0].analysis_name == params.analysis_name }
        ?: error("Analysis '${params.analysis_name}' not found in toolsheet")

    // Construct tool configurations
    def tool_config = toolsheet[0]
    return [
        getNormalizationConfig(tool_config),
        getDifferentialConfig(tool_config),
        getFunctionalConfig(tool_config)
    ]
}

/**
* Get default configurations from pipeline parameters
* @return A list of default tool configurations.
*/
def getDefaultConfigurations() {
    return [
        getNormalizationConfig(null),
        getDifferentialConfig(null),
        getFunctionalConfig(null)
    ]
}

/**
* Get the appropriate toolsheet path.
* @return The path to the toolsheet.
*/
def getToolsheetPath() {
    // use user-provided toolsheet if available
    if (params.toolsheet) {
        return params.toolsheet
    }
    // use default toolsheet based on study type
    switch (params.study_type) {
        case 'rnaseq':
            return "${projectDir}/assets/toolsheet_rnaseq.csv"
        case ['affy_array', 'geo_soft_file']:
            return "${projectDir}/assets/toolsheet_affy.csv"
        case 'maxquant':
            return "${projectDir}/assets/toolsheet_maxquant.csv"
        default:
            error("Invalid study_type. Must be one of: 'rnaseq', 'affy_array', 'geo_soft_file', 'maxquant'")
    }
}

/**
* Get normalization tool configuration
* @param tool_config The tool configuration.
* @return A map of normalization configuration.
*/
def getNormalizationConfig(tool_config) {
    // non-rnaseq studies use directly the normalized matrix output from the VALIDATOR module
    if (params.study_type != 'rnaseq') {
        return [method: 'validator', args: [:]]
    }
    // for rna-seq studies, use the differential method specified in the toolsheet
    // or the default differential method specified in the pipeline parameters
    def method = tool_config?.diff_method ?: params.differential_method
    def args = parseParams('differential', method, tool_config?.diff_args ?: null)
    return [method: method, args: args]
}

/**
* Get differential analysis tool configuration
* @param tool_config The tool configuration.
* @return A map of differential analysis configuration.
*/
def getDifferentialConfig(tool_config) {
    def method = tool_config?.diff_method ?: params.differential_method
    def args = parseParams('differential', method, tool_config?.diff_args ?: null)
    def fc_threshold = args.differential_min_fold_change
    def stat_threshold = args.differential_max_qval
    return [
        method: method,
        args: args,
        fc_threshold: fc_threshold,
        stat_threshold: stat_threshold
    ]
}

/**
* Get functional analysis tool configuration
* @param tool_config The tool configuration.
* @return A map of functional analysis configuration.
*/
def getFunctionalConfig(tool_config) {
    def method = tool_config?.func_method ?: params.functional_method
    def args = parseParams('functional', method, tool_config?.func_args ?: null)
    def input_type = method == 'gprofiler2' ? 'filtered' :
                    method == 'gsea' ? 'norm' :
                    null
    return [
        method: method,
        args: args,
        input_type: input_type
    ]
}
