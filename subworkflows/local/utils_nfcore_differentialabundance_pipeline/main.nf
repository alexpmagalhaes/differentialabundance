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

    // Check that params is available
    if (!params) {
        error("Pipeline parameters not initialized. This is a critical error.")
    }

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
    validate_params = (params.validate_params && !params.analysis_name && !params.toolsheet_custom) ? true : false  // the validation will be done through toolsheet, when provided
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

    // Define tool settings based on analysis name or default parameters
    ch_paramsets = Channel.fromList(getToolConfigurations())

    emit:
    paramsets = ch_paramsets
    versions  = ch_versions
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
    SUBWORKFLOW FOR GROUPING BY ANALYSIS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Toolsheet usage means that we sometimes group multiple analyses into a single run
where the input parameterisation would have been the same. There may also be multiple
results per group where we have multiple contrasts. To use these results downstream
we must 'expand' results back to the original groupings, but also group by the
individual analysis such that e.g. the results from multiple contrasts travel
through the workflow together.
*/

workflow GROUP_BY_ANALYSIS {

    take:
    channel_in
    paramsets_by_name

    main:

    channel_out = channel_in
        .map{meta, differential_results -> [meta.analysis_names, meta, differential_results]}
        .transpose()
        .join(paramsets_by_name)
        .map{analysis_name, meta, differential_results, paramset -> [paramset, meta, differential_results]}


    emit:
    with_study_meta = channel_out
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

    // Check file existence for input parameters
    def checkPathParamList = [ params.input ]
    for (param in checkPathParamList) {
        if (param) {
            file(param, checkIfExists: true)
        }
    }

    // Check mandatory parameters
    if (!params.input) {
        error("Input samplesheet not specified!")
    }

    // Validate study type specific parameters
    if (params.study_type == 'affy_array') {
        if (!params.affy_cel_files_archive) {
            error("CEL files archive not specified!")
        }
    } else if (params.study_type == 'maxquant') {
        if (params.functional_method) {
            error("Functional analysis is not yet possible with maxquant input data; please set --functional_method to null and rerun pipeline!")
        }
        if (!params.matrix) {
            error("Input matrix not specified!")
        }
    } else if (params.study_type == 'geo_soft_file') {
        if (!params.querygse || !params.features_metadata_cols) {
            error("Query GSE not specified or features metadata columns not specified")
        }
    } else if (params.study_type == "rnaseq") {
        if (!params.matrix) {
            error("Input matrix not specified!")
        }
    }

    // Validate functional analysis parameters
    if (params.functional_method) {
        if (params.functional_method == 'gsea' && !params.gene_sets_files) {
            error("GSEA activated but gene set file not specified!")
        } else if (params.functional_method == 'gprofiler2') {
            if (!params.gprofiler2_token && !params.gprofiler2_organism) {
                error("To run gprofiler2, please provide a run token, GMT file or organism!")
            }
            if (params.gene_sets_files && params.gene_sets_files.split(",").size() > 1) {
                error("gprofiler2 can currently only work with a single gene set file")
            }
        }
    }

    // Validate contrasts parameters
    if (params.contrasts_yml && params.contrasts) {
        error("Both '--contrasts' and '--contrasts_yml' parameters are set. Please specify only one of these options to define contrasts.")
    }
    if (!(params.contrasts_yml || params.contrasts)) {
        error("Either '--contrasts' and '--contrasts_yml' must be set. Please specify one of these options to define contrasts.")
    }
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

// Get tool configurations based on whether analysis_name is provided
def getToolConfigurations() {
    // Use toolsheet if either analysis_name or custom toolsheet is provided, otherwise use default params
    return (params.analysis_name || params.toolsheet_custom) ?
        getToolsheetConfigurations() :
        getDefaultConfigurations()
}

// Get configurations from toolsheet
def getToolsheetConfigurations() {
    // get toolsheet path
    def toolsheet_path = params.toolsheet_custom ?: "${projectDir}/assets/toolsheet.csv"
    // get toolsheet configurations - and validate
    def toolsheet = getToolsheetConfigurationsOnly(toolsheet_path)
    // get toolsheet configurations merged with pipeline params - and validate
    def toolsheet_full = getToolsheetConfigurationsFull(toolsheet)

    return toolsheet_full
}

// parse the params from toolsheet, and validate each of them
def getToolsheetConfigurationsOnly(toolsheet_path) {
    // Create temporary schema for validation - with only the fields in the toolsheet
    def schema_path = createToolsheetSchema(toolsheet_path)

    // Load toolsheet and validate each row against the transformed schema
    def raw_toolsheet = samplesheetToList(toolsheet_path, schema_path)

    def toolsheet = raw_toolsheet
        .findAll { row ->
            // If analysis_name is set, only keep matching rows
            if (params.analysis_name) {
                return row[0].analysis_name == params.analysis_name
            }
            return true
        }

    if (toolsheet.isEmpty()) {
        if (params.analysis_name) {
            error("No configuration found in toolsheet for analysis_name '${params.analysis_name}'")
        } else {
            error("No valid configurations found in toolsheet")
        }
    }

    return toolsheet
}

// merge the toolsheet params with pipeline params and validate them
def getToolsheetConfigurationsFull(toolsheet) {
    // gather toolsheet params with pipeline params - to create full paramsets
    paramsets = toolsheet
        .collect { row ->
            def ignore = ['help', 'help_full', 'show_hidden', 'genomes']  // these are not in the schema
            def paramsets = params + row[0]  // merge pipeline params with toolsheet params
            return paramsets.findAll { k,v ->
                // Only include parameters that are not in the ignore list
                !ignore.contains(k)
            }
        }

    // write the full paramsets to a temporary json file
    def tempFile = new File("${projectDir}/extended.json")
    tempFile.text = new groovy.json.JsonBuilder(paramsets).toPrettyString()
    def full_toolsheet_path = tempFile.absolutePath

    // create temporary toolsheet schema with all the parameters
    def full_schema_path = createToolsheetSchema(null)

    // validate the full paramsets
    def paramsets_validated = samplesheetToList(full_toolsheet_path, full_schema_path)

    // use only the static params, otherwise it causes issues with cache
    def staticparams = params.findAll { k, v -> k != 'trace_report_suffix' } as Map
    return paramsets_validated
        .collect{ row ->
            return ['id': staticparams.study_name] + staticparams + row[0]
        }
}

// Get default configurations from pipeline parameters
def getDefaultConfigurations() {
    return [['id': params.study_name] + params]
}

// Create a temporary schema file for toolsheet validation
// This schema is derived from the pipeline's nextflow_schema.json,
// to create a schema with toolsheet structure that can be used to
// run samplesheetToList. When the toolsheet is provided, it will
// also filter the properties to only include those present in the
// toolsheet. Concrete steps are:
// 1. Extracting properties from the pipeline schema
// 2. Determining which properties and required fields to keep, if
// toolsheet provided
// 3. Creating a schema that allows for an array of objects
// 4. Adding meta fields to all properties
// @return The absolute path to the temporary schema file
def createToolsheetSchema(toolsheet_path) {
    // Load and parse pipeline schema
    def pipeline_schema = new File("${projectDir}/nextflow_schema.json").text
    def schema_json = new groovy.json.JsonSlurper().parseText(pipeline_schema)

    // Determine which properties and required fields to keep, if toolsheet is provided
    def keep_properties = null
    def keep_required = null
    if (toolsheet_path) {
        def toolsheet_lines = new File(toolsheet_path).readLines()
        def headers = toolsheet_lines[0].split(',')
        // Ensure analysis_name is in headers
        if (!headers.contains('analysis_name')) {
            error("The toolsheet must contain an 'analysis_name' column")
        }
        keep_properties = headers
        keep_required = ['analysis_name']
    }

    // Extract properties and required fields from schema
    def all_properties = extractPropertiesFromSchema(schema_json, keep_properties)
    def all_required = extractRequiredFromSchema(schema_json, keep_required)

    // Create samplesheet schema with filtered properties and analysis_name as required
    def samplesheet_schema = [
        '$schema': 'https://json-schema.org/draft/2020-12/schema',
        'title': 'nf-core/differentialabundance - toolsheet schema',
        'description': 'Schema for validating the toolsheet configuration',
        'type': 'array',
        'items': [
            'type': 'object',
            'properties': all_properties,
            'required': ['analysis_name']
        ]
    ]

    // Write temporary schema file
    def temp_schema = File.createTempFile("samplesheet_schema", ".json")
    def schema_string = new groovy.json.JsonBuilder(samplesheet_schema).toPrettyString()
    temp_schema.text = schema_string

    return temp_schema.absolutePath
}

// extract properties from schema
// if keep is provided, keep only those properties
def extractPropertiesFromSchema(schema, keep) {
    def all_properties = [:]

    // Extract properties from schema
    schema.$defs.each { group_name, group_def ->
        if (group_def.properties) {
            group_def.properties.each { prop_name, prop_def ->
                // Only include properties that are in the toolsheet headers,
                // if toolsheet is provided
                if (!keep || keep.contains(prop_name)) {
                    // Add meta field to each property
                    def prop_with_meta = prop_def + [meta: [prop_name]]
                    all_properties[prop_name] = prop_with_meta
                }
            }
        }
    }

    // Ensure analysis_name property exists in schema
    if (!all_properties.containsKey('analysis_name')) {
        // Find analysis_name definition in pipeline schema
        def analysis_name_def = null
        schema.$defs.each { group_name, group_def ->
            if (group_def.properties && group_def.properties.containsKey('analysis_name')) {
                analysis_name_def = group_def.properties['analysis_name']
            }
        }

        if (!analysis_name_def) {
            error("Could not find analysis_name definition in pipeline schema")
        }

        // Add analysis_name to properties with meta field
        all_properties['analysis_name'] = analysis_name_def + [meta: ['analysis_name']]
    }

    return all_properties
}

// extract required fields from schema
def extractRequiredFromSchema(schema, keep) {
    def all_required = []
    schema.$defs.each { group_name, group_def ->
        if (group_def.required) {
            group_def.required.each { req_prop ->
                if (!keep || keep.contains(req_prop)) {
                    all_required << req_prop
                }
            }
        }
    }
    return all_required
}

// Split parameters into different sets based on prefixes and analysis type
def splitParameters(params) {
    // Base parameters are everything not in other sets
    def baseParams = [:]
    def differentialParams = [:]
    def functionalParams = [:]

    // First collect all differential parameters
    params.each { key, value ->
        if (key.startsWith('differential_') || key.startsWith('limma_') || key.startsWith('deseq2_')) {
            differentialParams[key] = value
        }
    }

    // Then collect functional parameters based on method
    if (params.functional_method == 'gprofiler2') {
        // For gprofiler2, include both base and differential params
        params.each { key, value ->
            if (key.startsWith('functional_') || key.startsWith('gprofiler2_')) {
                functionalParams[key] = value
            }
        }
        // Add all differential params to functional params for gprofiler2
        functionalParams.putAll(differentialParams)
    } else if (params.functional_method == 'gsea') {
        // For GSEA, only include functional params
        params.each { key, value ->
            if (key.startsWith('functional_') || key.startsWith('gsea_')) {
                functionalParams[key] = value
            }
        }
    }

    // Base params are everything not in other sets
    params.each { key, value ->
        if (!differentialParams.containsKey(key) && !functionalParams.containsKey(key)) {
            baseParams[key] = value
        }
    }

    // Add base params to both differential and functional
    differentialParams.putAll(baseParams)
    functionalParams.putAll(baseParams)

    return [
        base: baseParams,
        differential: differentialParams,
        functional: functionalParams
    ]
}

// Simplify meta map to only include parameters relevant for a specific analysis type
def subsetMeta(meta, analysis_type) {
    def relevantParams = [:]

    // Define static lists of parameter prefixes for different analysis types
    def differentialPrefixes = ['differential_', 'limma_', 'deseq2_']
    def functionalPrefixes = ['functional_']
    def gseaPrefixes = ['gsea_']
    def gprofiler2Prefixes = ['gprofiler2_']
    def basicFields = ['id', 'study_name', 'study_type', 'analysis_name']

    // Always include basic fields
    basicFields.each { field ->
        if (meta.containsKey(field)) {
            relevantParams[field] = meta[field]
        }
    }

    // Add parameters based on analysis type
    switch(analysis_type) {
        case 'differential':
            // Include differential-specific parameters
            meta.each { key, value ->
                if (differentialPrefixes.any { key.startsWith(it) }) {
                    relevantParams[key] = value
                }
            }
            break

        case 'functional':
            // Include functional-specific parameters
            meta.each { key, value ->
                if (functionalPrefixes.any { key.startsWith(it) } ||
                    (meta.functional_method == 'gsea' && gseaPrefixes.any { key.startsWith(it) }) ||
                    (meta.functional_method == 'gprofiler2' && gprofiler2Prefixes.any { key.startsWith(it) })) {
                    relevantParams[key] = value
                }
            }

            // For gprofiler2, also include differential parameters
            if (meta.functional_method == 'gprofiler2') {
                meta.each { key, value ->
                    if (differentialPrefixes.any { key.startsWith(it) }) {
                        relevantParams[key] = value
                    }
                }
            }
            break
    }
    return relevantParams
}

// Helper function to re-apply full parameters to deduplicated results
def reapplyParams(channel, paramsets_by_name) {
    channel
        .map{meta, differential_results -> [meta.analysis_names, meta, differential_results]}
        .transpose()
        .join(paramsets_by_name)
        .map{analysis_names, meta, differential_results, paramset -> [meta + ['study_meta': paramset], differential_results]}
}
