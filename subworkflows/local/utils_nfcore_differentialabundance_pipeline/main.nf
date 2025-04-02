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
    ch_tools = Channel.fromList(getToolConfigurations())

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

// Create a temporary schema file for toolsheet validation
// This schema is derived from the pipeline's nextflow_schema.json by:
// 1. Reading the toolsheet headers to determine which fields to include
// 2. Extracting only those properties from the pipeline schema
// 3. Creating a schema that allows for an array of objects
// 4. Making analysis_name the only required field
// 5. Adding meta fields to all properties
// @return The absolute path to the temporary schema file
def createToolsheetSchema() {
    // Load toolsheet and get headers
    def toolsheet_path = params.toolsheet_custom ?: "${projectDir}/assets/toolsheet.csv"
    def toolsheet_lines = new File(toolsheet_path).readLines()
    def headers = toolsheet_lines[0].split(',')

    // Ensure analysis_name is in headers
    if (!headers.contains('analysis_name')) {
        error("The toolsheet must contain an 'analysis_name' column")
    }

    // Load and parse pipeline schema
    def pipeline_schema = new File("${projectDir}/nextflow_schema.json").text
    def schema_json = new groovy.json.JsonSlurper().parseText(pipeline_schema)

    // Extract only the properties that are present in the toolsheet
    def all_properties = [:]
    schema_json.$defs.each { group_name, group_def ->
        if (group_def.properties) {
            group_def.properties.each { prop_name, prop_def ->
                // Only include properties that are in the toolsheet headers
                if (headers.contains(prop_name)) {
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
        schema_json.$defs.each { group_name, group_def ->
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

// Get configurations from toolsheet
def getToolsheetConfigurations() {
    // Load toolsheet
    def toolsheet_path = params.toolsheet_custom ?: "${projectDir}/assets/toolsheet.csv"

    // Create temporary schema for validation
    def schema_path = createToolsheetSchema()

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

    return toolsheet.collect{ row -> (params.findAll { k, v -> k != 'trace_report_suffix' } + row[0]) }
}

// Get default configurations from pipeline parameters
def getDefaultConfigurations() {
    return [params]
}
