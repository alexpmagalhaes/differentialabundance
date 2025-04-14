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
    // Get paramsets based on toolsheet or default parameters
    //
    paramsets = getToolConfigurations()
    ch_paramsets = Channel.fromList(paramsets)

    //
    // Custom validate input parameters
    //
    validateInputParameters(paramsets)

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
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters(paramsets) {

    genomeExistsError()

    // check the params in each of the paramsets
    paramsets.each { row ->

        // check file existence for input parameters
        if (row.input) {
            file(row.input, checkIfExists: true)
        }

        // Validate study type spepecific parameters
        if (row.study_type == 'affy_array') {
            if (!row.affy_cel_files_archive) {
                error("CEL files archive not specified!")
            }
        } else if (row.study_type == 'maxquant') {
            if (row.functional_method) {
                error("Functional analysis is not yet possible with maxquant input data; please set --functional_method to null and rerun pipeline!")
            }
            if (!row.matrix) {
                error("Input matrix not specified!")
            }
        } else if (row.study_type == 'geo_soft_file') {
            if (!row.querygse || !row.features_metadata_cols) {
                error("Query GSE not specified or features metadata columns not specified")
            }
        } else if (row.study_type == "rnaseq") {
            if (!row.matrix) {
                error("Input matrix not specified!")
            }
        }

        // Validate functional analysis parameters
        if (row.functional_method) {
            if (row.functional_method == 'gsea' && !row.gene_sets_files) {
                error("GSEA activated but gene set file not specified!")
            } else if (row.functional_method == 'gprofiler2') {
                if (!row.gprofiler2_token && !row.gprofiler2_organism) {
                    error("To run gprofiler2, please provide a run token, GMT file or organism!")
                }
                if (row.gene_sets_files && row.gene_sets_files.split(",").size() > 1) {
                    error("gprofiler2 can currently only work with a single gene set file")
                }
            }
        }

        // Validate contrasts parameters
        if (row.contrasts_yml && row.contrasts) {
            error("Both '--contrasts' and '--contrasts_yml' parameters are set. Please specify only one of these options to define contrasts.")
        }
        if (!(row.contrasts_yml || row.contrasts)) {
            error("Either '--contrasts' and '--contrasts_yml' must be set. Please specify one of these options to define contrasts.")
        }
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
// To be able to retrieve a full set of parameters considering both toolsheet rows
// (highest priority) and pipeline params scope and validate them all together,
// we need to do the following:
// 1. get and validate toolsheet configurations only
// 2. get and validate toolsheet configurations merged with pipeline params
// The second step involves creating a temporary json file with the full paramsets
// and a temporary schema file with all the parameters, and use samplesheetToList
// to validate them. In order to create this temporary json file, we need to first
// retrieve the toolsheet configurations with the proper schema. This is why this
// is done in two steps.
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
    def raw_toolsheet = samplesheetToList(toolsheet_path, schema_path).collect { it -> it[0] }

    def toolsheet = raw_toolsheet
        // remove empty values
        .collect { row ->
            return row.findAll { key, value -> value != [] }
        }
        // If analysis_name is set, only keep matching rows
        .findAll { row ->
            if (params.analysis_name) {
                return row.analysis_name == params.analysis_name
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
            // merge pipeline params with toolsheet params. Toolsheet params
            // override pipeline params as they have highest priority
            def paramset = params + row
            // Only include parameters that are not in the ignore list
            def ignore = ['help', 'help_full', 'show_hidden', 'genomes']
            return paramset.findAll { k,v -> !ignore.contains(k) } as Map
        }

    // write the full paramsets to a temporary json file
    def tempFile = new File("${projectDir}/extended.json")
    tempFile.text = new groovy.json.JsonBuilder(paramsets).toPrettyString()
    def full_toolsheet_path = tempFile.absolutePath

    // create temporary toolsheet schema with all the parameters
    def full_schema_path = createToolsheetSchema(null)

    // load and validate the full paramsets
    def paramsets_validated = samplesheetToList(full_toolsheet_path, full_schema_path).collect { it -> it[0] }

    return paramsets_validated
        .collect{ row ->
            // replace empty values [] by null
            def cleanparams = row.collectEntries { k, v -> [(k): (v == [] ? null : v)] }
            // sort toolsheet params based on pipeline params order
            def paramKeysList = params.keySet().toList()
            def sortedparams = cleanparams.sort { a, b ->
                    paramKeysList.indexOf(a.key) <=> paramKeysList.indexOf(b.key)
                } as Map
            // use only the static params, otherwise it causes issues with cache
            def staticparams = sortedparams.findAll { k, v -> k != 'trace_report_suffix' } as Map
            return staticparams
        }
}

// Get default configurations from pipeline parameters
def getDefaultConfigurations() {
    def ignore = ['help', 'help_full', 'show_hidden', 'genomes']
    def nonstaticparams = ['trace_report_suffix']
    return [params.findAll { k, v -> !(ignore + nonstaticparams).contains(k) } as Map]
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
    }

    // Extract properties and required fields from schema
    def all_properties = extractPropertiesFromSchema(schema_json, keep_properties)
    def all_required = (toolsheet_path) ? ['analysis_name'] : extractRequiredFromSchema(schema_json)

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
    // def temp_schema = File.createTempFile("samplesheet_schema", ".json")
    def temp_schema = new File("${projectDir}/schema_tools.json")
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
                // if keep list is provided, only keep those properties
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
def extractRequiredFromSchema(schema) {
    def all_required = []
    schema.$defs.each { group_name, group_def ->
        if (group_def.required) {
            group_def.required.each { req_prop ->
                all_required << req_prop
            }
        }
    }
    return all_required
}

// prepare the input for the module by keeping only the relevant params
// in the meta and group channel by simplified meta and unique inputs
// @param channel: the input channel
// @param module: the simplified module name
// @param size: the number of elements in the channel
def prepareModuleInput(channel, module='base', size=2) {
    return channel
        .map {
            // clean the meta by keeping only the relevant params
            def meta = getRelevantMeta(it[0], module)
            // and add analysis name
            [meta] + it.tail() + [it[0].analysis_name]
        }
        // groupTuple by simplified meta and unique inputs
        // Note that the same simplified meta relevant to the module may have
        // input files originated from different upstream processes/params
        .groupTuple(by: 0..(size-1))  // -1 because it is 0-based indexing
        .map {
            def meta = it[0] + [analysis_name: it[-1]]
            [meta] + it[1..-2]
        }
}

// prepare the output for the module by adding back the full paramsets
// based on the analysis_name
def prepareModuleOutput(channel, paramsets_by_name) {
    channel_by_name = channel
        // first parse the channel to have analysis_name as key
        .map { [it[0].analysis_name] + it }
        // transpose the list of analysis_name
        .transpose()
        .map {
            // remap single analysis_name string into the meta
            def meta = it[1] + [analysis_name: it[0]]
            [it[0], meta] + it[2..-1]  // [analysis_name, meta, files ...]
        }
    return paramsets_by_name
        // we use combine instead of join, because for the same analysis_name,
        // different meta and files may have been generated because of the different contrasts
        .combine(channel_by_name, by:0)
        .map {
            // use full paramsets as meta
            def meta_full = it[1] + it[2]
            [meta_full] + it[3..-1]    // [meta with full paramsets, files ...]
        }
}

// get the relevant meta for the module
// @params meta: the meta map
// @params module: the simplified module name
def getRelevantMeta (meta, module) {
    def relevantParams = [:]

    // define the base keys and prefix
    def keys_base = ['study_name', 'study_type']
    def keys_contrast = ['variable', 'reference', 'target', 'blocking', 'formula'] // also always keep contrast keys, if available
    def keys = ['id'] + keys_base + keys_contrast
    def prefix = []

    // add prefix and keys based on module
    switch(module) {
        case 'base':
            break
        case 'affy':
            prefix += ['affy_']
            keys += ['observations_id_col','observations_name_col']
            break
        case 'proteus':
            prefix += ['proteus_']
            keys += ['observations_id_col','features_id_col','report_round_digits']
            break
        case 'geoquery':
            keys += ['features_metadata_cols']
            break
        case 'gtf':
            prefix += ['features_gtf_']
            break
        case 'validator':
            keys += ['observations_id_col','features_id_col']
            break
        case 'matrixfilter':
            prefix += ['filtering_']
            keys += ['observations_id_col']
            break
        case 'differential':
            prefix += ['differential_', meta.differential_method+'_', 'exclude_samples_']
            keys += ['features_id_col', 'observations_id_col', 'report_round_digits', 'sizefactors_from_controls']
            break
        case 'functional':
            prefix += ['functional_', meta.functional_method+'_']
            keys += ['features_id_col', 'features_name_col']
            break
        case 'exploratory':
            prefix += ['exploratory_']
            keys += ['features_id_col', 'observations_id_col']
            break
        case 'plot_differential':
            prefix += ['differential_']
            keys += ['features_id_col']
            break
        case 'shiny':
            prefix += ['shinyngs_']
            break
        case 'report':
            prefix += ['report_']
            keys += ['logo_file', 'css_file', 'citations_file']
            break
        case 'immunedeconv':
            prefix += ['immunedeconv_']
            break
        default:
            error("Module '${module}' not recognized by getRelevantMeta.")
    }

    // get key, value pairs from meta that start with the prefix
    // or contain the keys in the keys list
    meta.each { k, v ->
        if ( (prefix.any { k.startsWith(it) }) || (keys.contains(k)) ) {
            relevantParams[k] = v
        }
    }
    return relevantParams
}

// get the meta excluding the keys related to contrast
def getMetaWithoutContrast(meta) {
    def key_contrast = ['id', 'variable', 'reference', 'target', 'blocking', 'formula']
    return meta.findAll { key, value ->
        ! key_contrast.contains(key)
    }
}
