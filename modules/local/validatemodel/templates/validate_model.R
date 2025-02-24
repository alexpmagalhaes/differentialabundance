#!/usr/bin/env Rscript

#' Model Validation for Differential Expression Analysis
#'
#' This script validates experimental design models for differential expression analysis.
#' It processes a YAML file containing model definitions and validates them against a
#' sample sheet (CSV/TSV format). The script performs several validation steps:
#'
#' 1. Validates model definitions from YAML
#' 2. Checks sample sheet compatibility
#' 3. Verifies model matrix rank
#' 4. Exports validated phenotype data and design matrices
#'
#' The script is designed to work within a Nextflow pipeline, handling template variables
#' appropriately.
#'
#' @author Alan Mobbs
#' @version 1.0
#' @date 2024

# LOAD LIBRARIES ----------------------------------------------------
suppressWarnings(suppressMessages({
    library(tidyverse)
    library(yaml)
    library(jsonlite)
}))

#' Parse Command Line Arguments
#'
#' Parses command line arguments in the format --arg value into a named list
#'
#' @param x String containing command line arguments
#' @return Named list of argument values
#' @examples
#' parse_args("--input file.txt --output result.txt")
parse_args <- function(x) {
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})
    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[ ( ! parsed_args %in%  c('', 'null')) & ! is.na(parsed_args)]
}

#' Process Model Definitions from YAML
#'
#' Processes model definitions from a YAML file to extract contrasts and variables
#'
#' @param models List containing model definitions from YAML
#' @return List containing contrasts_list and var components
#' @examples
#' models <- yaml::read_yaml("models.yml")
#' results <- process_models(models)
process_models <- function(models) {
    # Initialize lists to store the output
    contrasts_list <- list()
    var <- list()

    # Ensure models\$models is not null or empty
    if (is.null(models\$contrasts) || length(models\$contrasts) == 0) {
        stop("models\$contrasts is null or empty. Please provide valid input.")
    }

    # Temporary storage for gathering contrasts per variable
    temp_var <- list()
    blocking_vars <- c()

    for (contrast in models\$contrasts) {
        # Validate contrast\$id
        if (is.null(contrast\$id)) {
            stop("Missing contrast\$id detected.")
        }

        if (contrast\$id %in% names(contrasts_list)) {
            stop(paste0("Duplicate contrast\$id detected: ", contrast\$id, "."))
        }

        # Get column name: first comparison's element
        variable <- contrast\$comparison[1]

        # Check if name is valid
        if (is.null(variable) || variable == "") {
            stop("Invalid `variable` defined in `contrast\$comparison[1]`.")
        }

        # Get factors (rest of components)
        levels <- contrast\$comparison[2:length(contrast\$comparison)]

        # Populate contrasts_list for later model validation
        contrasts_list[[contrast\$id]] <- list(
            "variable" = variable,
            "contrast" = levels,
            "blocking_factors" = contrast\$blocking_factors
        )

        # Gather contrasts per variable into temporary storage
        if (variable %in% names(temp_var)) {
            temp_var[[variable]] <- unique(c(temp_var[[variable]], levels))
        } else {
            temp_var[[variable]] <- levels
        }

        # Gather blocking factors
        if (!is.null(contrast\$blocking_factors)) {
            blocking_vars <- unique(c(blocking_vars, contrast\$blocking_factors))
        }
    }

    # Consolidate gathered contrasts into `var`
    for (variable in names(temp_var)) {
        if (variable %in% names(var)) {
            var[[variable]] <- unique(c(var[[variable]], temp_var[[variable]]))
        } else {
            var[[variable]] <- temp_var[[variable]]
        }
    }

    # Add blocking variables to `var`
    var[[ "blocking_factors" ]] <- blocking_vars

    # Return the populated lists
    return(list(
        contrasts_list = contrasts_list,
        var = var
    ))
}

#' Validate Model Against Samplesheet
#'
#' Validates that the model definitions match the data in the samplesheet
#'
#' @param sample_column String name of the sample identifier column
#' @param variables List of variables and their expected levels
#' @param samplesheet Data frame containing the sample data
#' @return List containing validated phenotype table and warnings
validate_model <- function(sample_column, variables, samplesheet) {
    # Initialize vectors for errors, warnings and continuous variables
    errors <- character(0)
    warnings <- character(0)
    continuous <- character(0)

    # Regex for invalid characters
    undesired_chars <- regex("[^a-zA-Z0-9_.]")

    # Validate column names
    invalid_cols <- stringr::str_detect(names(samplesheet), pattern = undesired_chars)
    if (any(invalid_cols)) {
        errors <- c(errors,
            sprintf("The following columns contain undesired characters: %s",
                paste(names(samplesheet)[invalid_cols], collapse = " ")))
    }

    # Validate blocking factors if present
    if (!is.null(variables\$blocking_factors)) {
        blocking_factors <- variables\$blocking_factors
        for (var_name in blocking_factors) {
            if (var_name %in% colnames(samplesheet)) {
                na_rows <- which(is.na(samplesheet[[var_name]]))
                if (length(na_rows) > 0) {
                    errors <- c(errors,
                        sprintf("Blocking factor %s contains NA/s in rows: %s",
                            var_name, paste(na_rows, collapse = " ")))
                }
                if (is.numeric(samplesheet[[var_name]])) {
                    continuous <- c(continuous, var_name)
                }
            } else {
                errors <- c(errors,
                    sprintf("Blocking factor %s not present in sample sheet", var_name))
            }
        }
    }

    # Validate main variables
    main_vars <- setdiff(names(variables), "blocking_factors")
    for (var_name in main_vars) {
        if (var_name %in% colnames(samplesheet)) {
            # Check for continuous variables
            if (is.numeric(samplesheet[[var_name]])) {
                continuous <- c(continuous, var_name)
            }

            # Check for NAs
            na_rows <- which(is.na(samplesheet[[var_name]]))
            if (length(na_rows) > 0) {
                errors <- c(errors,
                    sprintf("Column %s contains NA/s in rows: %s",
                        var_name, paste(na_rows, collapse = " ")))
            }

            # Check for invalid characters
            if (any(stringr::str_detect(samplesheet[[var_name]], pattern = undesired_chars), na.rm = TRUE)) {
                errors <- c(errors,
                    sprintf("Column %s contains undesired characters", var_name))
            }

            # Validate levels
            expected_levels <- variables[[var_name]]
            actual_levels <- unique(samplesheet[[var_name]])
            missing_levels <- setdiff(expected_levels, actual_levels)
            extra_levels <- setdiff(actual_levels, expected_levels)

            if (length(missing_levels) > 0) {
                errors <- c(errors,
                    sprintf("Missing factor levels for variable '%s'. Present: %s. Missing: %s",
                        var_name, paste(setdiff(expected_levels, missing_levels), collapse = " "),
                        paste(missing_levels, collapse = " ")))
            }

            if (length(extra_levels) > 0) {
                warnings <- c(warnings,
                    sprintf("Labels present in samplesheet but not in models definition: %s",
                        paste(extra_levels, collapse = " ")))
            }
        } else {
            errors <- c(errors,
                sprintf("Column %s not present in sample sheet", var_name))
        }
    }

    # Handle errors
    if (length(errors) > 0) {
        stop(cat("\033[1;31mErrors in samplesheet validation:\n",
                paste(errors, collapse = "\n"), "\033[0m\n", sep = ""))
    }

    # Create validated phenotype table
    selected_columns <- c(sample_column, main_vars, variables\$blocking_factors)
    pheno_table <- samplesheet %>% dplyr::select(all_of(selected_columns))

    # Add continuous variable warning if needed
    if (length(continuous) > 0) {
        warnings <- c(warnings,
            sprintf("Continuous variables detected: %s",
                paste(continuous, collapse = " ")))
    }

    return(list(pheno_table = pheno_table, warnings = warnings))
}

#' Check Model Contrasts
#'
#' Validates that the model contrasts are full ranked
#'
#' @param contrasts_list List of contrast definitions
#' @param colData Data frame containing the phenotype data
#' @return List of design matrices and their rank status
check_model_contrasts <- function(contrasts_list, colData) {
    design_list <- list()

    for (model_name in names(contrasts_list)) {
        model <- contrasts_list[[model_name]]
        formula_terms <- c(model\$variable, model\$blocking_factors)
        dynamic_formula <- as.formula(paste("~", paste(formula_terms, collapse = "+")))

        cat("\nChecking:", model_name, "\nFormula:", deparse(dynamic_formula), "\n")

        design_matrix <- model.matrix(dynamic_formula, data = colData)
        rank <- qr(design_matrix)\$rank
        expected_rank <- ncol(design_matrix)

        design_list[[model_name]] <- list(
            formula = deparse(dynamic_formula),
            design_matrix = design_matrix,
            full_rank = rank == expected_rank
        )
    }

    # Print rank status
    for (design_name in names(design_list)) {
        cat(sprintf("Design %s %s full ranked\n",
                design_name,
                if(design_list[[design_name]]\$full_rank) "is" else "is not"))
    }

    return(design_list)
}

# PARSE ARGUMENTS ---------------------------------------------------
opt <- list(
    models_yml = '$models_yml',
    samplesheet = '$samplesheet',
    sample_id_col = "sample1"
)

# Process command line arguments
args_opt <- parse_args('$task.ext.args')
for (ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    }
    opt[[ao]] <- args_opt[[ao]]
}

# Validate required options
required_opts <- c("models_yml", "samplesheet")
missing_opts <- required_opts[!required_opts %in% names(opt) | sapply(opt[required_opts], function(x) is.null(x) || x == "null")]
if (length(missing_opts) > 0) {
    stop(paste("Missing required options:", paste(missing_opts, collapse = ", ")))
}
cat("DEBUG: Final options: ", paste(names(opt), opt, sep="=", collapse=", "), "\n")

# Set default sample_id_col if not provided
if (!("sample_id_col" %in% names(opt))) {
    cat("DEBUG: sample_id_col not found in options; defaulting to 'sample'\n")
    opt\$sample_id_col <- "sample"
}

# Validate input paths
if (is.null(opt\$models_yml) || is.null(opt\$samplesheet) ) {
    stop("'--yml' and '--samplesheet' arguments are required.", call. = FALSE)
}

path_yml         <- opt\$models_yml
path_samplesheet <- opt\$samplesheet
sample_column    <- opt\$sample_id_col

# LOAD FILES --------------------------------------------------------
# Load and validate models.yml file
tryCatch({
    if (!file.exists(path_yml)) {
        stop(sprintf("File '%s' does not exist.", path_yml), call. = FALSE)
    }
    if (!(endsWith(tolower(path_yml), ".yaml") || endsWith(tolower(path_yml), ".yml"))) {
        stop(sprintf("File '%s' is not a YAML file (expected .yaml or .yml extension).", path_yml), call. = FALSE)
    }
    models <- yaml::read_yaml(path_yml)
    cat("Loaded YML file successfully.\n")
}, error = function(e) {
    stop(sprintf("Error loading YML file '%s': %s", path_yml, e\$message), call. = FALSE)
})

# Load and validate samplesheet file
tryCatch({
    # Set the separator based on the file extension
    sep <- if (endsWith(tolower(path_samplesheet), ".csv")) {
        ","
    } else if (endsWith(tolower(path_samplesheet), ".tsv")) {
        "\t"
    } else {
        stop("Unsupported file extension. Please provide a .csv or .tsv file.")
    }

    # Read the file using the determined separator
    samplesheet <- read_delim(path_samplesheet, delim = sep, show_col_types = FALSE)
    if ( !sample_column %in% colnames(samplesheet) ) {
        stop(paste0("The sample identificator column '", sample_column, "' is not present in the samplesheet file") )
    }
    cat("Loaded samplesheet successfully.\n")
}, error = function(e) {
    stop("Error loading samplesheet CSV: ", e\$message)
})

# PARSE FACTORS/LEVELS FROM YML FILE --------------------------------
# Process model definitions and extract contrasts and variables
models_lists <- process_models(models)
contrasts_list <- models_lists[["contrasts_list"]]
var <- models_lists[["var"]]

# Print detected variables and their levels
for (var_idx in seq_len(length(var))) {
    cat(sprintf("\033[32mDetected '%s' variable with %s levels.\033[0m\n",
        names(var)[var_idx],
        paste(var[[var_idx]], collapse = " ")))
}

# VALIDATE SAMPLESHEET BASED ON YML FILE ----------------------------
# Validate model definitions against samplesheet data
phenotable_warnings <- validate_model(sample_column, var, samplesheet)

# Display any warnings from validation
cat("\033[1;33m", paste(phenotable_warnings[[2]], collapse = "\n"), "\033[0m\n")

# Extract validated phenotype table
pheno_table <- phenotable_warnings[[1]]

# Validate model matrix rank
design_results <- check_model_contrasts(contrasts_list, pheno_table)

# EXPORT DATA -------------------------------------------------------
# Export validated data and results
write_csv(pheno_table, "pheno_table.csv")
jsonlite::write_json(design_results,
    path = "designs.json",
    pretty = TRUE,
    auto_unbox = TRUE
)

# Export any warnings to JSON
if (!is.null(phenotable_warnings[[2]])) {
    jsonlite::write_json(phenotable_warnings[[2]],
        path = "warnings.json",
        pretty = TRUE
    )
}

# Save complete workspace
save.image("models.RData")

# SOFTWARE VERSIONS ------------------------------------------------
r_version <- strsplit(version[['version.string']], ' ')[[1]][3]
yaml_version <- as.character(packageVersion('yaml'))
jsonlite_version <- as.character(packageVersion('jsonlite'))
tidyverse_version <- as.character(packageVersion('tidyverse'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r_version),
        paste('    yaml:', yaml_version),
        paste('    jsonlite:', jsonlite_version),
        paste('    tidyverse:', tidyverse_version)
    ),
    'versions.yml'
)
