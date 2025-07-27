#' @title Naming Standards and Conventions
#' @description Functions to ensure consistent naming throughout the package
#' @author Deniz Akdemir

# =============================================================================
# NAMING CONVENTION UTILITIES
# =============================================================================

#' Convert function names to snake_case
#' 
#' @param name Function or variable name
#' @return Standardized snake_case name
#' @export
standardize_name <- function(name) {
  # Convert camelCase to snake_case
  name <- gsub("([a-z])([A-Z])", "\\1_\\2", name)
  
  # Convert to lowercase
  name <- tolower(name)
  
  # Replace multiple underscores with single
  name <- gsub("_{2,}", "_", name)
  
  # Remove leading/trailing underscores
  name <- gsub("^_|_$", "", name)
  
  return(name)
}

#' Check if parameter names follow conventions
#' 
#' @param param_names Vector of parameter names
#' @return List with standardized names and warnings
#' @export
check_parameter_naming <- function(param_names) {
  results <- list(
    original = param_names,
    standardized = sapply(param_names, standardize_name),
    warnings = character(0)
  )
  
  # Check for common issues
  non_standard <- param_names[param_names != results$standardized]
  if (length(non_standard) > 0) {
    results$warnings <- c(results$warnings, 
                         paste("Non-standard parameter names:", 
                               paste(non_standard, collapse = ", ")))
  }
  
  # Check for reserved words
  reserved_words <- c("if", "else", "repeat", "while", "function", "for", "in", "next", "break",
                      "TRUE", "FALSE", "NULL", "Inf", "NaN", "NA", "NA_integer_", "NA_real_", 
                      "NA_complex_", "NA_character_")
  
  reserved_used <- intersect(param_names, reserved_words)
  if (length(reserved_used) > 0) {
    results$warnings <- c(results$warnings,
                         paste("Reserved words used as parameters:", 
                               paste(reserved_used, collapse = ", ")))
  }
  
  return(results)
}

# =============================================================================
# PARAMETER STANDARDIZATION
# =============================================================================

#' Standardize common parameter names across functions
#' 
#' @param ... Named parameters to standardize
#' @return List with standardized parameter names and values
#' @export
standardize_parameters <- function(...) {
  params <- list(...)
  
  # Common parameter mappings
  standard_mappings <- list(
    # Core parameters
    "mc.cores" = "mc_cores",
    "mc.cleanup" = "mc_cleanup",
    "mc.silent" = "mc_silent",
    
    # Matrix parameters
    "lambda" = "ridge_param",  # Keep lambda for compatibility
    "K" = "kinship_matrix",
    "P" = "prediction_matrix",
    
    # Genetic algorithm parameters
    "npop" = "n_population",
    "ngen" = "n_generations", 
    "nelite" = "n_elite",
    "pmut" = "mutation_prob",
    "pcross" = "crossover_prob",
    
    # Criteria parameters
    "Vg" = "genetic_variance",
    "Ve" = "error_variance",
    "h2" = "heritability"
  )
  
  # Apply mappings while preserving original names for compatibility
  param_names <- names(params)
  for (i in seq_along(params)) {
    old_name <- param_names[i]
    if (old_name %in% names(standard_mappings)) {
      new_name <- standard_mappings[[old_name]]
      message("Parameter '", old_name, "' mapped to standard name '", new_name, "'")
    }
  }
  
  return(params)
}

# =============================================================================
# FUNCTION NAMING CONSISTENCY
# =============================================================================

#' Get consistent function name variants
#' 
#' @param base_name Base function name
#' @return List with name variants for different contexts
#' @export
get_function_variants <- function(base_name) {
  base_clean <- standardize_name(base_name)
  
  list(
    modern = base_clean,
    legacy = toupper(gsub("_", "", base_clean)),
    camel = gsub("_(.)", "\\U\\1", base_clean, perl = TRUE),
    abbreviated = abbreviate(base_clean, minlength = 6)
  )
}

#' Create function alias for backward compatibility
#' 
#' @param old_name Old function name
#' @param new_name New standardized function name
#' @param envir Environment to create alias in
#' @export
create_function_alias <- function(old_name, new_name, envir = .GlobalEnv) {
  if (exists(new_name, envir = envir)) {
    new_function <- get(new_name, envir = envir)
    
    # Create wrapper function with deprecation warning
    wrapper <- function(...) {
      .Deprecated(new_name, package = "STPGA", 
                  msg = paste0("Function '", old_name, "' is deprecated. Use '", 
                               new_name, "' instead."))
      do.call(new_function, list(...))
    }
    
    # Copy attributes from original function
    attributes(wrapper) <- attributes(new_function)
    
    assign(old_name, wrapper, envir = envir)
    message("Created alias: ", old_name, " -> ", new_name)
  } else {
    warning("Cannot create alias: function '", new_name, "' does not exist")
  }
}

# =============================================================================
# DOCUMENTATION HELPERS
# =============================================================================

#' Generate consistent parameter documentation
#' 
#' @param param_name Parameter name
#' @param param_type Parameter type
#' @param default_value Default value
#' @param description Parameter description
#' @return Formatted roxygen2 parameter documentation
#' @export
document_parameter <- function(param_name, param_type, default_value = NULL, description) {
  # Standardize the name
  std_name <- standardize_name(param_name)
  
  # Format type information
  type_info <- switch(param_type,
    "matrix" = "Matrix",
    "vector" = "Vector", 
    "numeric" = "Numeric scalar",
    "integer" = "Integer scalar",
    "logical" = "Logical scalar",
    "character" = "Character scalar",
    "list" = "List",
    param_type
  )
  
  # Include default if provided
  if (!is.null(default_value)) {
    type_info <- paste0(type_info, " (default: ", deparse(default_value), ")")
  }
  
  # Format for roxygen2
  paste0("#' @param ", std_name, " ", type_info, ". ", description)
}

#' Generate standard function documentation template
#' 
#' @param func_name Function name
#' @param purpose Function purpose
#' @param category Function category
#' @return Roxygen2 documentation template
#' @export
generate_doc_template <- function(func_name, purpose, category = "General") {
  std_name <- standardize_name(func_name)
  
  template <- paste0(
    "#' @title ", tools::toTitleCase(gsub("_", " ", std_name)), "\n",
    "#' @description ", purpose, "\n",
    "#' @author Deniz Akdemir\n",
    "#' @category ", category, "\n",
    "#' @export\n",
    std_name, " <- function() {\n",
    "  # Function implementation here\n",
    "}\n"
  )
  
  return(template)
}