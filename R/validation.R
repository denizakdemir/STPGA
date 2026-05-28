#' @title Input Validation and Utility Functions
#' @description Comprehensive validation functions for STPGA package
#' @author Deniz Akdemir

# =============================================================================
# VALIDATION HELPER FUNCTIONS
# =============================================================================

#' Validate matrix parameters for optimization criteria
#' 
#' @param P Prediction matrix (individuals x variables)
#' @param K Kinship matrix (optional, for mixed models)
#' @param lambda Ridge regularization parameter
#' @param train Training set names
#' @param test Test set names (optional)
#' @param C Contrast matrix (optional)
#' @return NULL if valid, throws error if invalid
#' @export
validate_matrix_params <- function(P, K = NULL, lambda = 1e-6, train = NULL, test = NULL, C = NULL) {
  
  # Check P matrix
  if (!is.matrix(P) && !is.data.frame(P)) {
    stop("P must be a matrix or data.frame")
  }
  
  if (any(!is.finite(P))) {
    stop("P contains non-finite values (NA, NaN, or Inf)")
  }
  
  # Check dimensions
  if (nrow(P) < 2 || ncol(P) < 1) {
    stop("P must have at least 2 rows and 1 column")
  }
  
  # Check row names
  if (is.null(rownames(P))) {
    stop("P must have row names to identify individuals")
  }
  
  # Check for duplicate row names
  if (any(duplicated(rownames(P)))) {
    stop("P has duplicate row names")
  }
  
  # Check lambda
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop("lambda must be a non-negative numeric scalar")
  }
  
  if (lambda > 1) {
    warning("Large lambda value (", lambda, ") may over-regularize the problem")
  }
  
  # Check training set
  if (!is.null(train)) {
    if (!all(train %in% rownames(P))) {
      missing_train <- train[!train %in% rownames(P)]
      stop("Training set contains individuals not in P: ", paste(missing_train, collapse = ", "))
    }
    
    if (length(train) < 2) {
      stop("Training set must contain at least 2 individuals")
    }
    
    if (length(unique(train)) != length(train)) {
      stop("Training set contains duplicate individuals")
    }
  }
  
  # Check test set
  if (!is.null(test)) {
    if (!all(test %in% rownames(P))) {
      missing_test <- test[!test %in% rownames(P)]
      stop("Test set contains individuals not in P: ", paste(missing_test, collapse = ", "))
    }
    
    if (length(unique(test)) != length(test)) {
      stop("Test set contains duplicate individuals")
    }
    
    # Check for overlap between train and test
    if (!is.null(train) && any(train %in% test)) {
      overlap <- intersect(train, test)
      stop("Training and test sets overlap: ", paste(overlap, collapse = ", "))
    }
  }
  
  # Check kinship matrix
  if (!is.null(K)) {
    validate_kinship_matrix(K, rownames(P))
  }
  
  # Check contrast matrix
  if (!is.null(C)) {
    if (!is.matrix(C)) {
      stop("C must be a matrix")
    }
    
    if (ncol(C) != ncol(P)) {
      stop("Contrast matrix C must have same number of columns as P")
    }
    
    if (any(!is.finite(C))) {
      stop("Contrast matrix C contains non-finite values")
    }
  }
  
  invisible(NULL)
}

#' Validate kinship matrix
#' 
#' @param K Kinship matrix
#' @param individual_names Names of individuals to check against
#' @return NULL if valid, throws error if invalid
#' @export
validate_kinship_matrix <- function(K, individual_names) {
  
  if (!is.matrix(K)) {
    stop("K must be a matrix")
  }
  
  # Check dimensions
  if (nrow(K) != ncol(K)) {
    stop("Kinship matrix K must be square")
  }
  
  # Check row/column names
  if (is.null(rownames(K)) || is.null(colnames(K))) {
    stop("Kinship matrix K must have row and column names")
  }
  
  if (!identical(rownames(K), colnames(K))) {
    stop("Kinship matrix K must have identical row and column names")
  }
  
  # Check for required individuals
  if (!all(individual_names %in% rownames(K))) {
    missing_individuals <- individual_names[!individual_names %in% rownames(K)]
    stop("Kinship matrix K missing individuals: ", paste(missing_individuals, collapse = ", "))
  }
  
  # Check for finite values
  if (any(!is.finite(K))) {
    stop("Kinship matrix K contains non-finite values")
  }
  
  # Check symmetry
  if (!isSymmetric(K, tol = 1e-10)) {
    stop("Kinship matrix K is not symmetric")
  }
  
  # Check positive semi-definiteness
  eigenvals <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  min_eigenval <- min(eigenvals)
  
  if (min_eigenval < -1e-10) {
    stop("Kinship matrix K is not positive semi-definite (min eigenvalue: ", 
         sprintf("%.2e", min_eigenval), ")")
  }
  
  # Warning for near-singular matrices
  condition_number <- max(eigenvals) / max(min_eigenval, 1e-15)
  if (condition_number > 1e12) {
    warning("Kinship matrix K is ill-conditioned (condition number: ", 
            sprintf("%.2e", condition_number), ")")
  }
  
  # Check diagonal values (should be close to 1 for relationship matrices)
  diag_vals <- diag(K)
  if (any(diag_vals < 0.5 | diag_vals > 2)) {
    warning("Kinship matrix diagonal values outside expected range [0.5, 2]")
  }
  
  invisible(NULL)
}

#' Check memory usage and warn for large computations
#' 
#' @param n_individuals Number of individuals
#' @param n_variables Number of variables/markers
#' @param operation Type of operation ("matrix_mult", "eigen", "svd", "chol")
#' @return NULL, but issues warnings for large computations
#' @export
check_memory_usage <- function(n_individuals, n_variables, operation = "matrix_mult") {
  
  # Estimate memory usage in bytes
  bytes_per_double <- 8
  
  memory_estimates <- switch(operation,
    "matrix_mult" = {
      # X'X computation: n_variables^2 * 8 bytes
      bytes_per_double * n_variables^2
    },
    "eigen" = {
      # Eigendecomposition: roughly 3 * n^2 * 8 bytes for workspace
      bytes_per_double * 3 * max(n_individuals, n_variables)^2
    },
    "svd" = {
      # SVD: roughly 5 * max(m,n) * min(m,n) * 8 bytes
      bytes_per_double * 5 * max(n_individuals, n_variables) * min(n_individuals, n_variables)
    },
    "chol" = {
      # Cholesky: n^2 * 8 bytes
      bytes_per_double * max(n_individuals, n_variables)^2
    },
    bytes_per_double * n_individuals * n_variables  # default
  )
  
  # Convert to readable units
  if (memory_estimates > 1e8) {  # > 100 MB
    memory_gb <- memory_estimates / 1e9
    warning("Large computation detected: estimated memory usage ", 
            sprintf("%.2f", memory_gb), " GB for ", operation, 
            ". Consider using subsets or parallel processing.")
  } else if (memory_estimates > 5e7) {  # > 50 MB
    memory_mb <- memory_estimates / 1e6
    message("Moderate computation: estimated memory usage ", 
            sprintf("%.0f", memory_mb), " MB for ", operation)
  }
  
  invisible(NULL)
}

#' Validate variance components
#' 
#' @param Vg Genetic variance component
#' @param Ve Error variance component  
#' @param h2 Heritability (optional, for cross-validation)
#' @return NULL if valid, throws error if invalid
#' @export
validate_variance_components <- function(Vg, Ve, h2 = NULL) {
  
  # Check individual components
  if (!is.numeric(Vg) || length(Vg) != 1 || Vg < 0) {
    stop("Vg must be a non-negative numeric scalar")
  }
  
  if (!is.numeric(Ve) || length(Ve) != 1 || Ve < 0) {
    stop("Ve must be a non-negative numeric scalar")
  }
  
  # Check that at least one is positive
  if (Vg == 0 && Ve == 0) {
    stop("At least one of Vg or Ve must be positive")
  }
  
  # Cross-validate with heritability if provided
  if (!is.null(h2)) {
    if (!is.numeric(h2) || length(h2) != 1 || h2 < 0 || h2 > 1) {
      stop("h2 must be a numeric scalar between 0 and 1")
    }
    
    total_var <- Vg + Ve
    computed_h2 <- Vg / total_var
    
    if (abs(computed_h2 - h2) > 1e-10) {
      stop("Inconsistent variance components: h2 = ", sprintf("%.6f", h2),
           " but Vg/(Vg+Ve) = ", sprintf("%.6f", computed_h2))
    }
  }
  
  invisible(NULL)
}

#' Get adaptive ridge parameter based on condition number
#' 
#' @param X Matrix to analyze
#' @param target_condition Target condition number (default: 1e8)
#' @param min_lambda Minimum lambda value (default: 1e-10)
#' @param max_lambda Maximum lambda value (default: 1e-2)
#' @return Adaptive ridge parameter
#' @export
get_adaptive_ridge <- function(X, target_condition = 1e8, min_lambda = 1e-10, max_lambda = 1e-2) {
  
  # Compute condition number of X'X
  if (is.matrix(X)) {
    xtx <- crossprod(X)
  } else {
    xtx <- X  # Assume X is already X'X
  }
  
  # Get eigenvalues for condition number
  eigenvals <- eigen(xtx, symmetric = TRUE, only.values = TRUE)$values
  max_eigenval <- max(eigenvals)
  min_eigenval <- max(min(eigenvals), 1e-15)  # Avoid division by zero
  
  current_condition <- max_eigenval / min_eigenval
  
  if (current_condition <= target_condition) {
    return(min_lambda)  # Matrix is well-conditioned
  }
  
  # Compute required ridge parameter
  # We want: (max_eigenval) / (min_eigenval + lambda) = target_condition
  # So: lambda = max_eigenval/target_condition - min_eigenval
  required_lambda <- max_eigenval / target_condition - min_eigenval
  
  # Clamp to reasonable range
  adaptive_lambda <- max(min_lambda, min(max_lambda, required_lambda))
  
  if (adaptive_lambda > 1e-4) {
    message("Using adaptive ridge parameter: ", sprintf("%.2e", adaptive_lambda),
            " (condition number: ", sprintf("%.2e", current_condition), ")")
  }
  
  return(adaptive_lambda)
}

# =============================================================================
# DIMENSION CHECKING UTILITIES
# =============================================================================

#' Check matrix dimensions for operations
#' 
#' @param ... Named matrices to check
#' @return NULL if dimensions are compatible, throws error otherwise
#' @export
check_matrix_dimensions <- function(...) {
  matrices <- list(...)
  
  if (length(matrices) < 2) {
    return(invisible(NULL))
  }
  
  for (i in 1:(length(matrices) - 1)) {
    mat1 <- matrices[[i]]
    mat2 <- matrices[[i + 1]]
    name1 <- names(matrices)[i]
    name2 <- names(matrices)[i + 1]
    
    if (is.null(name1)) name1 <- paste0("matrix_", i)
    if (is.null(name2)) name2 <- paste0("matrix_", i + 1)
    
    # Check if they can be multiplied or added
    if (is.matrix(mat1) && is.matrix(mat2)) {
      if (ncol(mat1) != nrow(mat2) && !identical(dim(mat1), dim(mat2))) {
        stop("Incompatible dimensions between ", name1, " (", 
             nrow(mat1), "x", ncol(mat1), ") and ", name2, " (",
             nrow(mat2), "x", ncol(mat2), ")")
      }
    }
  }
  
  invisible(NULL)
}
