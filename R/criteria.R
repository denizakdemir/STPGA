#' @title STPGA Optimization Criteria Functions
#' @description Literature-corrected optimization criteria for subset selection in experimental design
#' @author Deniz Akdemir
#' 
#' @references 
#' Fedorov, V.V. (1972). Theory of Optimal Experiments. Academic Press.
#' Henderson, C.R. (1984). Applications of Linear Models in Animal Breeding.
#' Kiefer, J. (1959). Optimum experimental designs. Journal of the Royal Statistical Society B.
#' Searle, S.R., Casella, G., McCulloch, C.E. (1992). Variance Components. Wiley.
#' Atkinson, A.C., Donev, A.N., Tobias, R.D. (2007). Optimum Experimental Designs. Oxford.

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Convert between variance components and heritability
#' 
#' @param h2 Heritability (proportion of variance due to genetics)
#' @param total_var Total variance (default: 1)
#' @return List with Vg (genetic variance) and Ve (error variance)
#' 
#' @examples
#' # Convert heritability to variance components
#' h2_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' 
#' for (h2 in h2_values) {
#'   var_comps <- h2_to_variances(h2)
#'   print(paste("h² =", h2, "-> Vg =", var_comps$Vg, ", Ve =", var_comps$Ve))
#' }
#' 
#' # Example with custom total variance
#' var_comps_scaled <- h2_to_variances(0.6, total_var = 2.5)
#' print(paste("h² = 0.6, total_var = 2.5 -> Vg =", var_comps_scaled$Vg, 
#'             ", Ve =", var_comps_scaled$Ve))
#' 
#' @export
h2_to_variances <- function(h2, total_var = 1) {
  list(
    Vg = h2 * total_var,
    Ve = (1 - h2) * total_var
  )
}

#' Safe matrix inversion with numerical stability
#' 
#' @param X Matrix to invert
#' @param lambda Ridge parameter for regularization (default: 1e-6)
#' @param tol Tolerance for rank detection in SVD (default: 1e-12)
#' @return Inverse matrix
#' 
#' @examples
#' # Example with well-conditioned matrix
#' set.seed(123)
#' X <- matrix(rnorm(25), 5, 5)
#' X <- t(X) %*% X  # Make positive definite
#' 
#' # Safe inversion
#' X_inv <- safe_matrix_inverse(X)
#' print("Matrix inversion successful")
#' print(paste("Condition number:", round(kappa(X), 2)))
#' 
#' # Verify inversion
#' identity_check <- X %*% X_inv
#' max_error <- max(abs(identity_check - diag(nrow(X))))
#' print(paste("Maximum error in X * X_inv - I:", max_error))
#' 
#' # Example with ill-conditioned matrix (uses ridge regularization)
#' X_ill <- matrix(c(1, 1, 1, 1.0001), 2, 2)
#' print(paste("Ill-conditioned matrix condition number:", round(kappa(X_ill), 2)))
#' 
#' X_ill_inv <- safe_matrix_inverse(X_ill, lambda = 1e-4)
#' print("Ill-conditioned matrix inverted with ridge regularization")
safe_matrix_inverse <- function(X, lambda = 1e-6, tol = 1e-12) {
  # Input validation
  if (!is.matrix(X)) {
    stop("X must be a matrix")
  }
  
  if (nrow(X) != ncol(X)) {
    stop("X must be a square matrix")
  }
  
  if (any(!is.finite(X))) {
    stop("X contains non-finite values")
  }
  
  # Add ridge regularization if needed
  if (lambda > 0) {
    X_reg <- X + lambda * diag(nrow(X))
  } else {
    X_reg <- X
  }
  
  # Try Cholesky first (fastest for positive definite)
  chol_result <- tryCatch({
    L <- chol(X_reg)
    list(success = TRUE, inverse = chol2inv(L))
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
  
  if (chol_result$success) {
    return(chol_result$inverse)
  }
  
  # Fall back to SVD for numerical stability
  tryCatch({
    svd_result <- svd(X_reg)
    
    # Determine numerical rank
    max_singular_val <- max(svd_result$d)
    rank_threshold <- max(tol, max_singular_val * tol)
    
    # Only invert non-zero singular values
    d_inv <- ifelse(svd_result$d > rank_threshold, 1/svd_result$d, 0)
    
    # Check if matrix is effectively singular
    numerical_rank <- sum(svd_result$d > rank_threshold)
    if (numerical_rank < nrow(X)) {
      warning("Matrix is rank deficient (rank ", numerical_rank, " < ", nrow(X), 
              "). Using pseudo-inverse.")
    }
    
    # Compute pseudo-inverse
    svd_result$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(svd_result$u)
    
  }, error = function(e) {
    stop("Matrix inversion failed even with SVD: ", e$message)
  })
}

# =============================================================================
# CORE PREDICTION FUNCTIONS
# =============================================================================

#' Compute prediction core for matrix operations
#' 
#' This function computes the core matrix operations used by all optimization criteria.
#' It handles numerical stability and implements the correct PEV formula from literature.
#' 
#' @param p_train Training prediction matrix (n_train x p)
#' @param p_test Test prediction matrix (n_test x p). If NULL, uses training set
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @param adaptive_ridge Whether to use adaptive ridge parameter (default: TRUE)
#' @return List containing core matrices and computations
#' @export
compute_prediction_core <- function(p_train, p_test = NULL, lambda = 1e-6, C = NULL, adaptive_ridge = TRUE) {
  
  # Input validation
  if (!is.matrix(p_train)) {
    p_train <- as.matrix(p_train)
  }
  
  if (!is.null(p_test) && !is.matrix(p_test)) {
    p_test <- as.matrix(p_test)
  }
  
  # Check dimensions
  if (nrow(p_train) < 2) {
    stop("Training set must have at least 2 individuals")
  }
  
  if (ncol(p_train) < 1) {
    stop("Prediction matrix must have at least 1 variable")
  }
  
  if (!is.null(p_test) && ncol(p_test) != ncol(p_train)) {
    stop("Test and training matrices must have same number of columns")
  }
  
  # Check for memory usage
  check_memory_usage(nrow(p_train), ncol(p_train), "matrix_mult")
  
  # Compute XtX with regularization
  xtx <- crossprod(p_train)
  
  # Use adaptive ridge if requested
  if (adaptive_ridge) {
    adaptive_lambda <- get_adaptive_ridge(xtx)
    final_lambda <- max(lambda, adaptive_lambda)
    if (final_lambda > lambda) {
      message("Using adaptive ridge parameter: ", sprintf("%.2e", final_lambda))
    }
  } else {
    final_lambda <- lambda
  }
  
  xtx_reg <- xtx + final_lambda * diag(ncol(p_train))
  
  # Compute inverse with numerical stability
  inv_xtx_reg <- safe_matrix_inverse(xtx_reg, lambda = 0)  # Already regularized
  
  # CORRECTED prediction error variance computation
  # Based on literature: PEV = Var(y - y_hat) = sigma^2[I + X(X'X)^(-1)X']
  # The identity matrix accounts for inherent variability of new observations
  if (!is.null(p_test)) {
    # Test set: PEV includes both model uncertainty AND observation variance
    # FIXED: Ensure matrices are compatible for addition
    prediction_var <- p_test %*% inv_xtx_reg %*% t(p_test)
    
    # Verify dimensions match before addition
    n_test <- nrow(p_test)
    if (nrow(prediction_var) != n_test || ncol(prediction_var) != n_test) {
      stop("Internal error: prediction variance matrix has wrong dimensions")
    }
    
    # Create identity matrix of correct size
    identity_matrix <- diag(n_test)
    pev_matrix <- identity_matrix + prediction_var
    
  } else {
    # Training set: leave-one-out cross-validation
    prediction_var <- p_train %*% inv_xtx_reg %*% t(p_train)
    
    # Verify dimensions
    n_train <- nrow(p_train)
    if (nrow(prediction_var) != n_train || ncol(prediction_var) != n_train) {
      stop("Internal error: prediction variance matrix has wrong dimensions")
    }
    
    # Create identity matrix of correct size
    identity_matrix <- diag(n_train)
    pev_matrix <- identity_matrix + prediction_var
  }
  
  # Handle contrast matrix if provided
  contrast_core <- if (!is.null(C)) {
    # Validate contrast matrix dimensions
    if (ncol(C) != ncol(p_train)) {
      stop("Contrast matrix must have same number of columns as prediction matrix")
    }
    
    c_inv_xtx_ct <- C %*% inv_xtx_reg %*% t(C)
    
    # Check condition number
    condition_num <- kappa(c_inv_xtx_ct)
    if (condition_num > 1e12) {
      warning("Contrast matrix leads to ill-conditioned system (condition number: ",
              sprintf("%.2e", condition_num), ")")
    }
    
    safe_matrix_inverse(c_inv_xtx_ct)
  } else {
    NULL
  }
  
  list(
    inv_xtx_reg = inv_xtx_reg,
    pev_matrix = pev_matrix,
    p_train = p_train,
    p_test = p_test,
    contrast_core = contrast_core,
    C = C,
    lambda_used = final_lambda
  )
}

# =============================================================================
# CLASSICAL OPTIMALITY CRITERIA (CORRECTLY IMPLEMENTED)
# =============================================================================

#' A-optimality criterion for experimental design
#' 
#' Computes the A-optimality criterion, which minimizes the trace of the 
#' inverse covariance matrix. Lower values indicate better designs.
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and variables as columns
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @return A-optimality value (trace of inverse covariance matrix)
#' 
#' @details
#' A-optimality minimizes the average prediction variance. It is equivalent to 
#' minimizing the trace of (X'X)^(-1) where X is the design matrix.
#' 
#' Mathematical formula: A = tr((X'X)^(-1))
#' 
#' @references
#' Fedorov, V.V. (1972). Theory of Optimal Experiments. Academic Press.
#' 
#' @examples
#' # Load wheat genomic data
#' data(WheatData)
#' 
#' # Create a subset for faster computation
#' set.seed(123)
#' subset_indices <- sample(1:nrow(Wheat.M), 80)
#' M_subset <- Wheat.M[subset_indices, 1:20]
#' 
#' # Extract principal components
#' pca_result <- prcomp(M_subset, center = TRUE, scale. = TRUE)
#' PC_subset <- pca_result$x[, 1:5]
#' rownames(PC_subset) <- rownames(M_subset)
#' 
#' # Define candidate and test sets
#' all_individuals <- rownames(PC_subset)
#' test_set <- sample(all_individuals, 15)
#' candidates <- setdiff(all_individuals, test_set)
#' train_set <- sample(candidates, 20)
#' 
#' # Compute A-optimality (minimizes average prediction variance)
#' aopt_value <- a_optimality(train_set, test_set, PC_subset)
#' print(paste("A-optimality value:", round(aopt_value, 6)))
#' 
#' # Compare different training set sizes
#' small_train <- sample(candidates, 10)
#' large_train <- sample(candidates, 30)
#' 
#' aopt_small <- a_optimality(small_train, test_set, PC_subset)
#' aopt_large <- a_optimality(large_train, test_set, PC_subset)
#' 
#' print(paste("Small training (n=10):", round(aopt_small, 6)))
#' print(paste("Large training (n=30):", round(aopt_large, 6)))
#' print(paste("Improvement factor:", round(aopt_small / aopt_large, 2)))
#' 
#' @export
a_optimality <- function(train, test = NULL, P, lambda = 1e-6, C = NULL) {
  # Input validation
  validate_matrix_params(P, train = train, test = test, lambda = lambda, C = C)
  
  p_train <- P[train, , drop = FALSE]
  p_test <- if (!is.null(test)) P[test, , drop = FALSE] else NULL
  core <- compute_prediction_core(p_train, p_test, lambda, C)
  
  if (!is.null(core$contrast_core)) {
    sum(diag(core$contrast_core))
  } else {
    sum(diag(core$inv_xtx_reg))
  }
}

#' D-optimality criterion for experimental design
#' 
#' Computes the D-optimality criterion, which minimizes the determinant of the 
#' inverse covariance matrix. Lower values indicate better designs.
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and variables as columns
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @return D-optimality value (negative log determinant)
#' 
#' @details
#' D-optimality minimizes the generalized variance of parameter estimates.
#' It is equivalent to maximizing the determinant of X'X where X is the design matrix.
#' 
#' Mathematical formula: D = log(det(X'X))
#' 
#' @references
#' Kiefer, J. (1959). Optimum experimental designs. Journal of the Royal Statistical Society B.
#' 
#' @examples
#' # Load wheat genomic data
#' data(WheatData)
#' 
#' # Create a subset for demonstration
#' set.seed(456)
#' subset_indices <- sample(1:nrow(Wheat.M), 60)
#' M_subset <- Wheat.M[subset_indices, 1:15]
#' 
#' # Extract principal components for D-optimality
#' pca_result <- prcomp(M_subset, center = TRUE, scale. = TRUE)
#' PC_subset <- pca_result$x[, 1:4]
#' rownames(PC_subset) <- rownames(M_subset)
#' 
#' # Define sets
#' all_individuals <- rownames(PC_subset)
#' test_set <- sample(all_individuals, 10)
#' candidates <- setdiff(all_individuals, test_set)
#' train_set <- sample(candidates, 15)
#' 
#' # Compute D-optimality (minimizes generalized variance)
#' dopt_value <- d_optimality(train_set, test_set, PC_subset)
#' print(paste("D-optimality value:", round(dopt_value, 4)))
#' 
#' # Compare with different training sets
#' train_set2 <- sample(candidates, 15)
#' dopt_value2 <- d_optimality(train_set2, test_set, PC_subset)
#' 
#' print(paste("Training set 1 D-opt:", round(dopt_value, 4)))
#' print(paste("Training set 2 D-opt:", round(dopt_value2, 4)))
#' 
#' if (dopt_value < dopt_value2) {
#'   print("Training set 1 is better (lower D-optimality)")
#' } else {
#'   print("Training set 2 is better (lower D-optimality)")
#' }
#' 
#' @export
d_optimality <- function(train, test = NULL, P, lambda = 1e-6, C = NULL) {
  p_train <- P[train, , drop = FALSE]
  p_test <- if (!is.null(test)) P[test, , drop = FALSE] else NULL
  core <- compute_prediction_core(p_train, p_test, lambda, C)
  
  if (!is.null(core$contrast_core)) {
    log_det <- determinant(core$contrast_core, logarithm = TRUE)$modulus
  } else {
    log_det <- determinant(core$inv_xtx_reg, logarithm = TRUE)$modulus
  }
  -as.numeric(log_det)
}

#' E-optimality criterion for experimental design
#' 
#' Computes the E-optimality criterion, which minimizes the maximum eigenvalue of the 
#' inverse covariance matrix. Lower values indicate better designs.
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and variables as columns
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @return E-optimality value (negative log of minimum eigenvalue)
#' 
#' @details
#' E-optimality minimizes the maximum variance of any linear combination of parameters.
#' It is equivalent to maximizing the minimum eigenvalue of X'X where X is the design matrix.
#' 
#' Mathematical formula: E = log(λ_max((X'X)^(-1)))
#' 
#' @references
#' Atkinson, A.C., Donev, A.N., Tobias, R.D. (2007). Optimum Experimental Designs. Oxford.
#' 
#' @examples
#' # Load wheat genomic data
#' data(WheatData)
#' 
#' # Create a subset for demonstration
#' set.seed(321)
#' subset_indices <- sample(1:nrow(Wheat.M), 80)
#' M_subset <- Wheat.M[subset_indices, 1:20]
#' 
#' # Extract principal components
#' pca_result <- prcomp(M_subset, center = TRUE, scale. = TRUE)
#' PC_subset <- pca_result$x[, 1:5]
#' rownames(PC_subset) <- rownames(M_subset)
#' 
#' # Define sets
#' all_individuals <- rownames(PC_subset)
#' test_set <- sample(all_individuals, 15)
#' candidates <- setdiff(all_individuals, test_set)
#' train_set <- sample(candidates, 20)
#' 
#' # Compute E-optimality (minimizes maximum variance)
#' eopt_value <- e_optimality(train_set, test_set, PC_subset)
#' print(paste("E-optimality value:", round(eopt_value, 6)))
#' 
#' # Compare with A and D optimality for same training set
#' aopt_value <- a_optimality(train_set, test_set, PC_subset)
#' dopt_value <- d_optimality(train_set, test_set, PC_subset)
#' 
#' print(paste("A-optimality:", round(aopt_value, 6)))
#' print(paste("D-optimality:", round(dopt_value, 6)))
#' print(paste("E-optimality:", round(eopt_value, 6)))
#' 
#' @export
e_optimality <- function(train, test = NULL, P, lambda = 1e-6, C = NULL) {
  p_train <- P[train, , drop = FALSE]
  p_test <- if (!is.null(test)) P[test, , drop = FALSE] else NULL
  core <- compute_prediction_core(p_train, p_test, lambda, C)
  
  if (!is.null(core$contrast_core)) {
    eigenvals <- eigen(core$contrast_core, symmetric = TRUE, only.values = TRUE)$values
  } else {
    eigenvals <- eigen(core$inv_xtx_reg, symmetric = TRUE, only.values = TRUE)$values
  }
  -log(min(eigenvals))
}

# =============================================================================
# PREDICTION ERROR VARIANCE CRITERIA (LITERATURE-CORRECTED)
# =============================================================================

#' Mean Prediction Error Variance (PEV) - Literature Corrected
#' 
#' Computes the mean prediction error variance for the test set given a training set.
#' This measures the average uncertainty in predictions.
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and variables as columns
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @param normalized Whether to normalize by trace for scale-invariance (default: FALSE)
#' @return Mean PEV value (lower is better)
#' 
#' @details
#' Corrected prediction error variance formula:
#' PEV = Var(y - ŷ) = σ²[I + X_test(X_train'X_train)^(-1)X_test']
#' 
#' The identity matrix accounts for the inherent variability of new observations.
#' This is the correct formula according to experimental design literature.
#' 
#' @examples
#' # Load wheat genomic data
#' data(WheatData)
#' 
#' # Create a subset for demonstration
#' set.seed(789)
#' subset_indices <- sample(1:nrow(Wheat.M), 80)
#' M_subset <- Wheat.M[subset_indices, 1:25]
#' 
#' # Extract principal components
#' pca_result <- prcomp(M_subset, center = TRUE, scale. = TRUE)
#' PC_subset <- pca_result$x[, 1:5]
#' rownames(PC_subset) <- rownames(M_subset)
#' 
#' # Define sets
#' all_individuals <- rownames(PC_subset)
#' test_set <- sample(all_individuals, 15)
#' candidates <- setdiff(all_individuals, test_set)
#' train_set <- sample(candidates, 20)
#' 
#' # Compute mean PEV (prediction accuracy measure)
#' pev_value <- pev_mean(train_set, test_set, PC_subset)
#' print(paste("Mean PEV:", round(pev_value, 6)))
#' 
#' # Compare normalized vs non-normalized
#' pev_norm <- pev_mean(train_set, test_set, PC_subset, normalized = TRUE)
#' print(paste("Normalized PEV:", round(pev_norm, 6)))
#' 
#' # Effect of training set size on PEV
#' small_train <- sample(candidates, 10)
#' large_train <- sample(candidates, 30)
#' 
#' pev_small <- pev_mean(small_train, test_set, PC_subset)
#' pev_large <- pev_mean(large_train, test_set, PC_subset)
#' 
#' print(paste("Small training PEV:", round(pev_small, 6)))
#' print(paste("Large training PEV:", round(pev_large, 6)))
#' print(paste("PEV reduction factor:", round(pev_small / pev_large, 2)))
#' 
#' @export
pev_mean <- function(train, test = NULL, P, lambda = 1e-6, C = NULL, normalized = FALSE) {
  p_train <- P[train, , drop = FALSE]
  p_test <- if (!is.null(test)) P[test, , drop = FALSE] else NULL
  core <- compute_prediction_core(p_train, p_test, lambda, C)
  
  # Extract diagonal elements (prediction error variances)
  pev_values <- diag(core$pev_matrix)
  
  # Apply normalization if requested
  if (normalized) {
    trace_val <- sum(pev_values)
    pev_values <- if (trace_val > 0) pev_values / trace_val else pev_values
  }
  
  mean(pev_values)
}

#' Maximum Prediction Error Variance (PEV)
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and variables as columns
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @param normalized Whether to normalize by trace for scale-invariance (default: FALSE)
#' @return Maximum PEV value (lower is better)
#' 
#' @examples
#' # Load wheat genomic data
#' data(WheatData)
#' 
#' # Create a subset for demonstration
#' set.seed(654)
#' subset_indices <- sample(1:nrow(Wheat.M), 80)
#' M_subset <- Wheat.M[subset_indices, 1:25]
#' 
#' # Extract principal components
#' pca_result <- prcomp(M_subset, center = TRUE, scale. = TRUE)
#' PC_subset <- pca_result$x[, 1:5]
#' rownames(PC_subset) <- rownames(M_subset)
#' 
#' # Define sets
#' all_individuals <- rownames(PC_subset)
#' test_set <- sample(all_individuals, 15)
#' candidates <- setdiff(all_individuals, test_set)
#' train_set <- sample(candidates, 20)
#' 
#' # Compute maximum PEV (worst-case prediction uncertainty)
#' pev_max_value <- pev_max(train_set, test_set, PC_subset)
#' pev_mean_value <- pev_mean(train_set, test_set, PC_subset)
#' 
#' print(paste("Maximum PEV:", round(pev_max_value, 6)))
#' print(paste("Mean PEV:", round(pev_mean_value, 6)))
#' print(paste("Max/Mean ratio:", round(pev_max_value / pev_mean_value, 2)))
#' 
#' # Compare normalized versions
#' pev_max_norm <- pev_max(train_set, test_set, PC_subset, normalized = TRUE)
#' pev_mean_norm <- pev_mean(train_set, test_set, PC_subset, normalized = TRUE)
#' 
#' print(paste("Normalized Max PEV:", round(pev_max_norm, 6)))
#' print(paste("Normalized Mean PEV:", round(pev_mean_norm, 6)))
#' 
#' @export
pev_max <- function(train, test = NULL, P, lambda = 1e-6, C = NULL, normalized = FALSE) {
  p_train <- P[train, , drop = FALSE]
  p_test <- if (!is.null(test)) P[test, , drop = FALSE] else NULL
  core <- compute_prediction_core(p_train, p_test, lambda, C)
  
  # Extract diagonal elements (prediction error variances)
  pev_values <- diag(core$pev_matrix)
  
  # Apply normalization if requested
  if (normalized) {
    trace_val <- sum(pev_values)
    pev_values <- if (trace_val > 0) pev_values / trace_val else pev_values
  }
  
  max(pev_values)
}

# =============================================================================
# COEFFICIENT OF DETERMINATION (R²) - LITERATURE-CORRECTED
# =============================================================================

#' Coefficient of Determination (R²) based on experimental design literature
#' 
#' Computes the coefficient of determination as the proportion of variance explained
#' by the model, using the hat matrix diagonal (leverage values).
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and variables as columns
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @param normalized Whether to normalize by trace for scale-invariance (default: FALSE)
#' @return Mean R² value (proportion of variance explained)
#' 
#' @details
#' Computes R² as the average leverage from the hat matrix, which represents
#' the proportion of variance explained by the model in experimental design.
#' 
#' Formula: R^2 = tr(H)/n where H = X(X'X)^(-1)X' is the hat matrix
#' 
#' This is the standard definition used when response data is not available.
#' 
#' @references
#' Fedorov, V.V. (1972). Theory of Optimal Experiments. Academic Press.
#' Atkinson, A.C., Donev, A.N., Tobias, R.D. (2007). Optimum Experimental Designs.
#' 
#' @export
cd_mean <- function(train, test = NULL, P, lambda = 1e-6, C = NULL, normalized = FALSE) {
  p_train <- P[train, , drop = FALSE]
  p_test <- if (!is.null(test)) P[test, , drop = FALSE] else NULL
  core <- compute_prediction_core(p_train, p_test, lambda, C)
  
  # Compute R² as average leverage from hat matrix
  # R² = proportion of variance explained = tr(H)/n
  
  # Leverage values (diagonal of hat matrix H = X(X'X)^(-1)X')
  leverage <- if (!is.null(core$p_test)) {
    diag(core$p_test %*% core$inv_xtx_reg %*% t(core$p_test))
  } else {
    diag(core$p_train %*% core$inv_xtx_reg %*% t(core$p_train))
  }
  
  # Apply normalization if requested
  if (normalized) {
    trace_val <- sum(leverage)
    leverage <- if (trace_val > 0) leverage / trace_val else leverage
  }
  
  # Return average R² (proportion of variance explained)
  mean(leverage)
}

# =============================================================================
# MIXED MODEL CRITERIA (HENDERSON'S BLUP - LITERATURE-CORRECTED)
# =============================================================================

#' Mean Prediction Error Variance for Mixed Models (Henderson's BLUP)
#' 
#' Computes the mean prediction error variance for mixed models accounting for
#' genetic relationships through a kinship matrix. Uses correct BLUP theory.
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and markers as columns
#' @param K Kinship/relationship matrix (symmetric, positive definite)
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @param Vg Genetic variance component (default: 1)
#' @param Ve Error variance component (default: 1)
#' @return Mean PEV value for mixed models (lower is better)
#' 
#' @details
#' Implements the correct BLUP prediction error variance formula from Henderson (1984):
#' PEV = sigma_u^2[G_22 - G_21 V_11^(-1) G_12 - G_21 V_11^(-1) X_1 (X_1' V_11^(-1) X_1)^(-1) X_1' V_11^(-1) G_12]
#' 
#' where V_11 = sigma_u^2 G_11 + sigma_e^2 I and both correction terms are SUBTRACTED.
#' 
#' @references
#' Henderson, C.R. (1984). Applications of Linear Models in Animal Breeding.
#' Searle, S.R., Casella, G., McCulloch, C.E. (1992). Variance Components.
#' 
#' @examples
#' # Load wheat genomic data
#' data(WheatData)
#' 
#' # Create a subset for mixed model demonstration
#' set.seed(987)
#' subset_indices <- sample(1:nrow(Wheat.M), 60)
#' M_subset <- Wheat.M[subset_indices, 1:50]
#' K_subset <- Wheat.K[subset_indices, subset_indices]
#' 
#' # Define train and test sets
#' all_individuals <- rownames(M_subset)
#' test_set <- sample(all_individuals, 12)
#' train_set <- sample(setdiff(all_individuals, test_set), 18)
#' 
#' # Mixed model PEV with different heritability scenarios
#' h2_values <- c(0.2, 0.5, 0.8)
#' 
#' print("PEV Mean MM - Different Heritabilities:")
#' for (h2 in h2_values) {
#'   var_comps <- h2_to_variances(h2)
#'   pev_mm <- pev_mean_mm(train_set, test_set, M_subset, K_subset, 
#'                         Vg = var_comps$Vg, Ve = var_comps$Ve)
#'   print(paste("h² =", h2, "-> PEV =", round(pev_mm, 6)))
#' }
#' 
#' # Compare with non-mixed model PEV using PCA
#' pca_result <- prcomp(M_subset, center = TRUE, scale. = TRUE)
#' PC_subset <- pca_result$x[, 1:5]
#' rownames(PC_subset) <- rownames(M_subset)
#' pev_regular <- pev_mean(train_set, test_set, PC_subset)
#' 
#' print(paste("Regular PEV (PCA):", round(pev_regular, 6)))
#' print("Mixed model accounts for genetic relationships via kinship matrix")
#' 
#' @export
pev_mean_mm <- function(train, test = NULL, P, K, lambda = 1e-6, C = NULL, 
                        Vg = 1, Ve = 1) {
  
  # Input validation
  validate_matrix_params(P, K = K, train = train, test = test, lambda = lambda, C = C)
  validate_variance_components(Vg, Ve)
  
  # Check memory usage for mixed model computation
  all_individuals <- unique(c(train, test))
  check_memory_usage(length(all_individuals), ncol(P), "eigen")
  
  # Calculate heritability for documentation
  h2 <- Vg / (Vg + Ve)
  total_var <- Vg + Ve
  
  # Get indices
  train_idx <- rownames(P) %in% train
  test_idx <- if (!is.null(test)) {
    rownames(P) %in% test
  } else {
    train_idx
  }
  
  # Extract relevant matrices according to Henderson's notation
  p_train <- P[train_idx, , drop = FALSE]        # X_1 (design matrix)
  k_train <- K[train_idx, train_idx, drop = FALSE]    # G_11 
  k_test_train <- K[test_idx, train_idx, drop = FALSE]  # G_21 
  k_test_test <- K[test_idx, test_idx, drop = FALSE]    # G_22 
  
  n_train <- nrow(p_train)
  
  # Variance matrix for training set: V_11 = sigma_u^2 * G_11 + sigma_e^2 * I
  v_train <- Vg * k_train + Ve * diag(n_train)
  
  # Compute V_11^(-1) with numerical stability
  inv_v_train <- safe_matrix_inverse(v_train)
  
  # Mixed model coefficient matrix: C = X_1' V_11^(-1) X_1 + lambda * I
  c_matrix <- t(p_train) %*% inv_v_train %*% p_train + lambda * diag(ncol(p_train))
  inv_c_matrix <- safe_matrix_inverse(c_matrix)
  
  # CORRECTED Henderson's BLUP PEV formula:
  # PEV = sigma_u^2[G_22 - G_21 V_11^(-1) G_12 - G_21 V_11^(-1) X_1 (X_1' V_11^(-1) X_1)^(-1) X_1' V_11^(-1) G_12]
  # CRITICAL: Both correction terms are SUBTRACTED
  
  # Term 1: sigma_u^2 * G_22 (genetic variance of test individuals)
  term1 <- Vg * k_test_test
  
  # Term 2: sigma_u^2 * G_21 V_11^(-1) G_12 (reduction due to genetic covariance with training)
  term2 <- Vg * k_test_train %*% inv_v_train %*% t(k_test_train)
  
  # Term 3: Fixed effects contribution (reduces PEV further)
  middle_part <- inv_v_train %*% p_train %*% inv_c_matrix %*% t(p_train) %*% inv_v_train
  term3 <- Vg * k_test_train %*% middle_part %*% t(k_test_train)
  
  # CORRECTED: Both terms reduce prediction error variance
  blup_pev <- term1 - term2 - term3
  
  # Check for numerical issues
  pev_diag <- diag(blup_pev)
  if (any(pev_diag < 0)) {
    warning("Negative PEV values detected. This may indicate numerical issues.")
  }
  
  # Return mean of diagonal elements
  mean(pmax(pev_diag, 0))  # Ensure non-negative
}

#' Coefficient of Determination (R²) for Mixed Models
#' 
#' Computes the coefficient of determination for mixed models, measuring
#' the expected prediction accuracy of random effects (breeding values).
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and markers as columns
#' @param K Kinship/relationship matrix (symmetric, positive definite)
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @param Vg Genetic variance component (default: 1)
#' @param Ve Error variance component (default: 1)
#' @return Mean R² value for mixed models (higher indicates better prediction accuracy)
#' 
#' @details
#' The Coefficient of Determination for BLUP of random effects is defined as:
#' R^2 = 1 - PEV/sigma_u^2
#' 
#' This represents the proportion of genetic variance that can be predicted,
#' or the reliability of breeding value prediction.
#' 
#' @references
#' Henderson, C.R. (1984). Applications of Linear Models in Animal Breeding.
#' 
#' @examples
#' # Load wheat genomic data
#' data(WheatData)
#' 
#' # Create a subset for mixed model demonstration
#' set.seed(159)
#' subset_indices <- sample(1:nrow(Wheat.M), 60)
#' M_subset <- Wheat.M[subset_indices, 1:50]
#' K_subset <- Wheat.K[subset_indices, subset_indices]
#' 
#' # Define train and test sets
#' all_individuals <- rownames(M_subset)
#' test_set <- sample(all_individuals, 12)
#' train_set <- sample(setdiff(all_individuals, test_set), 18)
#' 
#' # Mixed model R² (reliability) with different heritability scenarios
#' h2_values <- c(0.2, 0.5, 0.8)
#' 
#' print("R² (Coefficient of Determination) - Different Heritabilities:")
#' for (h2 in h2_values) {
#'   var_comps <- h2_to_variances(h2)
#'   cd_mm <- cd_mean_mm(train_set, test_set, M_subset, K_subset, 
#'                       Vg = var_comps$Vg, Ve = var_comps$Ve)
#'   print(paste("h² =", h2, "-> R² =", round(cd_mm, 4)))
#' }
#' 
#' # Compare training set sizes
#' small_train <- sample(setdiff(all_individuals, test_set), 10)
#' large_train <- sample(setdiff(all_individuals, test_set), 25)
#' 
#' var_comps <- h2_to_variances(0.5)  # Medium heritability
#' cd_small <- cd_mean_mm(small_train, test_set, M_subset, K_subset, 
#'                        Vg = var_comps$Vg, Ve = var_comps$Ve)
#' cd_large <- cd_mean_mm(large_train, test_set, M_subset, K_subset, 
#'                        Vg = var_comps$Vg, Ve = var_comps$Ve)
#' 
#' print(paste("Small training set (n=10): R² =", round(cd_small, 4)))
#' print(paste("Large training set (n=25): R² =", round(cd_large, 4)))
#' print(paste("Improvement with larger set:", round(cd_large - cd_small, 4)))
#' 
#' @export
cd_mean_mm <- function(train, test = NULL, P, K, lambda = 1e-6, C = NULL, 
                       Vg = 1, Ve = 1) {
  
  # Compute PEV using corrected Henderson's formula
  pev_value <- pev_mean_mm(train, test, P, K, lambda, C, Vg, Ve)
  
  # R^2 = 1 - PEV/sigma_u^2 (reliability of breeding value prediction)
  r_squared <- 1 - (pev_value / Vg)
  
  # Ensure R^2 is in valid [0,1] range
  max(0, min(1, r_squared))
}

# =============================================================================
# HERITABILITY-BASED WRAPPER FUNCTIONS
# =============================================================================

#' Mean Prediction Error Variance for Mixed Models using heritability
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and markers as columns
#' @param K Kinship/relationship matrix (symmetric, positive definite)
#' @param h2 Heritability (between 0 and 1)
#' @param lambda Ridge regularization parameter (default: 1e-6)
#' @param C Contrast matrix (default: NULL)
#' @param total_var Total variance (default: 1)
#' @param normalized If TRUE, return PEV as proportion of total variance
#' @return Mean PEV value
#' @export
pev_mean_mm_h2 <- function(train, test = NULL, P, K, h2, lambda = 1e-6, 
                           C = NULL, total_var = 1, normalized = FALSE) {
  
  # Convert heritability to variance components
  var_comp <- h2_to_variances(h2, total_var)
  
  # Compute PEV
  pev <- pev_mean_mm(train, test, P, K, lambda, C, var_comp$Vg, var_comp$Ve)
  
  # Normalize by total variance if requested
  if (normalized) {
    pev / total_var
  } else {
    pev
  }
}

#' Coefficient of Determination for Mixed Models using heritability
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and markers as columns
#' @param K Kinship/relationship matrix (symmetric, positive definite)
#' @param h2 Heritability (between 0 and 1)
#' @param lambda Ridge regularization parameter (default: 1e-6)
#' @param C Contrast matrix (default: NULL)
#' @param total_var Total variance (default: 1)
#' @return Mean R² value
#' @export
cd_mean_mm_h2 <- function(train, test = NULL, P, K, h2, lambda = 1e-6, 
                          C = NULL, total_var = 1) {
  
  # Convert heritability to variance components
  var_comp <- h2_to_variances(h2, total_var)
  
  # Compute R²
  cd_mean_mm(train, test, P, K, lambda, C, var_comp$Vg, var_comp$Ve)
}

# =============================================================================
# LEGACY COMPATIBILITY AND ADDITIONAL CRITERIA
# =============================================================================

#' Legacy influence measure (for backward compatibility)
#' 
#' Computes the leverage-based influence measure that was incorrectly labeled 
#' as "coefficient of determination" in the original implementation.
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and variables as columns
#' @param lambda Ridge regularization parameter (default: 1e-6)
#' @param C Contrast matrix (default: NULL)
#' @param normalized Whether to normalize (default: FALSE)
#' @return Mean influence measure value
#' 
#' @details
#' This function computes what the original cd_mean actually calculated:
#' Influence = (leverage_i / (1 - leverage_i)²) × p
#' 
#' This is related to Cook's distance without the residual component.
#' For proper R², use cd_mean() which computes average leverage.
#' 
#' @export
influence_measure_legacy <- function(train, test = NULL, P, lambda = 1e-6, C = NULL, 
                                    normalized = FALSE) {
  
  p_train <- P[train, , drop = FALSE]
  p_test <- if (!is.null(test)) P[test, , drop = FALSE] else NULL
  core <- compute_prediction_core(p_train, p_test, lambda, C)
  
  # Compute leverage values (diagonal of hat matrix)
  p <- ncol(core$p_train)
  
  leverage <- if (!is.null(core$p_test)) {
    diag(core$p_test %*% core$inv_xtx_reg %*% t(core$p_test))
  } else {
    diag(core$p_train %*% core$inv_xtx_reg %*% t(core$p_train))
  }
  
  # Original (incorrect) "coefficient of determination" formula
  influence_values <- (leverage / (1 - leverage)^2) * p
  
  # Apply normalization if requested
  if (normalized) {
    trace_val <- sum(influence_values)
    influence_values <- if (trace_val > 0) influence_values / trace_val else influence_values
  }
  
  mean(influence_values)
}

#' G-optimality criterion (minimizes maximum prediction variance)
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and variables as columns
#' @param lambda Ridge regularization parameter (default: 1e-6)
#' @return G-optimality value
#' 
#' @details
#' G-optimality minimizes the maximum entry in the diagonal of the hat matrix,
#' providing the best worst-case prediction precision.
#' 
#' @references
#' Kiefer, J. (1975). Construction and optimality of generalized Youden designs.
#' 
#' @examples
#' # Load wheat genomic data
#' data(WheatData)
#' 
#' # Create a subset for demonstration
#' set.seed(753)
#' subset_indices <- sample(1:nrow(Wheat.M), 80)
#' M_subset <- Wheat.M[subset_indices, 1:20]
#' 
#' # Extract principal components
#' pca_result <- prcomp(M_subset, center = TRUE, scale. = TRUE)
#' PC_subset <- pca_result$x[, 1:5]
#' rownames(PC_subset) <- rownames(M_subset)
#' 
#' # Define sets
#' all_individuals <- rownames(PC_subset)
#' test_set <- sample(all_individuals, 15)
#' candidates <- setdiff(all_individuals, test_set)
#' train_set <- sample(candidates, 20)
#' 
#' # Compute G-optimality (minimizes maximum prediction variance)
#' gopt_value <- g_optimality(train_set, test_set, PC_subset)
#' print(paste("G-optimality value:", round(gopt_value, 6)))
#' 
#' # Compare with other optimality criteria
#' aopt <- a_optimality(train_set, test_set, PC_subset)
#' dopt <- d_optimality(train_set, test_set, PC_subset)
#' eopt <- e_optimality(train_set, test_set, PC_subset)
#' 
#' print("Optimality Comparison:")
#' print(paste("A-optimality (average variance):", round(aopt, 6)))
#' print(paste("D-optimality (volume):", round(dopt, 6)))
#' print(paste("E-optimality (max eigenvalue):", round(eopt, 6)))
#' print(paste("G-optimality (max hat diagonal):", round(gopt_value, 6)))
#' 
#' @export
g_optimality <- function(train, test = NULL, P, lambda = 1e-6) {
  
  p_train <- P[train, , drop = FALSE]
  p_test <- if (!is.null(test)) P[test, , drop = FALSE] else p_train
  
  # Compute (X'X)^(-1)
  xtx <- crossprod(p_train)
  xtx_reg <- xtx + lambda * diag(ncol(p_train))
  inv_xtx <- safe_matrix_inverse(xtx_reg, lambda = 0)
  
  # Compute diagonal of hat matrix
  hat_diag <- diag(p_test %*% inv_xtx %*% t(p_test))
  
  # Return maximum (worst-case prediction variance)
  max(hat_diag)
}

#' I-optimality criterion (minimizes average prediction variance)
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and variables as columns
#' @param lambda Ridge regularization parameter (default: 1e-6)
#' @return I-optimality value
#' 
#' @details
#' I-optimality minimizes the average prediction variance over the design space.
#' 
#' @references
#' Fedorov, V.V. (1972). Theory of Optimal Experiments.
#' 
#' @examples
#' # Load wheat genomic data
#' data(WheatData)
#' 
#' # Create a subset for demonstration
#' set.seed(951)
#' subset_indices <- sample(1:nrow(Wheat.M), 80)
#' M_subset <- Wheat.M[subset_indices, 1:20]
#' 
#' # Extract principal components
#' pca_result <- prcomp(M_subset, center = TRUE, scale. = TRUE)
#' PC_subset <- pca_result$x[, 1:5]
#' rownames(PC_subset) <- rownames(M_subset)
#' 
#' # Define sets
#' all_individuals <- rownames(PC_subset)
#' test_set <- sample(all_individuals, 15)
#' candidates <- setdiff(all_individuals, test_set)
#' train_set <- sample(candidates, 20)
#' 
#' # Compute I-optimality (average prediction variance over design space)
#' iopt_value <- i_optimality(train_set, test_set, PC_subset)
#' print(paste("I-optimality value:", round(iopt_value, 6)))
#' 
#' # Compare I-optimality with A-optimality (both are average-based)
#' aopt_value <- a_optimality(train_set, test_set, PC_subset)
#' print(paste("A-optimality (parameter space):", round(aopt_value, 6)))
#' print(paste("I-optimality (design space):", round(iopt_value, 6)))
#' 
#' # Test different training set sizes
#' small_train <- sample(candidates, 12)
#' large_train <- sample(candidates, 28)
#' 
#' iopt_small <- i_optimality(small_train, test_set, PC_subset)
#' iopt_large <- i_optimality(large_train, test_set, PC_subset)
#' 
#' print(paste("Small training (n=12):", round(iopt_small, 6)))
#' print(paste("Large training (n=28):", round(iopt_large, 6)))
#' print(paste("Improvement factor:", round(iopt_small / iopt_large, 2)))
#' 
#' @export
i_optimality <- function(train, test = NULL, P, lambda = 1e-6) {
  
  p_train <- P[train, , drop = FALSE]
  p_test <- if (!is.null(test)) P[test, , drop = FALSE] else p_train
  
  # Compute (X'X)^(-1)
  xtx <- crossprod(p_train)
  xtx_reg <- xtx + lambda * diag(ncol(p_train))
  inv_xtx <- safe_matrix_inverse(xtx_reg, lambda = 0)
  
  # Compute average prediction variance
  pred_var_matrix <- p_test %*% inv_xtx %*% t(p_test)
  
  # Return average prediction variance
  mean(diag(pred_var_matrix))
}

# =============================================================================
# UNIFIED CRITERION FUNCTION
# =============================================================================

#' Unified criterion function for optimization
#' 
#' A unified interface to all optimization criteria. This function automatically
#' dispatches to the appropriate criterion based on the criterion name and parameters.
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and variables as columns
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @param criterion Character string specifying the criterion to compute
#' @param K Kinship matrix for mixed models (required for mixed model criteria)
#' @param Vg Genetic variance component for mixed models (default: 1 if K provided)
#' @param Ve Error variance component for mixed models (default: 1 if K provided)
#' @return Criterion value (interpretation depends on specific criterion)
#' 
#' @details
#' Available criteria:
#' 
#' **Classical Optimality:**
#' - "a_optimality": A-optimality (trace of inverse covariance)
#' - "d_optimality": D-optimality (log determinant)
#' - "e_optimality": E-optimality (minimum eigenvalue)
#' - "g_optimality": G-optimality (maximum prediction variance)
#' - "i_optimality": I-optimality (average prediction variance)
#' 
#' **Prediction-Based:**
#' - "pev_mean": Mean prediction error variance
#' - "pev_mean_normalized": Normalized mean PEV
#' - "pev_max": Maximum prediction error variance
#' - "pev_max_normalized": Normalized maximum PEV
#' - "cd_mean": Coefficient of determination (R²)
#' - "cd_mean_normalized": Normalized R²
#' 
#' **Mixed Models:**
#' - "pev_mean_mm": Mean PEV for mixed models (requires K)
#' - "cd_mean_mm": R² for mixed models (requires K)
#' 
#' **Legacy (backward compatibility):**
#' - "AOPT", "DOPT", "EOPT": Classical optimality criteria
#' - "PEVMEAN", "PEVMEAN2": PEV (normalized version)
#' - "PEVMAX", "PEVMAX2": Maximum PEV (normalized version)
#' - "CDMEAN", "CDMEAN2": R² (normalized version)
#' - "PEVMEANMM", "CDMEANMM": Mixed model versions
#' 
#' @export
criterion <- function(train, test = NULL, P, lambda = 1e-6, C = NULL,
                      criterion = "pev_mean", K = NULL, Vg = NULL, Ve = NULL) {
  
  # Handle mixed model criteria
  if (!is.null(K) || grepl("_mm$", criterion) || criterion %in% c("PEVMEANMM", "CDMEANMM")) {
    if (is.null(Vg)) Vg <- 1
    if (is.null(Ve)) Ve <- 1
    
    return(switch(criterion,
      "pev_mean_mm" = pev_mean_mm(train, test, P, K, lambda, C, Vg, Ve),
      "cd_mean_mm" = cd_mean_mm(train, test, P, K, lambda, C, Vg, Ve),
      "PEVMEANMM" = pev_mean_mm(train, test, P, K, lambda, C, Vg, Ve),
      "CDMEANMM" = cd_mean_mm(train, test, P, K, lambda, C, Vg, Ve),
      stop("Unknown mixed model criterion: ", criterion)
    ))
  }
  
  # Handle standard criteria
  switch(criterion,
    # Classical optimality
    "a_optimality" = a_optimality(train, test, P, lambda, C),
    "d_optimality" = d_optimality(train, test, P, lambda, C),
    "e_optimality" = e_optimality(train, test, P, lambda, C),
    "g_optimality" = g_optimality(train, test, P, lambda),
    "i_optimality" = i_optimality(train, test, P, lambda),
    
    # Prediction-based
    "pev_mean" = pev_mean(train, test, P, lambda, C, normalized = FALSE),
    "pev_mean_normalized" = pev_mean(train, test, P, lambda, C, normalized = TRUE),
    "pev_max" = pev_max(train, test, P, lambda, C, normalized = FALSE),
    "pev_max_normalized" = pev_max(train, test, P, lambda, C, normalized = TRUE),
    "cd_mean" = cd_mean(train, test, P, lambda, C, normalized = FALSE),
    "cd_mean_normalized" = cd_mean(train, test, P, lambda, C, normalized = TRUE),
    
    # Distance-based criteria (P must be a square distance matrix with individual names as row/col names)
    "max_to_test"      = distance_train_to_test_max(train, test, P, lambda, C),
    "mean_to_test"     = distance_train_to_test_mean(train, test, P, lambda, C),
    "neg_min_internal" = distance_internal_min(train, test, P, lambda, C),
    "neg_mean_internal"= distance_internal_mean(train, test, P, lambda, C),

    # Legacy support for old names
    "AOPT" = a_optimality(train, test, P, lambda, C),
    "DOPT" = d_optimality(train, test, P, lambda, C),
    "EOPT" = e_optimality(train, test, P, lambda, C),
    "PEVMEAN" = pev_mean(train, test, P, lambda, C, normalized = FALSE),
    "PEVMEAN2" = pev_mean(train, test, P, lambda, C, normalized = TRUE),
    "PEVMAX" = pev_max(train, test, P, lambda, C, normalized = FALSE),
    "PEVMAX2" = pev_max(train, test, P, lambda, C, normalized = TRUE),
    "CDMEAN" = cd_mean(train, test, P, lambda, C, normalized = FALSE),
    "CDMEAN2" = cd_mean(train, test, P, lambda, C, normalized = TRUE),
    
    stop("Unknown criterion: ", criterion)
  )
}