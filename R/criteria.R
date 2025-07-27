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
#' @return Inverse matrix
safe_matrix_inverse <- function(X, lambda = 1e-6) {
  # Add ridge regularization if needed
  if (lambda > 0) {
    X_reg <- X + lambda * diag(nrow(X))
  } else {
    X_reg <- X
  }
  
  # Try Cholesky first (fastest for positive definite)
  tryCatch({
    L <- chol(X_reg)
    chol2inv(L)
  }, error = function(e) {
    # Fall back to SVD for numerical stability
    svd_result <- svd(X_reg)
    d_inv <- ifelse(svd_result$d > 1e-12, 1/svd_result$d, 0)
    svd_result$v %*% diag(d_inv) %*% t(svd_result$u)
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
#' @return List containing core matrices and computations
#' @export
compute_prediction_core <- function(p_train, p_test = NULL, lambda = 1e-6, C = NULL) {
  
  # Compute XtX with regularization
  xtx <- crossprod(p_train)
  xtx_reg <- xtx + lambda * diag(ncol(p_train))
  
  # Compute inverse with numerical stability
  inv_xtx_reg <- safe_matrix_inverse(xtx_reg, lambda = 0)  # Already regularized
  
  # CORRECTED prediction error variance computation
  # Based on literature: PEV = Var(y - ŷ) = σ²[I + X(X'X)⁻¹X']
  # The identity matrix accounts for inherent variability of new observations
  if (!is.null(p_test)) {
    # Test set: PEV includes both model uncertainty AND observation variance
    prediction_var <- p_test %*% inv_xtx_reg %*% t(p_test)
    pev_matrix <- diag(nrow(p_test)) + prediction_var
  } else {
    # Training set: leave-one-out cross-validation
    prediction_var <- p_train %*% inv_xtx_reg %*% t(p_train)
    pev_matrix <- diag(nrow(p_train)) + prediction_var
  }
  
  # Handle contrast matrix if provided
  contrast_core <- if (!is.null(C)) {
    c_inv_xtx_ct <- C %*% inv_xtx_reg %*% t(C)
    
    if (kappa(c_inv_xtx_ct) > 1e12) {
      warning("Contrast matrix leads to ill-conditioned system")
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
    C = C
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
#' @export
a_optimality <- function(train, test = NULL, P, lambda = 1e-6, C = NULL) {
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
#' Formula: R² = tr(H)/n where H = X(X'X)⁻¹X' is the hat matrix
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
  
  # Leverage values (diagonal of hat matrix H = X(X'X)⁻¹X')
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
#' PEV = σ²ᵤ[G₂₂ - G₂₁V₁₁⁻¹G₁₂ - G₂₁V₁₁⁻¹X₁(X₁'V₁₁⁻¹X₁)⁻¹X₁'V₁₁⁻¹G₁₂]
#' 
#' where V₁₁ = σ²ᵤG₁₁ + σ²ₑI and both correction terms are SUBTRACTED.
#' 
#' @references
#' Henderson, C.R. (1984). Applications of Linear Models in Animal Breeding.
#' Searle, S.R., Casella, G., McCulloch, C.E. (1992). Variance Components.
#' 
#' @export
pev_mean_mm <- function(train, test = NULL, P, K, lambda = 1e-6, C = NULL, 
                        Vg = 1, Ve = 1) {
  
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
  p_train <- P[train_idx, , drop = FALSE]        # X₁ (design matrix)
  k_train <- K[train_idx, train_idx, drop = FALSE]    # G₁₁ 
  k_test_train <- K[test_idx, train_idx, drop = FALSE]  # G₂₁ 
  k_test_test <- K[test_idx, test_idx, drop = FALSE]    # G₂₂ 
  
  n_train <- nrow(p_train)
  
  # Variance matrix for training set: V₁₁ = σ²ᵤG₁₁ + σ²ₑI
  v_train <- Vg * k_train + Ve * diag(n_train)
  
  # Compute V₁₁⁻¹ with numerical stability
  inv_v_train <- safe_matrix_inverse(v_train)
  
  # Mixed model coefficient matrix: C = X₁'V₁₁⁻¹X₁ + λI
  c_matrix <- t(p_train) %*% inv_v_train %*% p_train + lambda * diag(ncol(p_train))
  inv_c_matrix <- safe_matrix_inverse(c_matrix)
  
  # CORRECTED Henderson's BLUP PEV formula:
  # PEV = σ²ᵤ[G₂₂ - G₂₁V₁₁⁻¹G₁₂ - G₂₁V₁₁⁻¹X₁(X₁'V₁₁⁻¹X₁)⁻¹X₁'V₁₁⁻¹G₁₂]
  # CRITICAL: Both correction terms are SUBTRACTED
  
  # Term 1: σ²ᵤG₂₂ (genetic variance of test individuals)
  term1 <- Vg * k_test_test
  
  # Term 2: σ²ᵤG₂₁V₁₁⁻¹G₁₂ (reduction due to genetic covariance with training)
  term2 <- Vg * k_test_train %*% inv_v_train %*% t(k_test_train)
  
  # Term 3: Fixed effects contribution (reduces PEV further)
  middle_part <- inv_v_train %*% p_train %*% inv_c_matrix %*% t(p_train) %*% inv_v_train
  term3 <- Vg * k_test_train %*% middle_part %*% t(k_test_train)
  
  # CORRECTED: Both terms reduce prediction error variance
  blup_pev <- term1 - term2 - term3
  
  # Return mean of diagonal elements
  mean(diag(blup_pev))
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
#' R² = 1 - PEV/σ²ᵤ
#' 
#' This represents the proportion of genetic variance that can be predicted,
#' or the reliability of breeding value prediction.
#' 
#' @references
#' Henderson, C.R. (1984). Applications of Linear Models in Animal Breeding.
#' 
#' @export
cd_mean_mm <- function(train, test = NULL, P, K, lambda = 1e-6, C = NULL, 
                       Vg = 1, Ve = 1) {
  
  # Compute PEV using corrected Henderson's formula
  pev_value <- pev_mean_mm(train, test, P, K, lambda, C, Vg, Ve)
  
  # R² = 1 - PEV/σ²ᵤ (reliability of breeding value prediction)
  r_squared <- 1 - (pev_value / Vg)
  
  # Ensure R² is in valid [0,1] range
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
#' @export
g_optimality <- function(train, test = NULL, P, lambda = 1e-6) {
  
  p_train <- P[train, , drop = FALSE]
  p_test <- if (!is.null(test)) P[test, , drop = FALSE] else p_train
  
  # Compute (X'X)⁻¹
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
#' @export
i_optimality <- function(train, test = NULL, P, lambda = 1e-6) {
  
  p_train <- P[train, , drop = FALSE]
  p_test <- if (!is.null(test)) P[test, , drop = FALSE] else p_train
  
  # Compute (X'X)⁻¹
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