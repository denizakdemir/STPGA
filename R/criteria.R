#' @title STPGA Optimization Criteria Functions
#' @description Modern, simplified optimization criteria for subset selection
#' @author Deniz Akdemir

#' Compute prediction core for matrix operations
#' 
#' This function computes the core matrix operations used by all optimization criteria.
#' It handles numerical stability through Cholesky decomposition with SVD fallback.
#' 
#' @param p_train Training prediction matrix (n_train x p)
#' @param p_test Test prediction matrix (n_test x p). If NULL, uses training set for validation
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @return List containing:
#'   \item{inv_xtx_reg}{Regularized inverse of X'X matrix}
#'   \item{pev_matrix}{Prediction error variance matrix}
#'   \item{p_train}{Training prediction matrix}
#'   \item{p_test}{Test prediction matrix (or NULL)}
#'   \item{contrast_core}{Contrast-specific computations (if C provided)}
#'   \item{C}{Contrast matrix (if provided)}
#' @examples
#' # Create example data
#' P <- matrix(rnorm(100), 20, 5)
#' rownames(P) <- paste0("ind", 1:20)
#' train_idx <- 1:15
#' test_idx <- 16:20
#' 
#' # Compute core matrices
#' core <- compute_prediction_core(P[train_idx, ], P[test_idx, ])
#' str(core)
#' @export
compute_prediction_core <- function(p_train, p_test = NULL, lambda = 1e-6, C = NULL) {
  
  # Compute XtX with regularization
  xtx <- crossprod(p_train)
  xtx_reg <- xtx + lambda * diag(ncol(p_train))
  
  # Use Cholesky decomposition for numerical stability
  inv_xtx_reg <- tryCatch({
    L <- chol(xtx_reg)
    chol2inv(L)
  }, error = function(e) {
    warning("Cholesky decomposition failed, falling back to SVD")
    svd_result <- svd(xtx_reg)
    svd_result$v %*% diag(1 / svd_result$d) %*% t(svd_result$u)
  })
  
  # Compute prediction variance matrix
  if (!is.null(p_test)) {
    # Test set provided - compute cross-prediction variances
    # Fix matrix dimension issue: use proper prediction variance formula
    pev_matrix <- p_test %*% inv_xtx_reg %*% t(p_test)
  } else {
    # No test set - use training set for internal validation
    pev_matrix <- p_train %*% inv_xtx_reg %*% t(p_train)
  }
  
  # Handle contrast matrix if provided
  contrast_core <- if (!is.null(C)) {
    c_inv_xtx_ct <- C %*% inv_xtx_reg %*% t(C)
    
    # Check if matrix is invertible
    if (kappa(c_inv_xtx_ct) > 1e12) {
      warning("Contrast matrix leads to ill-conditioned system")
    }
    
    tryCatch({
      l_c <- chol(c_inv_xtx_ct)
      chol2inv(l_c)
    }, error = function(e) {
      svd_c <- svd(c_inv_xtx_ct)
      svd_c$v %*% diag(1 / svd_c$d) %*% t(svd_c$u)
    })
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
#' @examples
#' # Create example data
#' P <- matrix(rnorm(200), 40, 5)
#' rownames(P) <- paste0("ind", 1:40)
#' colnames(P) <- paste0("var", 1:5)
#' 
#' # Define training and test sets
#' train <- paste0("ind", 1:20)
#' test <- paste0("ind", 21:30)
#' 
#' # Compute A-optimality
#' aopt_value <- a_optimality(train, test, P)
#' print(paste("A-optimality:", round(aopt_value, 4)))
#' 
#' # Compare different training set sizes
#' aopt_small <- a_optimality(train[1:10], test, P)
#' aopt_large <- a_optimality(train, test, P)
#' print(paste("Small training set:", round(aopt_small, 4)))
#' print(paste("Large training set:", round(aopt_large, 4)))
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
#' @examples
#' # Create example data
#' P <- matrix(rnorm(150), 30, 5)
#' rownames(P) <- paste0("ind", 1:30)
#' 
#' # Define training and test sets
#' train <- paste0("ind", 1:15)
#' test <- paste0("ind", 16:25)
#' 
#' # Compute D-optimality
#' dopt_value <- d_optimality(train, test, P)
#' print(paste("D-optimality:", round(dopt_value, 4)))
#' 
#' # Compare with A-optimality
#' aopt_value <- a_optimality(train, test, P)
#' print(paste("A-optimality:", round(aopt_value, 4)))
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
#' @examples
#' # Create example data
#' P <- matrix(rnorm(120), 24, 5)
#' rownames(P) <- paste0("ind", 1:24)
#' 
#' # Define training set
#' train <- paste0("ind", 1:12)
#' test <- paste0("ind", 13:20)
#' 
#' # Compute E-optimality
#' eopt_value <- e_optimality(train, test, P)
#' print(paste("E-optimality:", round(eopt_value, 4)))
#' 
#' # Compare all three optimality criteria
#' aopt <- a_optimality(train, test, P)
#' dopt <- d_optimality(train, test, P)
#' eopt <- e_optimality(train, test, P)
#' 
#' criteria_comparison <- data.frame(
#'   Criterion = c("A-optimality", "D-optimality", "E-optimality"),
#'   Value = c(aopt, dopt, eopt)
#' )
#' print(criteria_comparison)
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

#' Mean Prediction Error Variance (PEV)
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
#' PEV measures the expected squared prediction error for new observations.
#' It accounts for both model uncertainty and noise. Lower values indicate
#' more precise predictions.
#' 
#' @examples
#' # Create example genomic data
#' set.seed(123)
#' n_individuals <- 50
#' n_markers <- 100
#' P <- matrix(rnorm(n_individuals * n_markers), n_individuals, n_markers)
#' rownames(P) <- paste0("ind", 1:n_individuals)
#' colnames(P) <- paste0("snp", 1:n_markers)
#' 
#' # Define training and test sets
#' train <- paste0("ind", 1:30)
#' test <- paste0("ind", 31:45)
#' 
#' # Compute PEV
#' pev_unnorm <- pev_mean(train, test, P, normalized = FALSE)
#' pev_norm <- pev_mean(train, test, P, normalized = TRUE)
#' 
#' print(paste("PEV (unnormalized):", round(pev_unnorm, 6)))
#' print(paste("PEV (normalized):", round(pev_norm, 6)))
#' 
#' # Compare different training set sizes
#' train_small <- paste0("ind", 1:15)
#' train_large <- paste0("ind", 1:35)
#' 
#' pev_small <- pev_mean(train_small, test, P)
#' pev_large <- pev_mean(train_large, test, P)
#' 
#' print(paste("PEV with small training set:", round(pev_small, 6)))
#' print(paste("PEV with large training set:", round(pev_large, 6)))
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
#' Computes the maximum prediction error variance for the test set given a training set.
#' This measures the worst-case uncertainty in predictions.
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and variables as columns
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @param normalized Whether to normalize by trace for scale-invariance (default: FALSE)
#' @return Maximum PEV value (lower is better)
#' 
#' @details
#' Maximum PEV identifies the individual with the highest prediction uncertainty.
#' This is useful for robust design where you want to minimize worst-case scenarios.
#' 
#' @examples
#' # Create example data
#' set.seed(456)
#' P <- matrix(rnorm(300), 50, 6)
#' rownames(P) <- paste0("ind", 1:50)
#' 
#' train <- paste0("ind", 1:25)
#' test <- paste0("ind", 26:40)
#' 
#' # Compare mean vs maximum PEV
#' pev_mean_val <- pev_mean(train, test, P)
#' pev_max_val <- pev_max(train, test, P)
#' 
#' print(paste("Mean PEV:", round(pev_mean_val, 6)))
#' print(paste("Max PEV:", round(pev_max_val, 6)))
#' print(paste("Ratio (max/mean):", round(pev_max_val/pev_mean_val, 2)))
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

#' Mean Cook's Distance for influential observation detection
#' 
#' Computes the mean Cook's Distance, which measures the influence of observations
#' on model predictions. Higher values indicate more influential observations.
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and variables as columns
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @param normalized Whether to normalize by trace for scale-invariance (default: FALSE)
#' @return Mean Cook's Distance value (higher indicates more influence)
#' 
#' @details
#' Cook's Distance measures how much the fitted values change when an observation
#' is removed. It combines leverage (how far an observation is from the center)
#' with residual size. Values > 1 are traditionally considered influential.
#' 
#' @examples
#' # Create example data with some outliers
#' set.seed(789)
#' P <- matrix(rnorm(200), 40, 5)
#' # Add some outliers
#' P[c(5, 15, 25), ] <- P[c(5, 15, 25), ] + 3
#' rownames(P) <- paste0(\"ind\", 1:40)
#' 
#' train <- paste0(\"ind\", 1:20)
#' test <- paste0(\"ind\", 21:35)
#' 
#' # Compute Cook's Distance
#' cd_mean_val <- cd_mean(train, test, P)
#' cd_norm_val <- cd_mean(train, test, P, normalized = TRUE)
#' 
#' print(paste(\"Mean Cook's Distance:\", round(cd_mean_val, 4)))
#' print(paste(\"Normalized Cook's Distance:\", round(cd_norm_val, 4)))
#' 
#' # Compare with PEV
#' pev_val <- pev_mean(train, test, P)
#' print(paste(\"PEV for comparison:\", round(pev_val, 6)))
#' 
#' @export
cd_mean <- function(train, test = NULL, P, lambda = 1e-6, C = NULL, normalized = FALSE) {
  p_train <- P[train, , drop = FALSE]
  p_test <- if (!is.null(test)) P[test, , drop = FALSE] else NULL
  core <- compute_prediction_core(p_train, p_test, lambda, C)
  
  # Compute Cook's distances
  p <- ncol(core$p_train)
  
  # Leverage values (diagonal of hat matrix)
  leverage <- if (!is.null(core$p_test)) {
    diag(core$p_test %*% core$inv_xtx_reg %*% t(core$p_test))
  } else {
    diag(core$p_train %*% core$inv_xtx_reg %*% t(core$p_train))
  }
  
  # Cook's distance approximation
  cd_values <- (leverage / (1 - leverage)^2) * p
  
  # Apply normalization if requested
  if (normalized) {
    trace_val <- sum(cd_values)
    cd_values <- if (trace_val > 0) cd_values / trace_val else cd_values
  }
  
  mean(cd_values)
}

#' Mean Prediction Error Variance for Mixed Models
#' 
#' Computes the mean prediction error variance for mixed models accounting for
#' genetic relationships through a kinship matrix. Used in genomic selection.
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
#' This function implements BLUP (Best Linear Unbiased Prediction) theory for
#' mixed models. It accounts for both genetic relationships (through K matrix)
#' and fixed effects (through P matrix). The variance components control the
#' relative importance of genetic vs environmental effects.
#' 
#' @examples
#' # Create example genomic data
#' set.seed(101)
#' n <- 30
#' p <- 50
#' P <- matrix(rnorm(n * p), n, p)
#' rownames(P) <- paste0("ind", 1:n)
#' 
#' # Create kinship matrix (simulated)
#' K <- matrix(0.1, n, n) + diag(0.9, n)
#' rownames(K) <- colnames(K) <- rownames(P)
#' 
#' train <- paste0("ind", 1:20)
#' test <- paste0("ind", 21:25)
#' 
#' # Compute mixed model PEV
#' pev_mm <- pev_mean_mm(train, test, P, K)
#' print(paste("Mixed model PEV:", round(pev_mm, 6)))
#' 
#' # Compare with different variance components
#' pev_high_h2 <- pev_mean_mm(train, test, P, K, Vg = 0.8, Ve = 0.2)  # High heritability
#' pev_low_h2 <- pev_mean_mm(train, test, P, K, Vg = 0.2, Ve = 0.8)   # Low heritability
#' 
#' print(paste("High heritability PEV:", round(pev_high_h2, 6)))
#' print(paste("Low heritability PEV:", round(pev_low_h2, 6)))
#' 
#' @export
pev_mean_mm <- function(train, test = NULL, P, K, lambda = 1e-6, C = NULL, 
                        Vg = 1, Ve = 1) {
  
  # Get indices
  train_idx <- rownames(P) %in% train
  test_idx <- if (!is.null(test)) {
    rownames(P) %in% test
  } else {
    train_idx  # Use training set for validation
  }
  
  # Extract relevant matrices
  p_train <- P[train_idx, , drop = FALSE]
  k_train <- K[train_idx, train_idx, drop = FALSE]
  k_test_train <- K[test_idx, train_idx, drop = FALSE]
  
  # Mixed model equations
  n_train <- nrow(p_train)
  
  # Variance matrix for training set
  v_train <- Vg * k_train + Ve * diag(n_train)
  
  # Compute inverse of V_train
  inv_v_train <- tryCatch({
    l_v <- chol(v_train)
    chol2inv(l_v)
  }, error = function(e) {
    svd_v <- svd(v_train)
    svd_v$v %*% diag(1 / svd_v$d) %*% t(svd_v$u)
  })
  
  # Mixed model core computations
  ptvinvp <- t(p_train) %*% inv_v_train %*% p_train + lambda * diag(ncol(p_train))
  
  inv_ptvinvp <- tryCatch({
    l_pvp <- chol(ptvinvp)
    chol2inv(l_pvp)
  }, error = function(e) {
    svd_pvp <- svd(ptvinvp)
    svd_pvp$v %*% diag(1 / svd_pvp$d) %*% t(svd_pvp$u)
  })
  
  # Prediction error variance for mixed models
  k_test_test <- K[test_idx, test_idx, drop = FALSE]
  
  # BLUP prediction variance
  blup_var <- Vg * (k_test_test - k_test_train %*% inv_v_train %*% t(k_test_train) +
                    k_test_train %*% inv_v_train %*% p_train %*% inv_ptvinvp %*% 
                    t(p_train) %*% inv_v_train %*% t(k_test_train))
  
  # Return mean of diagonal elements
  mean(diag(blup_var))
}

#' Mean Cook's Distance for Mixed Models
#' 
#' Computes the mean Cook's Distance for mixed models, measuring the influence
#' of observations on BLUP predictions while accounting for genetic relationships.
#' 
#' @param train Vector of training set individual names
#' @param test Vector of test set individual names (if NULL, uses training set)
#' @param P Prediction matrix with individuals as rows and markers as columns
#' @param K Kinship/relationship matrix (symmetric, positive definite)
#' @param lambda Ridge regularization parameter for numerical stability (default: 1e-6)
#' @param C Contrast matrix for specific linear combinations (default: NULL)
#' @param Vg Genetic variance component (default: 1)
#' @param Ve Error variance component (default: 1)
#' @return Mean Cook's Distance for mixed models (higher indicates more influence)
#' 
#' @details
#' This extends Cook's Distance to mixed models by incorporating the genetic
#' relationship matrix. It identifies individuals whose inclusion/exclusion
#' significantly affects BLUP predictions.
#' 
#' @examples
#' # Create example data
#' set.seed(202)
#' n <- 25
#' P <- matrix(rnorm(n * 40), n, 40)
#' rownames(P) <- paste0("ind", 1:n)
#' 
#' # Kinship matrix with some strong relationships
#' K <- matrix(0.05, n, n)
#' diag(K) <- 1
#' # Add some close relationships
#' K[1:5, 1:5] <- 0.7
#' diag(K[1:5, 1:5]) <- 1
#' rownames(K) <- colnames(K) <- rownames(P)
#' 
#' train <- paste0("ind", 1:15)
#' test <- paste0("ind", 16:20)
#' 
#' # Compute mixed model Cook's Distance
#' cd_mm <- cd_mean_mm(train, test, P, K)
#' print(paste("Mixed model Cook's Distance:", round(cd_mm, 6)))
#' 
#' # Compare with regular Cook's Distance
#' cd_regular <- cd_mean(train, test, P)
#' print(paste("Regular Cook's Distance:", round(cd_regular, 6)))
#' 
#' @export
cd_mean_mm <- function(train, test = NULL, P, K, lambda = 1e-6, C = NULL, 
                       Vg = 1, Ve = 1) {
  
  # Use same computation as pev_mean_mm but normalize by Vg
  pev_value <- pev_mean_mm(train, test, P, K, lambda, C, Vg, Ve)
  
  # For mixed models, CD is typically PEV normalized by genetic variance
  pev_value / Vg
}

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
#' **Modern names:**
#' - "a_optimality": A-optimality (trace of inverse covariance)
#' - "d_optimality": D-optimality (log determinant)
#' - "e_optimality": E-optimality (minimum eigenvalue)
#' - "pev_mean": Mean prediction error variance
#' - "pev_mean_normalized": Normalized mean PEV
#' - "pev_max": Maximum prediction error variance
#' - "pev_max_normalized": Normalized maximum PEV
#' - "cd_mean": Mean Cook's Distance
#' - "cd_mean_normalized": Normalized mean Cook's Distance
#' - "pev_mean_mm": Mean PEV for mixed models (requires K)
#' - "cd_mean_mm": Mean Cook's Distance for mixed models (requires K)
#' 
#' **Legacy names (for backward compatibility):**
#' - "AOPT", "DOPT", "EOPT": Classical optimality criteria
#' - "PEVMEAN", "PEVMEAN2": PEV (normalized version)
#' - "PEVMAX", "PEVMAX2": Maximum PEV (normalized version)
#' - "CDMEAN", "CDMEAN2": Cook's Distance (normalized version)
#' - "PEVMEANMM", "CDMEANMM": Mixed model versions
#' 
#' @examples
#' # Create example data
#' set.seed(123)
#' P <- matrix(rnorm(200), 40, 5)
#' rownames(P) <- paste0("ind", 1:40)
#' 
#' train <- paste0("ind", 1:20)
#' test <- paste0("ind", 21:30)
#' 
#' # Compare different criteria
#' criteria_results <- list(
#'   a_opt = criterion(train, test, P, criterion = "a_optimality"),
#'   d_opt = criterion(train, test, P, criterion = "d_optimality"),
#'   pev_mean = criterion(train, test, P, criterion = "pev_mean"),
#'   pev_norm = criterion(train, test, P, criterion = "pev_mean_normalized"),
#'   cd_mean = criterion(train, test, P, criterion = "cd_mean")
#' )
#' 
#' print("Criteria comparison:")
#' for(i in seq_along(criteria_results)) {
#'   cat(sprintf("%s: %.6f\\n", names(criteria_results)[i], criteria_results[[i]]))
#' }
#' 
#' # Mixed model example
#' K <- matrix(0.1, 40, 40) + diag(0.9, 40)
#' rownames(K) <- colnames(K) <- rownames(P)
#' 
#' pev_mm <- criterion(train, test, P, criterion = "pev_mean_mm", K = K)
#' print(paste("Mixed model PEV:", round(pev_mm, 6)))
#' 
#' # Legacy compatibility
#' aopt_legacy <- criterion(train, test, P, criterion = "AOPT")
#' aopt_modern <- criterion(train, test, P, criterion = "a_optimality")
#' print(paste("Legacy vs Modern A-opt:", aopt_legacy == aopt_modern))
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
    "a_optimality" = a_optimality(train, test, P, lambda, C),
    "d_optimality" = d_optimality(train, test, P, lambda, C),
    "e_optimality" = e_optimality(train, test, P, lambda, C),
    "pev_mean" = pev_mean(train, test, P, lambda, C, normalized = FALSE),
    "pev_mean_normalized" = pev_mean(train, test, P, lambda, C, normalized = TRUE),
    "pev_max" = pev_max(train, test, P, lambda, C, normalized = FALSE),
    "pev_max_normalized" = pev_max(train, test, P, lambda, C, normalized = TRUE),
    "cd_mean" = cd_mean(train, test, P, lambda, C, normalized = FALSE),
    "cd_mean_normalized" = cd_mean(train, test, P, lambda, C, normalized = TRUE),
    # Legacy support for old names (will be removed in future versions)
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