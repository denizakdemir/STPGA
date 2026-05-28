#' @title STPGA Matrix Operations
#' @description Matrix utilities for genetic relationship computations
#' @author Deniz Akdemir

#' Compute A matrix in pieces for large marker datasets
#' @param M Marker matrix coded as -1,0,1 (individuals x markers)
#' @param pieces Number of pieces to split computation (default: 10)
#' @param mc.cores Number of cores for parallel processing (default: 1)
#' @return Combined genomic relationship matrix (A matrix)
compute_amatrix <- function(M, pieces = 10, mc.cores = 1) {
  
  # Input validation
  if (missing(M) || !is.matrix(M)) {
    stop("M must be a matrix")
  }
  if (missing(pieces)) {
    pieces <- min(10, ncol(M))
  }
  if (pieces <= 0 || pieces > ncol(M)) {
    stop("pieces must be positive and <= ncol(M)")
  }
  
  # Helper function for computing A matrix pieces
  compute_amatrix_piece <- function(M_subset) {
    # Vectorized computation - more efficient than apply
    pvec <- colMeans(M_subset + 1) / 2
    
    # More efficient matrix operations - fix dimensions
    n_individuals <- nrow(M_subset)
    ones_vector <- rep(1, n_individuals)
    
    # Center the markers: M_centered = (M + 1) - 2 * 1 * p'
    M_centered <- (M_subset + 1) - 2 * outer(ones_vector, pvec)
    MMt <- tcrossprod(M_centered)
    
    return(list(allele_freq = matrix(pvec, ncol = 1), cross_product = MMt))
  }
  
  # Helper function for combining A matrix pieces
  combine_amatrix_pieces <- function(pieces_list) {
    n_pieces <- length(pieces_list)
    
    # Initialize with first piece
    total_cross_product <- pieces_list[[1]]$cross_product
    total_denominator <- sum(2 * pieces_list[[1]]$allele_freq * (1 - pieces_list[[1]]$allele_freq))
    
    # Add remaining pieces
    for (i in 2:n_pieces) {
      total_cross_product <- total_cross_product + pieces_list[[i]]$cross_product
      total_denominator <- total_denominator + 
        sum(2 * pieces_list[[i]]$allele_freq * (1 - pieces_list[[i]]$allele_freq))
    }
    
    if (!is.finite(total_denominator) || total_denominator <= .Machine$double.eps) {
      return(matrix(0, nrow(total_cross_product), ncol(total_cross_product)))
    }
    
    return(total_cross_product / total_denominator)
  }
  
  # Split markers into pieces
  marker_indices <- 1:ncol(M)
  chunk_indices <- split(marker_indices, 
                        factor(sort(rank(marker_indices) %% pieces)))
  
  # Function for parallel processing
  process_chunk <- function(indices) {
    return(compute_amatrix_piece(M[, indices, drop = FALSE]))
  }
  
  # Parallel computation of pieces
  if (mc.cores > 1) {
    pieces_output <- mclapply(chunk_indices, process_chunk, mc.cores = mc.cores)
  } else {
    pieces_output <- lapply(chunk_indices, process_chunk)
  }
  
  # Combine results
  A_matrix <- combine_amatrix_pieces(pieces_output)
  
  # Set row and column names if available
  if (!is.null(rownames(M))) {
    rownames(A_matrix) <- colnames(A_matrix) <- rownames(M)
  }
  
  return(A_matrix)
}

#' Compute genomic relationship matrix using different methods
#' @param M Marker matrix coded as -1,0,1
#' @param method Method for computation: "vanraden", "pieces", "standard" (default: "vanraden")
#' @param pieces Number of pieces for piece-wise computation (default: 10)
#' @param mc.cores Number of cores for parallel processing (default: 1)
#' @return Genomic relationship matrix
genomic_relationship_matrix <- function(M, method = "vanraden", pieces = 10, mc.cores = 1) {
  
  method <- match.arg(method, c("vanraden", "pieces", "standard"))
  
  switch(method,
    "vanraden" = {
      # VanRaden (2008) method
      # Center markers by allele frequencies
      p <- colMeans(M + 1) / 2  # Allele frequencies
      W <- M + 1 - 2 * rep(1, nrow(M)) %o% p
      
      # Compute relationship matrix
      denom <- sum(2 * p * (1 - p))
      G <- if (!is.finite(denom) || denom <= .Machine$double.eps) {
        matrix(0, nrow(M), nrow(M))
      } else {
        tcrossprod(W) / denom
      }
      
      # Set names
      if (!is.null(rownames(M))) {
        rownames(G) <- colnames(G) <- rownames(M)
      }
      
      return(G)
    },
    
    "pieces" = {
      return(compute_amatrix(M, pieces, mc.cores))
    },
    
    "standard" = {
      # Simple centered cross-product
      M_centered <- scale(M, center = TRUE, scale = FALSE)
      G <- tcrossprod(M_centered) / ncol(M)
      
      if (!is.null(rownames(M))) {
        rownames(G) <- colnames(G) <- rownames(M)
      }
      
      return(G)
    }
  )
}

#' Efficient matrix inversion with numerical stability
#' @param X Matrix to invert
#' @param lambda Ridge parameter for regularization (default: 0)
#' @param method Inversion method: "cholesky", "svd", "eigen" (default: "cholesky")
#' @return Inverse matrix
#' @export
safe_matrix_inverse <- function(X, lambda = 0, method = "cholesky") {
  
  method <- match.arg(method, c("cholesky", "svd", "eigen"))
  
  if (!is.matrix(X)) {
    stop("X must be a matrix")
  }
  if (nrow(X) != ncol(X)) {
    stop("X must be a square matrix")
  }
  if (any(!is.finite(X))) {
    stop("X contains non-finite values")
  }
  
  # Add ridge regularization
  if (lambda > 0) {
    X_reg <- X + lambda * diag(nrow(X))
  } else {
    X_reg <- X
  }
  
  switch(method,
    "cholesky" = {
      stability <- matrix_stability_check(X_reg)
      if (isTRUE(stability$rank_deficient) || !isTRUE(stability$is_well_conditioned)) {
        warning("Matrix is rank deficient or ill-conditioned, falling back to SVD")
        return(safe_matrix_inverse(X, lambda, "svd"))
      }
      
      tryCatch({
        L <- chol(X_reg)
        return(chol2inv(L))
      }, error = function(e) {
        warning("Cholesky decomposition failed, falling back to SVD")
        return(safe_matrix_inverse(X, lambda, "svd"))
      })
    },
    
    "svd" = {
      svd_result <- svd(X_reg)
      d <- svd_result$d
      
      # Handle near-zero eigenvalues
      rank_threshold <- max(1e-12, max(d) * 1e-12)
      if (sum(d > rank_threshold) < nrow(X_reg)) {
        warning("Matrix is rank deficient. Using pseudo-inverse.")
      }
      d_inv <- ifelse(d > rank_threshold, 1/d, 0)
      
      return(svd_result$v %*% diag(d_inv) %*% t(svd_result$u))
    },
    
    "eigen" = {
      eigen_result <- eigen(X_reg, symmetric = TRUE)
      values <- eigen_result$values
      vectors <- eigen_result$vectors
      
      # Handle near-zero eigenvalues
      rank_threshold <- max(1e-12, max(abs(values)) * 1e-12)
      if (sum(abs(values) > rank_threshold) < nrow(X_reg)) {
        warning("Matrix is rank deficient. Using pseudo-inverse.")
      }
      values_inv <- ifelse(abs(values) > rank_threshold, 1/values, 0)
      
      return(vectors %*% diag(values_inv) %*% t(vectors))
    }
  )
}

#' Compute condition number and numerical stability metrics
#' @param X Matrix to analyze
#' @return List with stability metrics
matrix_stability_check <- function(X) {
  
  # Compute eigenvalues for symmetric matrices
  if (isSymmetric(X)) {
    eigenvals <- eigen(X, symmetric = TRUE, only.values = TRUE)$values
    scale <- max(abs(eigenvals))
    threshold <- max(1e-12, scale * 1e-12)
    numerical_rank <- sum(abs(eigenvals) > threshold)
    condition_number <- if (numerical_rank < min(dim(X))) {
      Inf
    } else {
      max(abs(eigenvals)) / min(abs(eigenvals))
    }
  } else {
    # Use SVD for non-symmetric matrices
    svd_result <- svd(X)
    singular_vals <- svd_result$d
    scale <- max(singular_vals)
    threshold <- max(1e-12, scale * 1e-12)
    numerical_rank <- sum(singular_vals > threshold)
    condition_number <- if (numerical_rank < min(dim(X))) {
      Inf
    } else {
      max(singular_vals) / min(singular_vals)
    }
  }
  
  # Check for potential issues
  is_positive_definite <- if (isSymmetric(X)) {
    all(eigenvals > threshold)
  } else {
    NA
  }
  
  is_well_conditioned <- is.finite(condition_number) && condition_number < 1e10
  
  list(
    condition_number = condition_number,
    is_positive_definite = is_positive_definite,
    is_well_conditioned = is_well_conditioned,
    rank_deficient = numerical_rank < min(dim(X)),
    numerical_rank = numerical_rank
  )
}

#' Ridge regression with optimal lambda selection
#' @param y Response vector
#' @param X Design matrix
#' @param lambda_seq Sequence of lambda values to test (default: NULL for automatic)
#' @param cv_folds Number of cross-validation folds (default: 5)
#' @param criterion Selection criterion: "mse", "aic", "bic" (default: "mse")
#' @return List with optimal lambda and fitted model
ridge_regression_cv <- function(y, X, lambda_seq = NULL, cv_folds = 5, 
                               criterion = "mse") {
  
  n <- length(y)
  p <- ncol(X)
  
  # Default lambda sequence
  if (is.null(lambda_seq)) {
    lambda_seq <- 10^seq(-8, 2, length.out = 50)
  }
  
  criterion <- match.arg(criterion, c("mse", "aic", "bic"))
  
  if (criterion == "mse" && cv_folds > 1) {
    # Cross-validation
    fold_ids <- sample(rep(1:cv_folds, length.out = n))
    cv_errors <- matrix(NA, cv_folds, length(lambda_seq))
    
    for (fold in 1:cv_folds) {
      train_idx <- fold_ids != fold
      test_idx <- fold_ids == fold
      
      X_train <- X[train_idx, , drop = FALSE]
      y_train <- y[train_idx]
      X_test <- X[test_idx, , drop = FALSE]
      y_test <- y[test_idx]
      
      for (i in seq_along(lambda_seq)) {
        lambda <- lambda_seq[i]
        
        # Fit ridge regression
        XtX_reg <- crossprod(X_train) + lambda * diag(p)
        Xty <- crossprod(X_train, y_train)
        
        tryCatch({
          beta <- solve(XtX_reg, Xty)
          y_pred <- X_test %*% beta
          cv_errors[fold, i] <- mean((y_test - y_pred)^2)
        }, error = function(e) {
          cv_errors[fold, i] <- Inf
        })
      }
    }
    
    # Select optimal lambda
    mean_cv_errors <- colMeans(cv_errors, na.rm = TRUE)
    optimal_lambda <- lambda_seq[which.min(mean_cv_errors)]
    
  } else {
    # Information criteria
    ic_values <- numeric(length(lambda_seq))
    
    for (i in seq_along(lambda_seq)) {
      lambda <- lambda_seq[i]
      
      # Fit model
      XtX_reg <- crossprod(X) + lambda * diag(p)
      Xty <- crossprod(X, y)
      
      tryCatch({
        beta <- solve(XtX_reg, Xty)
        y_pred <- X %*% beta
        mse <- mean((y - y_pred)^2)
        
        # Effective degrees of freedom for ridge regression
        H <- X %*% solve(XtX_reg) %*% t(X)
        df <- sum(diag(H))
        
        ic_values[i] <- switch(criterion,
          "aic" = n * log(mse) + 2 * df,
          "bic" = n * log(mse) + log(n) * df,
          "mse" = mse
        )
      }, error = function(e) {
        ic_values[i] <- Inf
      })
    }
    
    optimal_lambda <- lambda_seq[which.min(ic_values)]
  }
  
  # Fit final model with optimal lambda
  XtX_reg <- crossprod(X) + optimal_lambda * diag(p)
  Xty <- crossprod(X, y)
  beta_optimal <- solve(XtX_reg, Xty)
  
  list(
    lambda_optimal = optimal_lambda,
    beta = beta_optimal,
    lambda_seq = lambda_seq,
    cv_errors = if (exists("mean_cv_errors")) mean_cv_errors else ic_values
  )
}

# Legacy wrapper for backward compatibility
Amat.pieces <- function(M, pieces = 10, mc.cores = 1) {
  compute_amatrix(M, pieces, mc.cores)
}
