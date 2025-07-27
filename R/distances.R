#' @title STPGA Distance Functions
#' @description Distance-based criteria and utilities for subset selection
#' @author Deniz Akdemir

#' Distance to ideal point calculation
#' @param X Matrix of solutions (rows) by objectives (columns)
#' @param method Distance method: "euclidean", "manhattan", "chebyshev" (default: "euclidean")
#' @param handle_zeros How to handle zero ranges: "warning", "error", "ignore" (default: "warning")
#' @return Vector of distances to ideal point (0,0,...,0)
distance_to_ideal <- function(X, method = "euclidean", handle_zeros = "warning") {
  
  # Input validation
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data.frame")
  }
  
  if (nrow(X) == 0 || ncol(X) == 0) {
    return(numeric(0))
  }
  
  # Convert to matrix if data.frame
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  
  # Efficient normalization using vectorized operations
  col_mins <- apply(X, 2, min)
  col_maxs <- apply(X, 2, max)
  col_ranges <- col_maxs - col_mins
  
  # Handle zero ranges (constant columns)
  zero_range_cols <- col_ranges == 0
  if (any(zero_range_cols)) {
    message_text <- sprintf("Columns %s have zero range and will be set to 0", 
                           paste(which(zero_range_cols), collapse = ", "))
    
    switch(handle_zeros,
      "error" = stop(message_text),
      "warning" = warning(message_text),
      "ignore" = {}  # Do nothing
    )
    
    col_ranges[zero_range_cols] <- 1  # Avoid division by zero
  }
  
  # Vectorized normalization
  Xc <- sweep(X, 2, col_mins, "-")
  Xc <- sweep(Xc, 2, col_ranges, "/")
  
  # Set zero-range columns to 0 (ideal value)
  if (any(zero_range_cols)) {
    Xc[, zero_range_cols] <- 0
  }
  
  # Calculate distances to ideal point using specified method
  distances <- switch(method,
    "euclidean" = sqrt(rowSums(Xc^2)),
    "manhattan" = rowSums(abs(Xc)),
    "chebyshev" = apply(abs(Xc), 1, max),
    stop("Unknown distance method: ", method)
  )
  
  return(distances)
}

#' Maximum distance from training to test set
#' @param train Training set individual names
#' @param test Test set individual names
#' @param distance_matrix Distance matrix
#' @param lambda Ridge parameter (unused, kept for compatibility)
#' @param C Contrast matrix (unused, kept for compatibility)
#' @return Maximum distance value
distance_train_to_test_max <- function(train, test, distance_matrix, lambda = NULL, C = NULL) {
  
  # Input validation
  if (length(train) == 0 || length(test) == 0) {
    return(Inf)  # Infinite distance for empty sets
  }
  
  # Extract submatrix efficiently
  train_idx <- rownames(distance_matrix) %in% train
  test_idx <- colnames(distance_matrix) %in% test
  
  if (!any(train_idx) || !any(test_idx)) {
    warning("No matching individuals found in distance matrix")
    return(Inf)
  }
  
  sub_matrix <- distance_matrix[train_idx, test_idx, drop = FALSE]
  return(max(sub_matrix))
}

#' Mean distance from training to test set
#' @param train Training set individual names
#' @param test Test set individual names
#' @param distance_matrix Distance matrix
#' @param lambda Ridge parameter (unused, kept for compatibility)
#' @param C Contrast matrix (unused, kept for compatibility)
#' @return Mean distance value
distance_train_to_test_mean <- function(train, test, distance_matrix, lambda = NULL, C = NULL) {
  
  # Input validation
  if (length(train) == 0 || length(test) == 0) {
    return(Inf)
  }
  
  train_idx <- rownames(distance_matrix) %in% train
  test_idx <- colnames(distance_matrix) %in% test
  
  if (!any(train_idx) || !any(test_idx)) {
    warning("No matching individuals found in distance matrix")
    return(Inf)
  }
  
  sub_matrix <- distance_matrix[train_idx, test_idx, drop = FALSE]
  return(mean(sub_matrix))
}

#' Negative minimum distance within training set
#' @param train Training set individual names
#' @param test Test set individual names (unused, kept for compatibility)
#' @param distance_matrix Distance matrix
#' @param lambda Ridge parameter (unused, kept for compatibility)
#' @param C Contrast matrix (unused, kept for compatibility)
#' @return Negative minimum internal distance
distance_internal_min <- function(train, test = NULL, distance_matrix, lambda = NULL, C = NULL) {
  
  if (length(train) < 2) {
    return(0)  # No internal distances with < 2 individuals
  }
  
  train_idx <- rownames(distance_matrix) %in% train
  if (sum(train_idx) < 2) {
    warning("Less than 2 training individuals found in distance matrix")
    return(0)
  }
  
  Dt <- distance_matrix[train_idx, train_idx, drop = FALSE]
  
  # Extract lower triangular values efficiently
  lower_tri_values <- Dt[lower.tri(Dt, diag = FALSE)]
  
  if (length(lower_tri_values) == 0) {
    return(0)
  }
  
  return(-min(lower_tri_values))
}

#' Negative mean distance within training set
#' @param train Training set individual names
#' @param test Test set individual names (unused, kept for compatibility)
#' @param distance_matrix Distance matrix
#' @param lambda Ridge parameter (unused, kept for compatibility)
#' @param C Contrast matrix (unused, kept for compatibility)
#' @return Negative mean internal distance
distance_internal_mean <- function(train, test = NULL, distance_matrix, lambda = NULL, C = NULL) {
  
  if (length(train) < 2) {
    return(0)
  }
  
  train_idx <- rownames(distance_matrix) %in% train
  if (sum(train_idx) < 2) {
    warning("Less than 2 training individuals found in distance matrix")
    return(0)
  }
  
  Dt <- distance_matrix[train_idx, train_idx, drop = FALSE]
  
  # Extract lower triangular values efficiently
  lower_tri_values <- Dt[lower.tri(Dt, diag = FALSE)]
  
  if (length(lower_tri_values) == 0) {
    return(0)
  }
  
  return(-mean(lower_tri_values))
}

#' Unified distance criterion function
#' @param train Training set individual names
#' @param test Test set individual names
#' @param distance_matrix Distance matrix
#' @param criterion Type of distance criterion
#' @param lambda Ridge parameter (for compatibility)
#' @param C Contrast matrix (for compatibility)
#' @return Distance criterion value
distance_criterion <- function(train, test, distance_matrix, 
                              criterion = c("max_to_test", "mean_to_test", 
                                          "neg_min_internal", "neg_mean_internal"),
                              lambda = NULL, C = NULL) {
  
  criterion <- match.arg(criterion)
  
  switch(criterion,
    "max_to_test" = distance_train_to_test_max(train, test, distance_matrix, lambda, C),
    "mean_to_test" = distance_train_to_test_mean(train, test, distance_matrix, lambda, C),
    "neg_min_internal" = distance_internal_min(train, test, distance_matrix, lambda, C),
    "neg_mean_internal" = distance_internal_mean(train, test, distance_matrix, lambda, C)
  )
}

#' Calculate population diversity metrics
#' @param solutions List of solutions in the population
#' @param method Distance calculation method: "hamming", "jaccard", "sorensen" (default: "jaccard")
#' @return Distance matrix or diversity summary statistics
population_distances <- function(solutions, method = c("hamming", "jaccard", "sorensen")) {
  
  method <- match.arg(method)
  n <- length(solutions)
  
  if (n <= 1) return(matrix(0, n, n))
  
  # Convert to binary matrix for efficient computation
  all_items <- sort(unique(unlist(solutions)))
  binary_matrix <- matrix(FALSE, n, length(all_items))
  
  for (i in 1:n) {
    binary_matrix[i, ] <- all_items %in% solutions[[i]]
  }
  
  # Vectorized distance computation
  distances <- switch(method,
    "hamming" = {
      # Hamming distance (proportion of differing positions)
      dist_matrix <- matrix(0, n, n)
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          dist_matrix[i,j] <- dist_matrix[j,i] <- 
            sum(binary_matrix[i,] != binary_matrix[j,]) / length(all_items)
        }
      }
      dist_matrix
    },
    
    "jaccard" = {
      # Jaccard distance (1 - intersection/union)
      dist_matrix <- matrix(0, n, n)
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          intersection <- sum(binary_matrix[i,] & binary_matrix[j,])
          union_size <- sum(binary_matrix[i,] | binary_matrix[j,])
          jaccard_sim <- if (union_size > 0) intersection / union_size else 0
          dist_matrix[i,j] <- dist_matrix[j,i] <- 1 - jaccard_sim
        }
      }
      dist_matrix
    },
    
    "sorensen" = {
      # Sørensen-Dice distance
      dist_matrix <- matrix(0, n, n)
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          intersection <- sum(binary_matrix[i,] & binary_matrix[j,])
          sum_sizes <- sum(binary_matrix[i,]) + sum(binary_matrix[j,])
          sorensen_sim <- if (sum_sizes > 0) 2 * intersection / sum_sizes else 0
          dist_matrix[i,j] <- dist_matrix[j,i] <- 1 - sorensen_sim
        }
      }
      dist_matrix
    }
  )
  
  return(distances)
}

#' Population diversity summary
#' @param solutions List of solutions
#' @param method Distance method for calculation
#' @return List with diversity statistics
diversity_summary <- function(solutions, method = "jaccard") {
  
  if (length(solutions) <= 1) {
    return(list(
      mean_distance = 0,
      min_distance = 0,
      max_distance = 0,
      diversity_index = 0
    ))
  }
  
  dist_matrix <- population_distances(solutions, method)
  
  # Extract upper triangular distances (excluding diagonal)
  upper_tri_distances <- dist_matrix[upper.tri(dist_matrix)]
  
  list(
    mean_distance = mean(upper_tri_distances),
    min_distance = min(upper_tri_distances),
    max_distance = max(upper_tri_distances),
    diversity_index = mean(upper_tri_distances),  # Average pairwise distance
    distance_matrix = dist_matrix
  )
}

# Legacy wrapper functions for backward compatibility
disttoideal <- function(X, method = "euclidean", handle_zeros = "warning") {
  distance_to_ideal(X, method, handle_zeros)
}

dist_to_test <- function(train, test, Dst, lambda = NULL, C = NULL) {
  distance_train_to_test_max(train, test, Dst, lambda, C)
}

dist_to_test2 <- function(train, test, Dst, lambda = NULL, C = NULL) {
  distance_train_to_test_mean(train, test, Dst, lambda, C)
}

neg_dist_in_train <- function(train, test = NULL, Dst, lambda = NULL, C = NULL) {
  distance_internal_min(train, test, Dst, lambda, C)
}

neg_dist_in_train2 <- function(train, test = NULL, Dst, lambda = NULL, C = NULL) {
  distance_internal_mean(train, test, Dst, lambda, C)
}

unified_distance_criterion <- function(train, test, Dst, 
                                      criterion = c("max_to_test", "mean_to_test", 
                                                   "neg_min_internal", "neg_mean_internal"),
                                      lambda = NULL, C = NULL) {
  distance_criterion(train, test, Dst, criterion, lambda, C)
}

calculate_population_distances <- function(solutions, method = c("hamming", "jaccard", "sorensen")) {
  population_distances(solutions, method)
}