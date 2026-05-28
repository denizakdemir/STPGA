#' @title STPGA Population Evaluation Functions
#' @description Efficient population evaluation and fitness computation
#' @author Deniz Akdemir

#' Optimized population evaluation with smart caching
#' @param population List of candidate solutions (training sets)
#' @param P Prediction matrix
#' @param test Test set individual names (can be NULL)
#' @param criterion Optimization criterion name (default: "PEVMEAN2")
#' @param lambda Ridge regression parameter (default: 1e-6)
#' @param C Contrast matrix (default: NULL)
#' @param K Kinship matrix for mixed models (default: NULL)
#' @param Vg Genetic variance matrix (default: NULL)
#' @param Ve Error variance matrix (default: NULL)
#' @param mc.cores Number of cores for parallel processing (default: 1)
#' @param use_cache Whether to use computation cache (default: TRUE)
#' @return Vector of fitness values for population
evaluate_population <- function(population, P, test = NULL, criterion = "PEVMEAN2",
                               lambda = 1e-6, C = NULL, K = NULL, 
                               Vg = NULL, Ve = NULL, mc.cores = 1, use_cache = TRUE) {
  
  # Initialize cache if requested
  if (use_cache) {
    cache <- create_evaluation_cache()
  }
  
  # No need to parse criterion - criterion() function handles all cases
  
  # Evaluation function for single solution
  evaluate_single <- function(i) {
    solution <- population[[i]]
    
    # Check cache first
    if (use_cache) {
      cached_result <- cache$get(solution, test, criterion, lambda)
      if (!is.null(cached_result)) {
        return(cached_result)
      }
    }
    
    # Compute fitness
    tryCatch({
      fitness <- criterion(solution, test, P, lambda, C, criterion, K, Vg, Ve)
      
      # Store in cache
      if (use_cache) {
        cache$store(solution, test, criterion, lambda, fitness)
      }
      
      return(fitness)
    }, error = function(e) {
      warning("Error in population evaluation: ", e$message, call. = FALSE)
      return(Inf)
    })
  }
  
  # Parallel or sequential evaluation
  if (mc.cores > 1 && length(population) > 10) {
    fitness_values <- mclapply(1:length(population), evaluate_single, 
                              mc.cores = mc.cores, mc.preschedule = FALSE)
    fitness_values <- unlist(fitness_values)
  } else {
    fitness_values <- sapply(1:length(population), evaluate_single)
  }
  
  return(fitness_values)
}

#' Smart population evaluation with adaptive batch processing
#' @param population List of solutions
#' @param P Prediction matrix
#' @param test Test set names
#' @param criterion Criterion name
#' @param lambda Ridge parameter
#' @param C Contrast matrix
#' @param K Kinship matrix
#' @param Vg Genetic variance
#' @param Ve Error variance
#' @param mc.cores Number of cores
#' @param batch_size Batch size for processing (default: NULL for automatic)
#' @return Vector of fitness values
evaluate_population_smart <- function(population, P, test = NULL, criterion = "PEVMEAN2",
                                     lambda = 1e-6, C = NULL, K = NULL,
                                     Vg = NULL, Ve = NULL, mc.cores = 1, 
                                     batch_size = NULL) {
  
  n_solutions <- length(population)
  
  # Determine optimal batch size
  if (is.null(batch_size)) {
    batch_size <- determine_optimal_batch_size(n_solutions, mc.cores)
  }
  
  # Group similar solutions for batch processing
  solution_groups <- group_similar_solutions(population, batch_size)
  
  # Process each group
  all_fitness <- numeric(n_solutions)
  processed_count <- 0
  
  for (group in solution_groups) {
    group_indices <- group$indices
    group_solutions <- population[group_indices]
    
    # Evaluate group
    group_fitness <- evaluate_solution_batch(group_solutions, P, test, criterion,
                                            lambda, C, K, Vg, Ve, mc.cores)
    
    # Store results
    all_fitness[group_indices] <- group_fitness
    processed_count <- processed_count + length(group_indices)
  }
  
  return(all_fitness)
}

#' Determine optimal batch size for evaluation
#' @param n_solutions Number of solutions
#' @param mc.cores Number of cores
#' @return Optimal batch size
determine_optimal_batch_size <- function(n_solutions, mc.cores) {
  
  # Base batch size on available cores and problem size
  if (mc.cores == 1) {
    # Sequential processing
    return(min(50, n_solutions))
  } else {
    # Parallel processing
    base_batch <- max(10, n_solutions %/% (mc.cores * 2))
    return(min(100, base_batch))
  }
}

#' Group similar solutions for batch processing
#' @param population List of solutions
#' @param target_batch_size Target size for each batch
#' @return List of solution groups
group_similar_solutions <- function(population, target_batch_size) {
  
  n_solutions <- length(population)
  if (n_solutions <= target_batch_size) {
    return(list(list(indices = 1:n_solutions, solutions = population)))
  }
  
  # Simple grouping by solution size
  solution_sizes <- sapply(population, length)
  size_order <- order(solution_sizes)
  
  # Create groups
  groups <- list()
  current_group <- integer(0)
  current_size <- 0
  
  for (i in size_order) {
    current_group <- c(current_group, i)
    current_size <- current_size + 1
    
    if (current_size >= target_batch_size) {
      groups[[length(groups) + 1]] <- list(
        indices = current_group,
        solutions = population[current_group]
      )
      current_group <- integer(0)
      current_size <- 0
    }
  }
  
  # Add remaining solutions to last group
  if (current_size > 0) {
    if (length(groups) > 0) {
      # Add to last group
      last_group <- groups[[length(groups)]]
      groups[[length(groups)]] <- list(
        indices = c(last_group$indices, current_group),
        solutions = population[c(last_group$indices, current_group)]
      )
    } else {
      # Create new group
      groups[[1]] <- list(
        indices = current_group,
        solutions = population[current_group]
      )
    }
  }
  
  return(groups)
}

#' Evaluate a batch of solutions efficiently
#' @param solutions List of solutions in the batch
#' @param P Prediction matrix
#' @param test Test set
#' @param criterion Criterion name
#' @param lambda Ridge parameter
#' @param C Contrast matrix
#' @param K Kinship matrix
#' @param Vg Genetic variance
#' @param Ve Error variance
#' @param mc.cores Number of cores
#' @return Vector of fitness values
evaluate_solution_batch <- function(solutions, P, test, criterion, lambda, C, K, Vg, Ve, mc.cores) {
  
  # Check for identical solutions in batch to avoid redundant computation
  solution_keys <- vapply(lapply(solutions, sort), paste, collapse = "\r", FUN.VALUE = character(1))
  unique_keys <- unique(solution_keys)
  unique_indices <- match(solution_keys, unique_keys)
  unique_solutions <- lapply(unique_keys, function(key) {
    solutions[[which(solution_keys == key)[1]]]
  })
  
  # Evaluate unique solutions only
  unique_fitness <- numeric(length(unique_solutions))
  
  evaluate_unique <- function(i) {
    solution <- unique_solutions[[i]]
    tryCatch({
      criterion(solution, test, P, lambda, C, criterion, K, Vg, Ve)
    }, error = function(e) {
      warning(paste("Error in batch evaluation:", e$message))
      return(Inf)
    })
  }
  
  if (mc.cores > 1 && length(unique_solutions) > 3) {
    unique_fitness <- unlist(mclapply(1:length(unique_solutions), evaluate_unique,
                                     mc.cores = mc.cores, mc.preschedule = FALSE))
  } else {
    unique_fitness <- sapply(1:length(unique_solutions), evaluate_unique)
  }
  
  # Map back to original solutions
  return(unique_fitness[unique_indices])
}

#' Create evaluation cache for storing computed fitness values
#' @return Cache object with get and store methods
create_evaluation_cache <- function() {
  
  cache_data <- new.env(hash = TRUE)
  
  # Generate cache key from solution and parameters
  generate_key <- function(solution, test, criterion, lambda) {
    solution_key <- paste(sort(solution), collapse = "_")
    test_key <- if (is.null(test)) "NULL" else paste(sort(test), collapse = "_")
    paste(solution_key, test_key, criterion, lambda, sep = "|")
  }
  
  list(
    get = function(solution, test, criterion, lambda) {
      key <- generate_key(solution, test, criterion, lambda)
      if (exists(key, envir = cache_data)) {
        return(get(key, envir = cache_data))
      }
      return(NULL)
    },
    
    store = function(solution, test, criterion, lambda, fitness) {
      key <- generate_key(solution, test, criterion, lambda)
      assign(key, fitness, envir = cache_data)
    },
    
    clear = function() {
      rm(list = ls(cache_data), envir = cache_data)
    },
    
    size = function() {
      length(ls(cache_data))
    }
  )
}

#' Fitness transformation and normalization utilities
#' @param fitness Vector of fitness values
#' @param method Transformation method: "linear", "rank", "tournament" (default: "linear")
#' @param scaling_factor Scaling factor for transformations (default: 2)
#' @return Transformed fitness values
transform_fitness <- function(fitness, method = "linear", scaling_factor = 2) {
  
  method <- match.arg(method, c("linear", "rank", "tournament"))
  
  switch(method,
    "linear" = {
      # Linear scaling to positive values
      min_fitness <- min(fitness, na.rm = TRUE)
      if (min_fitness < 0) {
        return(fitness - min_fitness + 1)
      } else {
        return(fitness)
      }
    },
    
    "rank" = {
      # Rank-based transformation
      ranks <- rank(-fitness, ties.method = "random", na.last = "keep")  # Negative for minimization
      return(as.numeric(ranks))
    },
    
    "tournament" = {
      # Tournament pressure scaling
      min_fitness <- min(fitness, na.rm = TRUE)
      max_fitness <- max(fitness, na.rm = TRUE)
      
      if (max_fitness == min_fitness) {
        return(rep(1, length(fitness)))
      }
      
      # Scale to [1, scaling_factor]
      normalized <- (fitness - min_fitness) / (max_fitness - min_fitness)
      return(1 + normalized * (scaling_factor - 1))
    }
  )
}

#' Population statistics computation
#' @param population List of solutions
#' @param fitness Vector of fitness values
#' @return List with population statistics
compute_population_stats <- function(population, fitness) {
  
  valid_indices <- !is.infinite(fitness) & !is.na(fitness)
  valid_fitness <- fitness[valid_indices]
  valid_population <- population[valid_indices]
  
  if (length(valid_fitness) == 0) {
    return(list(
      best_fitness = Inf,
      worst_fitness = Inf,
      mean_fitness = Inf,
      fitness_std = 0,
      diversity = 0,
      convergence_metric = 0,
      valid_solutions = 0,
      invalid_solutions = length(fitness)
    ))
  }
  
  # Fitness statistics
  best_fitness <- min(valid_fitness)
  worst_fitness <- max(valid_fitness)
  mean_fitness <- mean(valid_fitness)
  fitness_std <- sd(valid_fitness)
  
  # Population diversity
  diversity <- calculate_diversity(valid_population)
  
  # Convergence metric (coefficient of variation)
  convergence_metric <- if (mean_fitness > 0) fitness_std / mean_fitness else 0
  
  list(
    best_fitness = best_fitness,
    worst_fitness = worst_fitness,
    mean_fitness = mean_fitness,
    fitness_std = fitness_std,
    diversity = diversity,
    convergence_metric = convergence_metric,
    valid_solutions = length(valid_fitness),
    invalid_solutions = length(fitness) - length(valid_fitness)
  )
}

# Legacy wrapper functions for backward compatibility
evaluate_population_optimized <- function(population, P, Test, errorstat = "PEVMEAN2", 
                                         lambda = 1e-6, C = NULL, mc.cores = 1, 
                                         K = NULL, Vg = NULL, Ve = NULL) {
  evaluate_population(population, P, Test, errorstat, lambda, C, K, Vg, Ve, mc.cores)
}

evaluate_population_legacy <- function(population, errorstat, Test, P, K = NULL, 
                                      lambda = 1e-6, C = NULL, mc.cores = 1,
                                      Vg = NULL, Ve = NULL) {
  evaluate_population(population, P, Test, errorstat, lambda, C, K, Vg, Ve, mc.cores)
}
