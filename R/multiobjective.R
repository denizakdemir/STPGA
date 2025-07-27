#' @title STPGA Multi-Objective Genetic Algorithm
#' @description Multi-objective optimization for subset selection
#' @author Deniz Akdemir

#' Multi-objective genetic algorithm for subset selection with test set
#' @param Pcs Principal components matrix (can be NULL)
#' @param Dist Distance matrix (can be NULL)
#' @param Kernel Kernel matrix (can be NULL)
#' @param candidates Vector of candidate individual names
#' @param test Vector of test individual names
#' @param ntoselect Number of individuals to select for training
#' @param criteria Vector of optimization criteria names
#' @param criteria_types Vector indicating minimize/maximize for each criterion
#' @param plot_directions Vector for plotting directions (1 or -1 for each criterion)
#' @param npop Population size (default: 100)
#' @param mutprob Mutation probability (default: 0.8)
#' @param mutintensity Mutation intensity (default: 1)
#' @param niterations Maximum iterations (default: 500)
#' @param lambda Ridge regression parameter (default: 1e-6)
#' @param plot_iterations Whether to plot during iterations (default: FALSE)
#' @param mc.cores Number of cores for parallel processing (default: 1)
#' @param initial_population Initial population (default: NULL)
#' @param C Contrast matrix (default: NULL)
#' @param axes_labels Labels for plot axes (default: NULL)
#' @param archive_size Maximum size of Pareto archive (default: 200)
#' @param diversity_maintenance Whether to maintain diversity (default: TRUE)
#' @return List containing Pareto front and optimization results
subset_ga_multiobjective <- function(Pcs = NULL, Dist = NULL, Kernel = NULL, 
                                    candidates, test, ntoselect, criteria, 
                                    criteria_types, plot_directions,
                                    npop = 100, mutprob = 0.8, mutintensity = 1, 
                                    niterations = 500, lambda = 1e-6, 
                                    plot_iterations = FALSE, mc.cores = 1, 
                                    initial_population = NULL, C = NULL, 
                                    axes_labels = NULL, archive_size = 200,
                                    diversity_maintenance = TRUE) {
  
  # Input validation
  if (is.null(Pcs) && is.null(Dist) && is.null(Kernel)) {
    stop("At least one of Pcs, Dist, or Kernel must be provided")
  }
  
  if (length(criteria) != length(criteria_types)) {
    stop("Length of criteria must match length of criteria_types")
  }
  
  # Determine primary matrix for optimization
  P <- determine_primary_matrix(Pcs, Dist, Kernel)
  
  # Initialize population
  if (is.null(initial_population)) {
    population <- initialize_population(candidates, ntoselect, npop)
  } else {
    population <- initial_population
    # Fill remaining slots if needed
    while (length(population) < npop) {
      population <- c(population, list(sample(candidates, ntoselect, replace = FALSE)))
    }
    # Truncate if too large
    if (length(population) > npop) {
      population <- population[1:npop]
    }
  }
  
  # Initialize archive for Pareto solutions
  pareto_archive <- list()
  archive_fitness <- matrix(NA, 0, length(criteria))
  
  # Track optimization history
  generation_stats <- list()
  convergence_history <- matrix(NA, niterations, length(criteria))
  
  # Main evolution loop
  for (generation in 1:niterations) {
    
    # Evaluate population
    population_fitness <- evaluate_multiobjective_population(
      population, P, test, criteria, lambda, C, mc.cores
    )
    
    # Update Pareto archive
    archive_result <- update_pareto_archive(
      population, population_fitness, pareto_archive, archive_fitness,
      archive_size, criteria_types
    )
    pareto_archive <- archive_result$archive
    archive_fitness <- archive_result$fitness
    
    # Track convergence
    convergence_history[generation, ] <- compute_convergence_metrics(archive_fitness)
    
    # Non-dominated sorting
    dominance_ranks <- non_dominated_sorting(population_fitness, criteria_types)
    
    # Selection for next generation
    selected_indices <- multiobjective_selection(
      population_fitness, dominance_ranks, npop * 0.5
    )
    selected_population <- population[selected_indices]
    
    # Generate offspring
    offspring <- generate_multiobjective_offspring(
      selected_population, candidates, ntoselect, npop - length(selected_population),
      mutprob, mutintensity, mc.cores
    )
    
    # Combine parents and offspring
    combined_population <- c(selected_population, offspring)
    combined_fitness <- evaluate_multiobjective_population(
      combined_population, P, test, criteria, lambda, C, mc.cores
    )
    
    # Environmental selection (keep best npop individuals)
    if (diversity_maintenance) {
      survival_indices <- diversity_preserving_selection(
        combined_population, combined_fitness, npop, criteria_types
      )
    } else {
      combined_ranks <- non_dominated_sorting(combined_fitness, criteria_types)
      survival_indices <- environmental_selection(combined_fitness, combined_ranks, npop)
    }
    
    population <- combined_population[survival_indices]
    
    # Store generation statistics
    generation_stats[[generation]] <- compute_generation_stats(
      population, combined_fitness[survival_indices, ], archive_fitness
    )
    
    # Plotting
    if (plot_iterations && generation %% 10 == 0) {
      plot_pareto_progress(archive_fitness, criteria, plot_directions, 
                          generation, axes_labels)
    }
    
    # Convergence check
    if (generation > 50 && check_convergence(convergence_history, generation)) {
      message(paste("Converged at generation", generation))
      break
    }
  }
  
  # Final results
  final_pareto_front <- extract_pareto_front(pareto_archive, archive_fitness, criteria_types)
  
  list(
    pareto_solutions = final_pareto_front$solutions,
    pareto_fitness = final_pareto_front$fitness,
    all_archive = pareto_archive,
    convergence_history = convergence_history[1:generation, , drop = FALSE],
    generation_stats = generation_stats[1:generation],
    parameters = list(
      npop = npop,
      niterations = generation,
      criteria = criteria,
      criteria_types = criteria_types
    )
  )
}

#' Multi-objective GA without test set (single population)
#' @param Pcs Principal components matrix
#' @param Dist Distance matrix
#' @param Kernel Kernel matrix
#' @param candidates Candidate individuals
#' @param ntoselect Number to select
#' @param criteria Optimization criteria
#' @param criteria_types Minimize/maximize indicators
#' @param plot_directions Plot directions
#' @param ... Additional parameters
#' @return Multi-objective optimization results
subset_ga_multiobjective_single <- function(Pcs = NULL, Dist = NULL, Kernel = NULL,
                                           candidates, ntoselect, criteria,
                                           criteria_types, plot_directions, ...) {
  
  # Call main multiobjective function with test = NULL
  subset_ga_multiobjective(
    Pcs = Pcs, Dist = Dist, Kernel = Kernel,
    candidates = candidates, test = NULL, ntoselect = ntoselect,
    criteria = criteria, criteria_types = criteria_types,
    plot_directions = plot_directions, ...
  )
}

#' Determine primary matrix for optimization
#' @param Pcs Principal components matrix
#' @param Dist Distance matrix
#' @param Kernel Kernel matrix
#' @return Primary matrix to use for optimization
determine_primary_matrix <- function(Pcs, Dist, Kernel) {
  
  if (!is.null(Pcs)) {
    return(Pcs)
  } else if (!is.null(Kernel)) {
    return(Kernel)
  } else if (!is.null(Dist)) {
    return(Dist)
  } else {
    stop("No valid matrix provided")
  }
}

#' Initialize random population
#' @param candidates Candidate individuals
#' @param ntoselect Number to select per solution
#' @param npop Population size
#' @return List of initial solutions
initialize_population <- function(candidates, ntoselect, npop) {
  
  population <- list()
  for (i in 1:npop) {
    population[[i]] <- sample(candidates, ntoselect, replace = FALSE)
  }
  return(population)
}

#' Evaluate population for multiple objectives
#' @param population List of solutions
#' @param P Primary matrix
#' @param test Test set
#' @param criteria Vector of criteria names
#' @param lambda Ridge parameter
#' @param C Contrast matrix
#' @param mc.cores Number of cores
#' @return Matrix of fitness values (solutions x criteria)
evaluate_multiobjective_population <- function(population, P, test, criteria, 
                                              lambda, C, mc.cores) {
  
  n_solutions <- length(population)
  n_criteria <- length(criteria)
  fitness_matrix <- matrix(NA, n_solutions, n_criteria)
  
  # Evaluate each criterion
  for (j in 1:n_criteria) {
    criterion_name <- criteria[j]
    
    # Evaluate population for this criterion
    fitness_values <- evaluate_population(
      population, P, test, criterion_name, lambda, C, mc.cores = mc.cores
    )
    
    fitness_matrix[, j] <- fitness_values
  }
  
  colnames(fitness_matrix) <- criteria
  return(fitness_matrix)
}

#' Update Pareto archive with new solutions
#' @param population Current population
#' @param fitness Fitness matrix
#' @param archive Current archive
#' @param archive_fitness Archive fitness matrix
#' @param max_size Maximum archive size
#' @param criteria_types Optimization types (minimize/maximize)
#' @return Updated archive and fitness
update_pareto_archive <- function(population, fitness, archive, archive_fitness,
                                 max_size, criteria_types) {
  
  # Combine current archive with new solutions
  combined_solutions <- c(archive, population)
  combined_fitness <- rbind(archive_fitness, fitness)
  
  # Find non-dominated solutions
  non_dominated_indices <- find_non_dominated(combined_fitness, criteria_types)
  
  # Update archive
  new_archive <- combined_solutions[non_dominated_indices]
  new_archive_fitness <- combined_fitness[non_dominated_indices, , drop = FALSE]
  
  # Limit archive size if necessary
  if (length(new_archive) > max_size) {
    # Use crowding distance to select most diverse solutions
    selected_indices <- select_diverse_subset(new_archive_fitness, max_size)
    new_archive <- new_archive[selected_indices]
    new_archive_fitness <- new_archive_fitness[selected_indices, , drop = FALSE]
  }
  
  list(archive = new_archive, fitness = new_archive_fitness)
}

#' Non-dominated sorting for multi-objective optimization
#' @param fitness Fitness matrix
#' @param criteria_types Vector indicating minimize (TRUE) or maximize (FALSE)
#' @return Vector of dominance ranks
non_dominated_sorting <- function(fitness, criteria_types) {
  
  n_solutions <- nrow(fitness)
  dominance_count <- numeric(n_solutions)  # Number of solutions that dominate this one
  dominated_list <- vector("list", n_solutions)  # Solutions dominated by this one
  ranks <- numeric(n_solutions)
  
  # Initialize dominated lists
  for (i in 1:n_solutions) {
    dominated_list[[i]] <- integer(0)
  }
  
  # Calculate dominance relationships
  for (i in 1:(n_solutions-1)) {
    for (j in (i+1):n_solutions) {
      dominance_result <- check_dominance(fitness[i, ], fitness[j, ], criteria_types)
      
      if (dominance_result == 1) {  # i dominates j
        dominated_list[[i]] <- c(dominated_list[[i]], j)
        dominance_count[j] <- dominance_count[j] + 1
      } else if (dominance_result == -1) {  # j dominates i
        dominated_list[[j]] <- c(dominated_list[[j]], i)
        dominance_count[i] <- dominance_count[i] + 1
      }
    }
  }
  
  # Assign ranks
  current_rank <- 1
  front <- which(dominance_count == 0)
  
  while (length(front) > 0) {
    ranks[front] <- current_rank
    next_front <- integer(0)
    
    for (i in front) {
      for (j in dominated_list[[i]]) {
        dominance_count[j] <- dominance_count[j] - 1
        if (dominance_count[j] == 0) {
          next_front <- c(next_front, j)
        }
      }
    }
    
    front <- next_front
    current_rank <- current_rank + 1
  }
  
  return(ranks)
}

#' Check dominance between two solutions
#' @param solution1 First solution fitness
#' @param solution2 Second solution fitness
#' @param criteria_types Optimization types
#' @return 1 if solution1 dominates, -1 if solution2 dominates, 0 if non-dominated
check_dominance <- function(solution1, solution2, criteria_types) {
  
  better_count <- 0
  worse_count <- 0
  
  for (i in 1:length(solution1)) {
    if (criteria_types[i]) {  # Minimize
      if (solution1[i] < solution2[i]) {
        better_count <- better_count + 1
      } else if (solution1[i] > solution2[i]) {
        worse_count <- worse_count + 1
      }
    } else {  # Maximize
      if (solution1[i] > solution2[i]) {
        better_count <- better_count + 1
      } else if (solution1[i] < solution2[i]) {
        worse_count <- worse_count + 1
      }
    }
  }
  
  if (better_count > 0 && worse_count == 0) {
    return(1)  # solution1 dominates
  } else if (worse_count > 0 && better_count == 0) {
    return(-1)  # solution2 dominates
  } else {
    return(0)  # Non-dominated
  }
}

#' Find non-dominated solutions
#' @param fitness Fitness matrix
#' @param criteria_types Optimization types
#' @return Indices of non-dominated solutions
find_non_dominated <- function(fitness, criteria_types) {
  
  ranks <- non_dominated_sorting(fitness, criteria_types)
  return(which(ranks == 1))
}

#' Multi-objective selection based on dominance ranks
#' @param fitness Fitness matrix
#' @param ranks Dominance ranks
#' @param n_select Number of solutions to select
#' @return Indices of selected solutions
multiobjective_selection <- function(fitness, ranks, n_select) {
  
  selected_indices <- integer(0)
  current_rank <- 1
  
  while (length(selected_indices) < n_select) {
    rank_indices <- which(ranks == current_rank)
    
    if (length(selected_indices) + length(rank_indices) <= n_select) {
      # Add entire rank
      selected_indices <- c(selected_indices, rank_indices)
    } else {
      # Partial selection based on crowding distance
      remaining_slots <- n_select - length(selected_indices)
      crowding_distances <- calculate_crowding_distance(fitness[rank_indices, , drop = FALSE])
      
      # Select solutions with highest crowding distance
      selected_from_rank <- rank_indices[order(crowding_distances, decreasing = TRUE)[1:remaining_slots]]
      selected_indices <- c(selected_indices, selected_from_rank)
    }
    
    current_rank <- current_rank + 1
    
    # Safety check
    if (current_rank > max(ranks)) break
  }
  
  return(selected_indices)
}

#' Calculate crowding distance for diversity preservation
#' @param fitness Fitness matrix for a single front
#' @return Vector of crowding distances
calculate_crowding_distance <- function(fitness) {
  
  n_solutions <- nrow(fitness)
  n_objectives <- ncol(fitness)
  
  if (n_solutions <= 2) {
    return(rep(Inf, n_solutions))
  }
  
  crowding_distance <- numeric(n_solutions)
  
  for (obj in 1:n_objectives) {
    # Sort by objective value
    sorted_indices <- order(fitness[, obj])
    
    # Set boundary points to infinite distance
    crowding_distance[sorted_indices[1]] <- Inf
    crowding_distance[sorted_indices[n_solutions]] <- Inf
    
    # Calculate crowding distance for intermediate points
    obj_range <- fitness[sorted_indices[n_solutions], obj] - fitness[sorted_indices[1], obj]
    
    if (obj_range > 0) {
      for (i in 2:(n_solutions-1)) {
        crowding_distance[sorted_indices[i]] <- crowding_distance[sorted_indices[i]] +
          (fitness[sorted_indices[i+1], obj] - fitness[sorted_indices[i-1], obj]) / obj_range
      }
    }
  }
  
  return(crowding_distance)
}

#' Generate offspring for multi-objective optimization
#' @param parents Parent solutions
#' @param candidates Candidate pool
#' @param ntoselect Solution size
#' @param n_offspring Number of offspring to generate
#' @param mutprob Mutation probability
#' @param mutintensity Mutation intensity
#' @param mc.cores Number of cores
#' @return List of offspring solutions
generate_multiobjective_offspring <- function(parents, candidates, ntoselect, n_offspring,
                                             mutprob, mutintensity, mc.cores) {
  
  offspring <- list()
  
  generate_single_offspring <- function(i) {
    # Select two parents
    parent_indices <- sample(1:length(parents), 2, replace = TRUE)
    parent1 <- parents[[parent_indices[1]]]
    parent2 <- parents[[parent_indices[2]]]
    
    # Crossover
    child <- crossover(parent1, parent2, ntoselect, candidates, "uniform")
    
    # Mutation
    if (runif(1) < mutprob) {
      child <- mutate_solution(child, candidates, mutintensity, ntoselect)
    }
    
    return(child)
  }
  
  # Generate offspring
  if (mc.cores > 1 && n_offspring > 10) {
    offspring <- mclapply(1:n_offspring, generate_single_offspring,
                         mc.cores = mc.cores, mc.preschedule = FALSE)
  } else {
    offspring <- lapply(1:n_offspring, generate_single_offspring)
  }
  
  return(offspring)
}

#' Plot Pareto front progress
#' @param fitness Fitness matrix
#' @param criteria Criteria names
#' @param directions Plot directions
#' @param generation Current generation
#' @param axes_labels Axis labels
plot_pareto_progress <- function(fitness, criteria, directions, generation, axes_labels) {
  
  if (is.null(fitness) || nrow(fitness) == 0) return()
  
  # Simple 2D plot for first two objectives
  if (ncol(fitness) >= 2) {
    plot(fitness[, 1] * directions[1], fitness[, 2] * directions[2],
         main = paste("Pareto Front - Generation", generation),
         xlab = if (!is.null(axes_labels)) axes_labels[1] else criteria[1],
         ylab = if (!is.null(axes_labels)) axes_labels[2] else criteria[2],
         pch = 16, col = "blue")
  }
}

#' Select diverse subset using crowding distance
#' @param fitness Fitness matrix
#' @param max_size Maximum number of solutions to select
#' @return Indices of selected solutions
select_diverse_subset <- function(fitness, max_size) {
  n_solutions <- nrow(fitness)
  
  if (n_solutions <= max_size) {
    return(1:n_solutions)
  }
  
  # Calculate crowding distances
  crowding_distances <- calculate_crowding_distance(fitness)
  
  # Select solutions with highest crowding distance
  selected_indices <- order(crowding_distances, decreasing = TRUE)[1:max_size]
  
  return(selected_indices)
}

#' Compute convergence metrics for multi-objective optimization
#' @param fitness Archive fitness matrix
#' @return Vector of convergence metrics
compute_convergence_metrics <- function(fitness) {
  if (is.null(fitness) || nrow(fitness) == 0) {
    return(rep(Inf, ncol(fitness)))
  }
  
  # Return range of each objective as convergence metric
  apply(fitness, 2, function(x) max(x) - min(x))
}

#' Extract final Pareto front
#' @param archive Archive of solutions
#' @param fitness Archive fitness matrix
#' @param criteria_types Optimization types
#' @return List with Pareto solutions and fitness
extract_pareto_front <- function(archive, fitness, criteria_types) {
  if (length(archive) == 0) {
    return(list(solutions = list(), fitness = matrix(NA, 0, length(criteria_types))))
  }
  
  # Find non-dominated solutions one more time
  non_dominated_indices <- find_non_dominated(fitness, criteria_types)
  
  list(
    solutions = archive[non_dominated_indices],
    fitness = fitness[non_dominated_indices, , drop = FALSE]
  )
}

#' Check convergence based on convergence history
#' @param convergence_history Matrix of convergence metrics over generations
#' @param current_generation Current generation number
#' @param window_size Window size for checking convergence (default: 20)
#' @param tolerance Convergence tolerance (default: 1e-6)
#' @return TRUE if converged, FALSE otherwise
check_convergence <- function(convergence_history, current_generation, 
                             window_size = 20, tolerance = 1e-6) {
  
  if (current_generation < window_size) {
    return(FALSE)
  }
  
  # Check if convergence metrics have stabilized
  recent_window <- convergence_history[(current_generation - window_size + 1):current_generation, , drop = FALSE]
  
  # Convergence if all objectives show little change
  for (obj in 1:ncol(recent_window)) {
    obj_range <- max(recent_window[, obj]) - min(recent_window[, obj])
    if (obj_range > tolerance) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

#' Compute generation statistics for multi-objective optimization
#' @param population Current population
#' @param fitness Current fitness matrix
#' @param archive_fitness Archive fitness matrix
#' @return List with generation statistics
compute_generation_stats <- function(population, fitness, archive_fitness) {
  
  # Population diversity
  diversity <- calculate_diversity(population)
  
  # Hypervolume approximation (for 2D case)
  hypervolume <- if (ncol(fitness) == 2 && nrow(fitness) > 0) {
    # Simple hypervolume calculation for 2D
    sorted_indices <- order(fitness[, 1])
    sorted_fitness <- fitness[sorted_indices, ]
    
    hv <- 0
    for (i in 1:nrow(sorted_fitness)) {
      if (i == 1) {
        width <- sorted_fitness[i, 1]
      } else {
        width <- sorted_fitness[i, 1] - sorted_fitness[i-1, 1]
      }
      height <- sorted_fitness[i, 2]
      hv <- hv + width * height
    }
    hv
  } else {
    NA  # Hypervolume not computed for higher dimensions
  }
  
  list(
    diversity = diversity,
    hypervolume = hypervolume,
    archive_size = if (!is.null(archive_fitness)) nrow(archive_fitness) else 0,
    population_size = length(population),
    fitness_ranges = apply(fitness, 2, function(x) max(x) - min(x))
  )
}

#' Diversity-preserving selection for environmental selection
#' @param population Combined population
#' @param fitness Combined fitness matrix
#' @param target_size Target population size
#' @param criteria_types Optimization types
#' @return Indices of selected solutions
diversity_preserving_selection <- function(population, fitness, target_size, criteria_types) {
  
  n_solutions <- length(population)
  if (n_solutions <= target_size) {
    return(1:n_solutions)
  }
  
  # Non-dominated sorting
  ranks <- non_dominated_sorting(fitness, criteria_types)
  
  selected_indices <- integer(0)
  current_rank <- 1
  
  while (length(selected_indices) < target_size && current_rank <= max(ranks)) {
    rank_indices <- which(ranks == current_rank)
    
    if (length(selected_indices) + length(rank_indices) <= target_size) {
      # Add entire rank
      selected_indices <- c(selected_indices, rank_indices)
    } else {
      # Select most diverse solutions from this rank
      remaining_slots <- target_size - length(selected_indices)
      
      if (length(rank_indices) > 1) {
        # Use crowding distance for selection
        rank_fitness <- fitness[rank_indices, , drop = FALSE]
        crowding_distances <- calculate_crowding_distance(rank_fitness)
        
        # Select solutions with highest crowding distance
        diverse_indices <- rank_indices[order(crowding_distances, decreasing = TRUE)[1:remaining_slots]]
        selected_indices <- c(selected_indices, diverse_indices)
      } else {
        # Only one solution in this rank
        selected_indices <- c(selected_indices, rank_indices[1])
      }
    }
    
    current_rank <- current_rank + 1
  }
  
  return(selected_indices)
}

#' Environmental selection for NSGA-II
#' @param fitness Combined fitness matrix
#' @param ranks Dominance ranks
#' @param target_size Target population size
#' @return Indices of selected solutions
environmental_selection <- function(fitness, ranks, target_size) {
  
  selected_indices <- integer(0)
  current_rank <- 1
  
  while (length(selected_indices) < target_size && current_rank <= max(ranks)) {
    rank_indices <- which(ranks == current_rank)
    
    if (length(selected_indices) + length(rank_indices) <= target_size) {
      # Add entire rank
      selected_indices <- c(selected_indices, rank_indices)
    } else {
      # Partial selection based on crowding distance
      remaining_slots <- target_size - length(selected_indices)
      
      if (length(rank_indices) > 1) {
        rank_fitness <- fitness[rank_indices, , drop = FALSE]
        crowding_distances <- calculate_crowding_distance(rank_fitness)
        
        # Select solutions with highest crowding distance
        best_crowding <- rank_indices[order(crowding_distances, decreasing = TRUE)[1:remaining_slots]]
        selected_indices <- c(selected_indices, best_crowding)
      } else {
        selected_indices <- c(selected_indices, rank_indices[1])
      }
    }
    
    current_rank <- current_rank + 1
  }
  
  return(selected_indices)
}

# Legacy wrapper functions for backward compatibility
GenAlgForSubsetSelectionMO <- function(...) {
  .Deprecated("subset_ga_multiobjective", package = "STPGA", 
              msg = "GenAlgForSubsetSelectionMO is deprecated. Use subset_ga_multiobjective instead.")
  
  args <- list(...)
  # Convert old parameter names to new ones
  if (!is.null(args$selectionstats)) args$criteria <- args$selectionstats
  if (!is.null(args$selectionstatstypes)) args$criteria_types <- args$selectionstatstypes
  
  do.call(subset_ga_multiobjective, args)
}

GenAlgForSubsetSelectionMONoTest <- function(...) {
  .Deprecated("subset_ga_multiobjective_single", package = "STPGA", 
              msg = "GenAlgForSubsetSelectionMONoTest is deprecated. Use subset_ga_multiobjective_single instead.")
  
  args <- list(...)
  if (!is.null(args$selectionstats)) args$criteria <- args$selectionstats
  if (!is.null(args$selectionstatstypes)) args$criteria_types <- args$selectionstatstypes
  
  do.call(subset_ga_multiobjective_single, args)
}