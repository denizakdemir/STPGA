#' @title STPGA Genetic Algorithm Operators
#' @description Crossover, mutation, and selection operators for genetic algorithms
#' @author Deniz Akdemir

#' Crossover operation for subset selection
#' @param parent1 First parent solution
#' @param parent2 Second parent solution  
#' @param ntoselect Target size of offspring
#' @param candidates Pool of candidate individuals
#' @param method Crossover method: "uniform", "twopoint", "proportional" (default: "uniform")
#' @param bias Bias toward parent1 (0.5 = equal, >0.5 favors parent1) (default: 0.5)
#' @return Offspring solution
crossover <- function(parent1, parent2, ntoselect, candidates, 
                     method = "uniform", bias = 0.5) {
  
  # Input validation
  if (length(candidates) < ntoselect) {
    stop("candidates must contain at least ntoselect individuals")
  }
  
  if (length(parent1) == 0 || length(parent2) == 0) {
    return(sample(candidates, ntoselect, replace = FALSE))
  }
  
  # Get union and intersection efficiently
  union_parents <- union(parent1, parent2)
  intersect_parents <- intersect(parent1, parent2)
  unique_p1 <- setdiff(parent1, parent2)
  unique_p2 <- setdiff(parent2, parent1)
  
  offspring <- switch(method,
    "uniform" = {
      # Uniform crossover - each gene has equal probability from either parent
      result <- intersect_parents
      
      # Add from unique elements based on bias
      remaining_slots <- ntoselect - length(result)
      if (remaining_slots > 0) {
        unique_all <- c(unique_p1, unique_p2)
        if (length(unique_all) >= remaining_slots) {
          # Bias selection toward parent1
          p1_select <- round(remaining_slots * bias)
          p2_select <- remaining_slots - p1_select
          
          selected <- c()
          if (p1_select > 0 && length(unique_p1) > 0) {
            selected <- c(selected, sample(unique_p1, min(p1_select, length(unique_p1))))
          }
          if (p2_select > 0 && length(unique_p2) > 0) {
            remaining_p2 <- remaining_slots - length(selected)
            if (remaining_p2 > 0) {
              selected <- c(selected, sample(unique_p2, min(remaining_p2, length(unique_p2))))
            }
          }
          result <- c(result, selected)
        } else {
          result <- c(result, unique_all)
        }
      }
      result
    },
    
    "twopoint" = {
      # Two-point crossover for ordered solutions
      all_candidates <- union_parents
      if (length(all_candidates) < ntoselect) {
        c(all_candidates, sample(setdiff(candidates, all_candidates), 
                                ntoselect - length(all_candidates)))
      } else {
        # Create ordered representation and perform two-point crossover
        ordered_candidates <- sample(all_candidates)
        cutpoint1 <- sample(1:max(1, length(ordered_candidates)-1), 1)
        cutpoint2 <- sample(cutpoint1:length(ordered_candidates), 1)
        
        # Select middle segment from parent1, sides from parent2
        middle_from_p1 <- intersect(ordered_candidates[cutpoint1:cutpoint2], parent1)
        sides_from_p2 <- setdiff(parent2, middle_from_p1)
        
        result <- c(middle_from_p1, sides_from_p2)
        if (length(result) > ntoselect) {
          sample(result, ntoselect)
        } else if (length(result) < ntoselect) {
          c(result, sample(setdiff(candidates, result), ntoselect - length(result)))
        } else {
          result
        }
      }
    },
    
    "proportional" = {
      # Proportional selection based on parent sizes
      prop1 <- length(parent1) / (length(parent1) + length(parent2))
      n_from_p1 <- round(ntoselect * prop1)
      n_from_p2 <- ntoselect - n_from_p1
      
      # Sample from each parent proportionally
      selected_p1 <- if (n_from_p1 > 0 && length(parent1) > 0) {
        sample(parent1, min(n_from_p1, length(parent1)))
      } else character(0)
      
      selected_p2 <- if (n_from_p2 > 0 && length(parent2) > 0) {
        remaining_needed <- ntoselect - length(selected_p1)
        available_p2 <- setdiff(parent2, selected_p1)
        if (length(available_p2) > 0) {
          sample(available_p2, min(remaining_needed, length(available_p2)))
        } else character(0)
      } else character(0)
      
      result <- c(selected_p1, selected_p2)
      
      # Fill remaining slots if needed
      if (length(result) < ntoselect) {
        remaining_candidates <- setdiff(candidates, result)
        if (length(remaining_candidates) > 0) {
          result <- c(result, sample(remaining_candidates, 
                                   min(ntoselect - length(result), 
                                       length(remaining_candidates))))
        }
      }
      result
    }
  )
  
  # Ensure exact size and valid candidates
  if (length(offspring) > ntoselect) {
    offspring <- sample(offspring, ntoselect)
  } else if (length(offspring) < ntoselect) {
    remaining <- setdiff(candidates, offspring)
    if (length(remaining) > 0) {
      offspring <- c(offspring, sample(remaining, 
                                     min(ntoselect - length(offspring), 
                                         length(remaining))))
    }
  }
  
  return(offspring)
}

#' Generate offspring from elite solutions
#' @param elites List of elite solutions
#' @param candidates Pool of candidate individuals  
#' @param npop Population size to generate
#' @param mutprob Mutation probability
#' @param mc.cores Number of cores for parallel processing (default: 1)
#' @param mutintensity Mutation intensity (default: 1)
#' @param tabu_memory Tabu memory for avoiding revisits (default: NULL)
#' @param ntoselect Target size of solutions (default: NULL, inferred from elites)
#' @param crossover_method Method for crossover operations (default: "uniform")
#' @param parallel_threshold Minimum population size for parallel processing (default: 50)
#' @return List of offspring solutions
generate_offspring <- function(elites, candidates, npop, mutprob, 
                               mc.cores = 1, mutintensity = 1, 
                               tabu_memory = NULL, ntoselect = NULL,
                               crossover_method = "uniform",
                               parallel_threshold = 50) {
  
  nelite <- length(elites)
  if (nelite == 0) {
    stop("No elite solutions provided")
  }
  
  # Infer ntoselect if not provided
  if (is.null(ntoselect)) {
    ntoselect <- length(elites[[1]])
  }
  
  # Check if solution exists in tabu memory
  check_tabu <- function(x) {
    if (is.null(tabu_memory)) return(FALSE)
    
    for (memory_set in tabu_memory) {
      if (!is.null(memory_set)) {
        for (solution in memory_set) {
          if (identical(sort(x), sort(solution))) {
            return(TRUE)
          }
        }
      }
    }
    return(FALSE)
  }
  
  # Generate single offspring function
  generate_single <- function(i) {
    max_attempts <- 10
    attempt <- 0
    
    repeat {
      attempt <- attempt + 1
      
      # Select parents (with replacement to allow self-crossing for small elite sets)
      parent_indices <- sample(1:nelite, 2, replace = TRUE)
      parent1 <- elites[[parent_indices[1]]]
      parent2 <- elites[[parent_indices[2]]]
      
      # Generate offspring through crossover
      offspring <- crossover(parent1, parent2, ntoselect, candidates, crossover_method)
      
      # Apply mutation
      if (runif(1) < mutprob) {
        offspring <- mutate_solution(offspring, candidates, mutintensity, ntoselect)
      }
      
      # Check tabu memory
      if (!check_tabu(offspring) || attempt >= max_attempts) {
        break
      }
      
      # If in tabu memory, try more radical mutation
      if (attempt < max_attempts) {
        offspring <- mutate_solution(offspring, candidates, 
                                   min(ceiling(ntoselect * 0.3), length(offspring)), 
                                   ntoselect)
      }
    }
    
    return(offspring)
  }
  
  # Generate population (parallel or sequential based on size)
  if (npop >= parallel_threshold && mc.cores > 1) {
    offspring_list <- mclapply(1:npop, generate_single, 
                              mc.cores = mc.cores, mc.preschedule = FALSE)
  } else {
    offspring_list <- lapply(1:npop, generate_single)
  }
  
  return(offspring_list)
}

#' Mutate a solution
#' @param solution Current solution to mutate
#' @param candidates Pool of candidates
#' @param intensity Mutation intensity
#' @param ntoselect Target solution size
#' @return Mutated solution
mutate_solution <- function(solution, candidates, intensity, ntoselect) {
  # Adaptive mutation intensity
  n_mutations <- if (intensity <= 1) {
    rbinom(1, ntoselect, intensity)
  } else {
    min(intensity, floor(ntoselect * 0.5))
  }
  
  if (n_mutations > 0) {
    # Select positions to mutate
    positions_to_mutate <- sample(1:ntoselect, n_mutations, replace = FALSE)
    
    for (pos in positions_to_mutate) {
      # Select new candidate not already in solution
      available_candidates <- setdiff(candidates, solution)
      if (length(available_candidates) > 0) {
        solution[pos] <- sample(available_candidates, 1)
      }
    }
  }
  
  return(solution)
}

#' Calculate population diversity
#' @param population List of solutions in the population
#' @return Diversity measure (0 to 1, higher = more diverse)
calculate_diversity <- function(population) {
  if (length(population) <= 1) return(0)
  
  # Calculate pairwise Jaccard distances
  n_pop <- length(population)
  distances <- numeric(0)
  
  for (i in 1:(n_pop-1)) {
    for (j in (i+1):n_pop) {
      intersection <- length(intersect(population[[i]], population[[j]]))
      union_size <- length(union(population[[i]], population[[j]]))
      jaccard_dist <- 1 - (intersection / union_size)
      distances <- c(distances, jaccard_dist)
    }
  }
  
  return(mean(distances))
}

#' Adaptive mutation rate based on population diversity
#' @param base_mutprob Base mutation probability
#' @param diversity Current population diversity (0-1)
#' @param generation Current generation number
#' @param max_generations Maximum generations
#' @param min_rate Minimum mutation rate (default: 0.01)
#' @param max_rate Maximum mutation rate (default: 0.95)
#' @return Adjusted mutation probability
adaptive_mutation_rate <- function(base_mutprob, diversity, generation, max_generations,
                                  min_rate = 0.01, max_rate = 0.95) {
  
  # Increase mutation when diversity is low
  diversity_factor <- (1 - diversity) * 0.5
  
  # Increase mutation early and late in evolution
  generation_progress <- generation / max_generations
  generation_factor <- if (generation_progress < 0.2 || generation_progress > 0.8) {
    0.2
  } else {
    0
  }
  
  adjusted_rate <- base_mutprob + diversity_factor + generation_factor
  
  # Ensure bounds
  return(pmax(min_rate, pmin(max_rate, adjusted_rate)))
}

#' Tournament selection
#' @param population List of solutions
#' @param fitness Vector of fitness values
#' @param tournament_size Size of tournament (default: 3)
#' @param n_select Number of individuals to select (default: 1)
#' @return Indices of selected individuals
tournament_selection <- function(population, fitness, tournament_size = 3, n_select = 1) {
  selected <- integer(n_select)
  n_pop <- length(population)
  
  for (i in 1:n_select) {
    # Select random individuals for tournament
    tournament_indices <- sample(1:n_pop, tournament_size, replace = FALSE)
    tournament_fitness <- fitness[tournament_indices]
    
    # Select best from tournament (assuming minimization)
    winner_idx <- tournament_indices[which.min(tournament_fitness)]
    selected[i] <- winner_idx
  }
  
  return(selected)
}

#' Elite selection
#' @param population List of solutions
#' @param fitness Vector of fitness values
#' @param n_select Number of elites to select
#' @return Indices of selected elite individuals
elite_selection <- function(population, fitness, n_select) {
  # Return indices of best solutions (assuming minimization)
  order(fitness)[1:n_select]
}

#' Rank-based selection with selection pressure control
#' @param population List of solutions
#' @param fitness Vector of fitness values
#' @param n_select Number of individuals to select
#' @param selection_pressure Selection pressure parameter (1.0-2.0, default: 1.5)
#' @param method Rank selection method: "linear", "exponential" (default: "linear")
#' @return Indices of selected individuals
#' @export
#' 
#' @examples
#' # Create a sample population and fitness values
#' set.seed(123)
#' population <- list(
#'   c("ind1", "ind2", "ind3"),
#'   c("ind2", "ind4", "ind5"),
#'   c("ind1", "ind3", "ind6"),
#'   c("ind4", "ind5", "ind6"),
#'   c("ind1", "ind4", "ind7")
#' )
#' 
#' # Fitness values (lower is better for minimization)
#' fitness <- c(10.5, 8.2, 12.1, 9.8, 7.5)
#' 
#' # Rank-based selection with moderate pressure
#' selected_indices <- rank_selection(population, fitness, n_select = 3, 
#'                                   selection_pressure = 1.4, method = "linear")
#' print("Selected indices:")
#' print(selected_indices)
#' print("Corresponding fitness values:")
#' print(fitness[selected_indices])
#' 
#' # Compare different selection pressures
#' low_pressure <- rank_selection(population, fitness, 3, 1.1, "linear")
#' high_pressure <- rank_selection(population, fitness, 3, 1.8, "linear")
#' 
#' print("Low pressure selection:")
#' print(fitness[low_pressure])
#' print("High pressure selection:")
#' print(fitness[high_pressure])
#' 
#' # Exponential ranking
#' exp_selection <- rank_selection(population, fitness, 3, 1.3, "exponential")
#' print("Exponential ranking selection:")
#' print(fitness[exp_selection])
rank_selection <- function(population, fitness, n_select, selection_pressure = 1.5, 
                          method = "linear") {
  
  # Validate selection pressure
  if (selection_pressure < 1.0 || selection_pressure > 2.0) {
    warning("Selection pressure should be between 1.0 and 2.0. Adjusting to valid range.")
    selection_pressure <- pmax(1.0, pmin(2.0, selection_pressure))
  }
  
  n_pop <- length(population)
  if (n_select > n_pop) {
    warning("Cannot select more individuals than population size")
    n_select <- n_pop
  }
  
  # Rank individuals (1 = best, n_pop = worst for minimization)
  ranks <- rank(fitness, ties.method = "random")
  
  # Calculate selection probabilities based on rank
  probabilities <- switch(method,
    "linear" = {
      # Linear ranking: P(i) = (2 - SP + 2 * (SP - 1) * (rank - 1) / (n - 1)) / n
      linear_probs <- (2 - selection_pressure + 
                      2 * (selection_pressure - 1) * (n_pop - ranks) / (n_pop - 1)) / n_pop
      pmax(linear_probs, 0)  # Ensure non-negative
    },
    
    "exponential" = {
      # Exponential ranking: P(i) = c * r^(rank-1)
      c_value <- (selection_pressure - 1) / (selection_pressure^n_pop - 1)
      exp_probs <- c_value * selection_pressure^(n_pop - ranks)
      exp_probs / sum(exp_probs)  # Normalize
    }
  )
  
  # Perform selection with replacement using calculated probabilities
  selected_indices <- sample(1:n_pop, n_select, replace = TRUE, prob = probabilities)
  
  return(selected_indices)
}

#' Validate and adjust selection parameters
#' @param selection_method Selection method name
#' @param tournament_size Tournament size (for tournament selection)
#' @param selection_pressure Selection pressure (for rank selection)
#' @param population_size Population size
#' @return List with validated parameters and warnings
#' @export
#' 
#' @examples
#' # Test parameter validation with valid parameters
#' validation_good <- validate_selection_parameters(
#'   selection_method = "rank",
#'   tournament_size = 3,
#'   selection_pressure = 1.5,
#'   population_size = 50
#' )
#' 
#' print("Valid parameters:")
#' print(paste("Method:", validation_good$selection_method))
#' print(paste("Tournament size:", validation_good$tournament_size))
#' print(paste("Selection pressure:", validation_good$selection_pressure))
#' print(paste("Warnings:", length(validation_good$warnings)))
#' 
#' # Test parameter validation with invalid parameters
#' validation_bad <- validate_selection_parameters(
#'   selection_method = "invalid_method",
#'   tournament_size = 100,  # Too large
#'   selection_pressure = 3.0,  # Out of range
#'   population_size = 25
#' )
#' 
#' print("Invalid parameters corrected:")
#' print(paste("Method corrected to:", validation_bad$selection_method))
#' print(paste("Tournament size corrected to:", validation_bad$tournament_size))
#' print(paste("Selection pressure corrected to:", validation_bad$selection_pressure))
#' print(paste("Number of warnings:", length(validation_bad$warnings)))
#' 
#' for (warning in validation_bad$warnings) {
#'   print(paste("Warning:", warning))
#' }
validate_selection_parameters <- function(selection_method, tournament_size = 3, 
                                        selection_pressure = 1.5, population_size) {
  
  warnings <- character(0)
  
  # Validate selection method
  valid_methods <- c("tournament", "elite", "rank", "hybrid")
  if (!selection_method %in% valid_methods) {
    warning(paste("Invalid selection method:", selection_method, 
                 ". Using 'tournament' instead."))
    selection_method <- "tournament"
    warnings <- c(warnings, "Invalid selection method corrected")
  }
  
  # Validate tournament size
  if (selection_method %in% c("tournament", "hybrid")) {
    if (tournament_size < 2 || tournament_size > population_size) {
      old_tournament_size <- tournament_size
      tournament_size <- max(2, min(tournament_size, population_size))
      warnings <- c(warnings, 
                   paste("Tournament size adjusted from", old_tournament_size, 
                        "to", tournament_size))
    }
  }
  
  # Validate selection pressure
  if (selection_method %in% c("rank", "hybrid")) {
    if (selection_pressure < 1.0 || selection_pressure > 2.0) {
      old_pressure <- selection_pressure
      selection_pressure <- max(1.0, min(2.0, selection_pressure))
      warnings <- c(warnings,
                   paste("Selection pressure adjusted from", old_pressure,
                        "to", selection_pressure))
    }
  }
  
  # Population size warnings
  if (population_size < 10) {
    warnings <- c(warnings, "Very small population size may lead to premature convergence")
  }
  
  list(
    selection_method = selection_method,
    tournament_size = tournament_size,
    selection_pressure = selection_pressure,
    warnings = warnings
  )
}

#' Crowding replacement for diversity preservation
#' @param offspring New offspring solutions
#' @param offspring_fitness Fitness of offspring
#' @param population Current population
#' @param population_fitness Current population fitness
#' @param crowding_factor Number of candidates to compare for replacement
#' @return Updated population after crowding replacement
crowding_replacement <- function(offspring, offspring_fitness, population, 
                                population_fitness, crowding_factor = 3) {
  
  new_population <- population
  new_fitness <- population_fitness
  
  for (i in 1:length(offspring)) {
    # Select random candidates for replacement
    candidates <- sample(1:length(population), 
                        min(crowding_factor, length(population)), 
                        replace = FALSE)
    
    # Find most similar individual among candidates
    similarities <- sapply(candidates, function(j) {
      intersection <- length(intersect(offspring[[i]], population[[j]]))
      union_size <- length(union(offspring[[i]], population[[j]]))
      intersection / union_size  # Jaccard similarity
    })
    
    # Replace most similar if offspring is better
    most_similar_idx <- candidates[which.max(similarities)]
    if (offspring_fitness[i] < population_fitness[most_similar_idx]) {
      new_population[[most_similar_idx]] <- offspring[[i]]
      new_fitness[most_similar_idx] <- offspring_fitness[i]
    }
  }
  
  return(list(population = new_population, fitness = new_fitness))
}

#' Fitness sharing for diversity preservation
#' @param fitness Vector of fitness values
#' @param population List of solutions
#' @param sharing_radius Radius for fitness sharing (default: 0.1)
#' @param alpha Sharing function exponent (default: 1)
#' @return Adjusted fitness values
fitness_sharing <- function(fitness, population, sharing_radius = 0.1, alpha = 1) {
  if (sharing_radius >= 1) {
    return(as.numeric(fitness))
  }
  
  n_pop <- length(population)
  shared_fitness <- numeric(n_pop)
  
  for (i in 1:n_pop) {
    sharing_sum <- 0
    
    for (j in 1:n_pop) {
      # Calculate distance between solutions
      intersection <- length(intersect(population[[i]], population[[j]]))
      min_size <- min(length(population[[i]]), length(population[[j]]))
      distance <- if (min_size > 0) 1 - (intersection / min_size) else 1
      
      # Apply sharing function
      if (distance < sharing_radius) {
        sharing_value <- 1 - (distance / sharing_radius)^alpha
        sharing_sum <- sharing_sum + sharing_value
      }
    }
    
    # Adjust fitness
    shared_fitness[i] <- fitness[i] * max(sharing_sum, 1)
  }
  
  return(shared_fitness)
}
