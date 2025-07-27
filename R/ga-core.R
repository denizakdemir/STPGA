#' @title STPGA Core Genetic Algorithm Functions
#' @description Clean, well-organized genetic algorithm functions for subset selection
#' @author Deniz Akdemir

#' Primary genetic algorithm for subset selection with training/test split
#' @param P Prediction matrix (individuals x markers/predictors)
#' @param Candidates Vector of candidate individual names
#' @param Test Vector of test individual names
#' @param ntoselect Number of individuals to select for training
#' @param npop Population size (default: 100)
#' @param nelite Number of elite individuals to keep (default: 5)
#' @param keepbest Whether to keep best solution across generations (default: TRUE)
#' @param tabu Whether to use tabu search (default: TRUE)
#' @param tabumemsize Tabu memory size (default: 1)
#' @param mutprob Mutation probability (default: 0.8)
#' @param mutintensity Mutation intensity (default: 1)
#' @param niterations Maximum iterations (default: 500)
#' @param minitbefstop Minimum iterations before stopping (default: 200)
#' @param niterreg Iterations for convergence check (default: 5)
#' @param lambda Ridge regression parameter (default: 1e-6)
#' @param plotiters Whether to plot during iterations (default: FALSE)
#' @param plottype Type of plot (default: 1)
#' @param criterion Optimization criterion (default: "PEVMEAN2")
#' @param C Contrast matrix (default: NULL)
#' @param mc.cores Number of cores for parallel processing (default: 1)
#' @param InitPop Initial population (default: NULL)
#' @param tolconv Convergence tolerance (default: 1e-7)
#' @param Vg Genetic variance matrix for mixed models (default: NULL)
#' @param Ve Error variance matrix for mixed models (default: NULL)
#' @param Fedorov Whether to use Fedorov exchange (default: FALSE)
#' @param adaptive_mutation Whether to use adaptive mutation rates (default: TRUE)
#' @param selection_method Selection method: "tournament", "elite", "hybrid" (default: "tournament")
#' @param tournament_size Tournament size for tournament selection (default: 3)
#' @param diversity_preservation Whether to preserve diversity (default: TRUE)
#' @param diversity_method Diversity method: "crowding", "sharing", "both" (default: "crowding")
#' @param crowding_factor Factor for crowding replacement (default: 3)
#' @param sharing_radius Radius for fitness sharing (default: 0.1)
#' @param crossover_method Crossover method (default: "adaptive")
#' @param diversity_target Target diversity level (default: 0.3)
#' @param verbose Whether to print verbose output (default: FALSE)
#' @return List containing best solution, fitness history, and statistics
subset_ga <- function(P, Candidates, Test, ntoselect, 
                     npop = 100, nelite = 5, keepbest = TRUE, tabu = TRUE, 
                     tabumemsize = 1, mutprob = 0.8, mutintensity = 1, 
                     niterations = 500, minitbefstop = 200, niterreg = 5,
                     lambda = 1e-6, plotiters = FALSE, plottype = 1, 
                     criterion = "PEVMEAN2", C = NULL, mc.cores = 1, 
                     InitPop = NULL, tolconv = 1e-7, Vg = NULL, Ve = NULL, 
                     Fedorov = FALSE, adaptive_mutation = TRUE,
                     selection_method = "tournament", tournament_size = 3,
                     diversity_preservation = TRUE, diversity_method = "crowding",
                     crowding_factor = 3, sharing_radius = 0.1,
                     crossover_method = "adaptive", diversity_target = 0.3,
                     verbose = FALSE) {
  
  # Input validation
  if ((ncol(P) + 1) > ntoselect) {
    warning("The algorithm does not work well with p>ntrain, perhaps use unsupervised dimension reduction on P.")
  }
  
  # Validate inputs
  validate_matrix_params(P, train = Candidates, test = Test, lambda = lambda)
  
  # Initialize population
  if (!is.null(InitPop)) {
    population <- InitPop
    if (length(population) < npop) {
      # Fill remaining slots with random solutions
      additional_pop <- replicate(npop - length(population), 
                                 sample(Candidates, ntoselect), 
                                 simplify = FALSE)
      population <- c(population, additional_pop)
    } else if (length(population) > npop) {
      population <- population[1:npop]
    }
  } else {
    population <- replicate(npop, sample(Candidates, ntoselect), simplify = FALSE)
  }
  
  # Initialize tracking variables
  fitness_history <- matrix(NA, niterations, 3)
  colnames(fitness_history) <- c("best", "mean", "worst")
  
  generation_stats <- list()
  tabu_memory <- if (tabu) vector("list", tabumemsize) else NULL
  
  best_ever_solution <- NULL
  best_ever_fitness <- Inf
  
  no_improvement_count <- 0
  convergence_history <- numeric(niterreg)
  
  # Create progress bar if verbose
  if (verbose) {
    pb <- create_progress_bar(niterations, "GA Evolution")
  }
  
  # Main evolution loop
  for (generation in 1:niterations) {
    
    # Evaluate population
    fitness_values <- evaluate_population(
      population, P, Test, criterion, lambda, C, 
      K = NULL, Vg = Vg, Ve = Ve, mc.cores = mc.cores
    )
    
    # Track best solution
    current_best_idx <- which.min(fitness_values)
    current_best_fitness <- fitness_values[current_best_idx]
    current_best_solution <- population[[current_best_idx]]
    
    # Update global best
    if (current_best_fitness < best_ever_fitness) {
      best_ever_fitness <- current_best_fitness
      best_ever_solution <- current_best_solution
      no_improvement_count <- 0
    } else {
      no_improvement_count <- no_improvement_count + 1
    }
    
    # Record fitness statistics
    fitness_history[generation, "best"] <- current_best_fitness
    fitness_history[generation, "mean"] <- mean(fitness_values)
    fitness_history[generation, "worst"] <- max(fitness_values)
    
    # Population statistics
    generation_stats[[generation]] <- compute_population_stats(population, fitness_values)
    
    # Adaptive mutation rate
    if (adaptive_mutation) {
      diversity <- generation_stats[[generation]]$diversity
      mutprob <- adaptive_mutation_rate(mutprob, diversity, generation, 
                                       niterations, min_rate = 0.01, max_rate = 0.95)
    }
    
    # Convergence check
    convergence_history[((generation - 1) %% niterreg) + 1] <- current_best_fitness
    
    if (generation >= minitbefstop && generation > niterreg) {
      recent_improvement <- max(convergence_history) - min(convergence_history)
      if (recent_improvement < tolconv) {
        if (verbose) {
          message(paste("Convergence achieved at generation", generation))
        }
        break
      }
    }
    
    # Early stopping if no improvement
    if (no_improvement_count > minitbefstop) {
      if (verbose) {
        message(paste("Early stopping at generation", generation, "- no improvement"))
      }
      break
    }
    
    # Selection
    if (selection_method == "tournament") {
      selected_indices <- tournament_selection(population, fitness_values, 
                                              tournament_size, nelite)
    } else if (selection_method == "elite") {
      selected_indices <- elite_selection(population, fitness_values, nelite)
    } else { # hybrid
      elite_indices <- elite_selection(population, fitness_values, nelite %/% 2)
      tournament_indices <- tournament_selection(population, fitness_values, 
                                                tournament_size, nelite - length(elite_indices))
      selected_indices <- c(elite_indices, tournament_indices)
    }
    
    elites <- population[selected_indices]
    
    # Update tabu memory
    if (tabu && generation > 1) {
      tabu_memory[[(generation - 1) %% tabumemsize + 1]] <- elites
    }
    
    # Generate offspring
    offspring <- generate_offspring(
      elites, Candidates, npop - nelite, mutprob, mc.cores, 
      mutintensity, tabu_memory, ntoselect, crossover_method
    )
    
    # Combine elites and offspring
    if (keepbest) {
      population <- c(elites, offspring)
    } else {
      population <- offspring
    }
    
    # Diversity preservation
    if (diversity_preservation) {
      if (diversity_method %in% c("crowding", "both")) {
        # Apply crowding replacement
        offspring_fitness <- evaluate_population(
          offspring, P, Test, criterion, lambda, C, 
          K = NULL, Vg = Vg, Ve = Ve, mc.cores = mc.cores
        )
        
        replacement_result <- crowding_replacement(
          offspring, offspring_fitness, population[1:length(elites)], 
          fitness_values[selected_indices], crowding_factor
        )
        
        population <- c(replacement_result$population, 
                       population[(length(elites)+1):length(population)])
      }
      
      if (diversity_method %in% c("sharing", "both")) {
        # Apply fitness sharing
        current_fitness <- evaluate_population(
          population, P, Test, criterion, lambda, C, 
          K = NULL, Vg = Vg, Ve = Ve, mc.cores = mc.cores
        )
        
        shared_fitness <- fitness_sharing(current_fitness, population, 
                                         sharing_radius, alpha = 1)
        # Use shared fitness for next selection
        fitness_values <- shared_fitness
      }
    }
    
    # Plotting
    if (plotiters && generation %% 10 == 0) {
      if (plottype == 1) {
        plot(1:generation, fitness_history[1:generation, "best"], 
             type = "l", col = "blue", lwd = 2,
             xlab = "Generation", ylab = "Fitness",
             main = paste("GA Progress - Generation", generation))
        lines(1:generation, fitness_history[1:generation, "mean"], col = "red")
        legend("topright", c("Best", "Mean"), col = c("blue", "red"), lty = 1)
      }
    }
    
    # Update progress bar
    if (verbose) {
      update_progress(pb)
    }
  }
  
  # Finish progress bar
  if (verbose) {
    finish_progress(pb)
  }
  
  # Final evaluation if needed
  if (is.null(best_ever_solution)) {
    final_fitness <- evaluate_population(
      population, P, Test, criterion, lambda, C, 
      K = NULL, Vg = Vg, Ve = Ve, mc.cores = mc.cores
    )
    best_idx <- which.min(final_fitness)
    best_ever_solution <- population[[best_idx]]
    best_ever_fitness <- final_fitness[best_idx]
  }
  
  # Return comprehensive results
  list(
    best_solution = best_ever_solution,
    best_fitness = best_ever_fitness,
    fitness_history = fitness_history[1:generation, , drop = FALSE],
    population_stats = generation_stats[1:generation],
    final_population = population,
    convergence_generation = generation,
    parameters = list(
      npop = npop,
      niterations = generation,
      criterion = criterion,
      selection_method = selection_method,
      adaptive_mutation = adaptive_mutation,
      diversity_preservation = diversity_preservation,
      diversity_method = diversity_method
    )
  )
}

#' Genetic algorithm for subset selection without test set (single objective)
#' @param P Prediction matrix (individuals x markers/predictors)
#' @param ntoselect Number of individuals to select
#' @param ... Additional parameters (same as subset_ga but without Test parameter)
#' @return List containing best solution and statistics
subset_ga_single <- function(P, ntoselect, ...) {
  # Call main GA function with Test = NULL
  args <- list(...)
  args$P <- P
  args$Candidates <- rownames(P)
  args$Test <- NULL
  args$ntoselect <- ntoselect
  
  do.call(subset_ga, args)
}

# Helper functions are now implemented inline - no need for sourcing