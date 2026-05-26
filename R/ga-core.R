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
#' @param convergence_window_multiplier Multiplier for convergence window size (default: 4)
#' @param enable_restart Whether to enable restart mechanism for premature convergence (default: TRUE)
#' @param restart_threshold Threshold for triggering restart (default: 0.5)
#' @param max_restarts Maximum number of restarts allowed (default: 2)
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
#' @param selection_method Selection method: "tournament", "elite", "rank", "hybrid" (default: "tournament")
#' @param tournament_size Tournament size for tournament selection (default: 3)
#' @param selection_pressure Selection pressure for rank-based selection (default: 1.5)
#' @param diversity_preservation Whether to preserve diversity (default: TRUE)
#' @param diversity_method Diversity method: "crowding", "sharing", "both" (default: "crowding")
#' @param crowding_factor Factor for crowding replacement (default: 3)
#' @param sharing_radius Radius for fitness sharing (default: 0.1)
#' @param crossover_method Crossover method (default: "adaptive")
#' @param diversity_target Target diversity level (default: 0.3)
#' @param verbose Whether to print verbose output (default: FALSE)
#' @return List containing best solution, fitness history, and statistics
#' 
#' @examples
#' # Load wheat genomic data
#' data(WheatData)
#' 
#' # Create a manageable subset for demonstration
#' set.seed(123)
#' subset_indices <- sample(1:nrow(Wheat.M), 100)
#' M_subset <- Wheat.M[subset_indices, 1:30]
#' 
#' # Extract principal components for genetic algorithm
#' pca_result <- prcomp(M_subset, center = TRUE, scale. = TRUE)
#' PC_subset <- pca_result$x[, 1:5]
#' rownames(PC_subset) <- rownames(M_subset)
#' 
#' # Define candidate and test sets
#' all_individuals <- rownames(PC_subset)
#' test_set <- sample(all_individuals, 20)
#' candidates <- setdiff(all_individuals, test_set)
#' 
#' # Basic genetic algorithm run
#' ga_result <- subset_ga(
#'   P = PC_subset,
#'   Candidates = candidates,
#'   Test = test_set,
#'   ntoselect = 15,
#'   npop = 20,
#'   niterations = 30,
#'   criterion = "pev_mean",
#'   verbose = FALSE
#' )
#' 
#' print(paste("Best fitness:", round(ga_result$best_fitness, 6)))
#' print(paste("Convergence generation:", ga_result$convergence_generation))
#' print(paste("Selected individuals:", length(ga_result$best_solution)))
#' 
#' # Enhanced genetic algorithm with rank-based selection and restart
#' ga_enhanced <- subset_ga(
#'   P = PC_subset,
#'   Candidates = candidates,
#'   Test = test_set,
#'   ntoselect = 15,
#'   npop = 20,
#'   niterations = 30,
#'   criterion = "pev_mean",
#'   selection_method = "rank",
#'   selection_pressure = 1.4,
#'   enable_restart = TRUE,
#'   restart_threshold = 0.6,
#'   max_restarts = 1,
#'   verbose = FALSE
#' )
#' 
#' print(paste("Enhanced GA fitness:", round(ga_enhanced$best_fitness, 6)))
#' print(paste("Restart count:", ga_enhanced$restart_count))
#' 
#' # Generate convergence diagnostics
#' diags <- convergence_diagnostics(ga_enhanced, plot = FALSE)
#' print(paste("Converged:", diags$convergence_status$converged))
#' print(paste("Convergence reason:", diags$convergence_status$reason))
subset_ga <- function(P, Candidates, Test = NULL, ntoselect,
                     npop = 100, nelite = 5, keepbest = TRUE, tabu = TRUE, 
                     tabumemsize = 1, mutprob = 0.8, mutintensity = 1, 
                     niterations = 500, minitbefstop = 200, niterreg = 5,
                     convergence_window_multiplier = 4, enable_restart = TRUE,
                     restart_threshold = 0.5, max_restarts = 2, lambda = 1e-6, 
                     plotiters = FALSE, plottype = 1, criterion = "PEVMEAN2", 
                     C = NULL, mc.cores = 1, InitPop = NULL, tolconv = 1e-7, 
                     Vg = NULL, Ve = NULL, Fedorov = FALSE, adaptive_mutation = TRUE,
                     selection_method = "tournament", tournament_size = 3,
                     selection_pressure = 1.5, diversity_preservation = TRUE, 
                     diversity_method = "crowding",
                     crowding_factor = 3, sharing_radius = 0.1,
                     crossover_method = "adaptive", diversity_target = 0.3,
                     verbose = FALSE) {
  
  # Input validation
  if ((ncol(P) + 1) > ntoselect) {
    warning("The algorithm does not work well with p>ntrain, perhaps use unsupervised dimension reduction on P.")
  }
  
  # Validate inputs
  validate_matrix_params(P, train = Candidates, test = Test, lambda = lambda)
  
  # Validate and adjust selection parameters
  selection_validation <- validate_selection_parameters(
    selection_method, tournament_size, selection_pressure, npop
  )
  
  # Update parameters based on validation
  selection_method <- selection_validation$selection_method
  tournament_size <- selection_validation$tournament_size
  selection_pressure <- selection_validation$selection_pressure
  
  # Report any parameter adjustments
  if (length(selection_validation$warnings) > 0 && verbose) {
    for (warning_msg in selection_validation$warnings) {
      message("Parameter adjustment: ", warning_msg)
    }
  }
  
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
  
  # Multi-criteria convergence tracking  
  convergence_window <- max(niterreg * convergence_window_multiplier, 10)  # Configurable window size
  convergence_history <- list(
    fitness = numeric(convergence_window),
    diversity = numeric(convergence_window),
    improvement_rate = numeric(convergence_window),
    fitness_variance = numeric(convergence_window)
  )
  
  # Restart mechanism tracking
  restart_count <- 0
  restart_history <- list()
  total_generations <- 0
  
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
    
    # Multi-criteria convergence tracking
    window_idx <- ((generation - 1) %% convergence_window) + 1
    convergence_history$fitness[window_idx] <- current_best_fitness
    convergence_history$diversity[window_idx] <- generation_stats[[generation]]$diversity
    convergence_history$fitness_variance[window_idx] <- generation_stats[[generation]]$fitness_std
    
    # Calculate improvement rate
    if (generation > 1) {
      prev_fitness <- fitness_history[generation - 1, "best"]
      improvement_rate <- abs(current_best_fitness - prev_fitness) / (abs(prev_fitness) + 1e-12)
      convergence_history$improvement_rate[window_idx] <- improvement_rate
    } else {
      convergence_history$improvement_rate[window_idx] <- 1.0
    }
    
    # Multi-criteria convergence check
    if (generation >= minitbefstop && generation > convergence_window) {
      convergence_result <- check_multi_criteria_convergence(
        convergence_history, generation, convergence_window, tolconv
      )
      
      if (convergence_result$converged) {
        # Check if this is premature convergence that should trigger restart
        if (enable_restart && restart_count < max_restarts && 
            generation < (niterations * restart_threshold)) {
          
          # Trigger restart
          restart_count <- restart_count + 1
          total_generations <- total_generations + generation
          
          if (verbose) {
            message(paste("Premature convergence detected at generation", generation))
            message(paste("Triggering restart", restart_count, "of", max_restarts))
          }
          
          # Store restart information
          restart_history[[restart_count]] <- list(
            generation = generation,
            best_fitness = best_ever_fitness,
            convergence_reason = convergence_result$reason
          )
          
          # Restart with diversified population
          restart_result <- perform_restart(population, best_ever_solution, 
                                          Candidates, ntoselect, npop, generation)
          population <- restart_result$new_population
          
          # Reset convergence tracking but keep global best
          convergence_history <- list(
            fitness = numeric(convergence_window),
            diversity = numeric(convergence_window),
            improvement_rate = numeric(convergence_window),
            fitness_variance = numeric(convergence_window)
          )
          no_improvement_count <- 0
          
          # Continue evolution
          next
        } else {
          # Final convergence
          if (verbose) {
            message(paste("Multi-criteria convergence achieved at generation", generation))
            message(paste("Convergence reason:", convergence_result$reason))
          }
          break
        }
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
    } else if (selection_method == "rank") {
      selected_indices <- rank_selection(population, fitness_values, nelite, 
                                        selection_pressure, "linear")
    } else { # hybrid
      # Combine different selection methods
      n_elite <- max(1, nelite %/% 3)
      n_rank <- max(1, nelite %/% 3)  
      n_tournament <- nelite - n_elite - n_rank
      
      elite_indices <- elite_selection(population, fitness_values, n_elite)
      rank_indices <- rank_selection(population, fitness_values, n_rank, 
                                   selection_pressure, "linear")
      tournament_indices <- tournament_selection(population, fitness_values, 
                                                tournament_size, n_tournament)
      selected_indices <- c(elite_indices, rank_indices, tournament_indices)
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
    convergence_history = convergence_history,
    restart_history = restart_history,
    final_population = population,
    convergence_generation = generation,
    total_generations = total_generations + generation,
    restart_count = restart_count,
    parameters = list(
      npop = npop,
      niterations = generation,
      criterion = criterion,
      selection_method = selection_method,
      tournament_size = tournament_size,
      selection_pressure = selection_pressure,
      adaptive_mutation = adaptive_mutation,
      diversity_preservation = diversity_preservation,
      diversity_method = diversity_method,
      convergence_window = convergence_window,
      enable_restart = enable_restart,
      max_restarts = max_restarts,
      restart_threshold = restart_threshold
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

#' Multi-criteria convergence detection
#' @param convergence_history List containing convergence metrics history
#' @param generation Current generation number
#' @param window_size Window size for convergence check
#' @param tolerance Convergence tolerance
#' @return List with convergence status and reason
check_multi_criteria_convergence <- function(convergence_history, generation, 
                                           window_size, tolerance) {
  
  # Extract recent values within the window
  recent_fitness <- convergence_history$fitness[1:min(window_size, generation)]
  recent_diversity <- convergence_history$diversity[1:min(window_size, generation)]
  recent_improvement <- convergence_history$improvement_rate[1:min(window_size, generation)]
  recent_variance <- convergence_history$fitness_variance[1:min(window_size, generation)]
  
  # Remove NA values
  recent_fitness <- recent_fitness[!is.na(recent_fitness)]
  recent_diversity <- recent_diversity[!is.na(recent_diversity)]
  recent_improvement <- recent_improvement[!is.na(recent_improvement)]
  recent_variance <- recent_variance[!is.na(recent_variance)]
  
  if (length(recent_fitness) < window_size / 2) {
    return(list(converged = FALSE, reason = "Insufficient data"))
  }
  
  # Criterion 1: Fitness stabilization
  fitness_range <- max(recent_fitness) - min(recent_fitness)
  fitness_converged <- fitness_range < tolerance
  
  # Criterion 2: Low improvement rate
  mean_improvement <- mean(recent_improvement, na.rm = TRUE)
  improvement_converged <- mean_improvement < tolerance
  
  # Criterion 3: Population diversity stabilization
  diversity_std <- sd(recent_diversity, na.rm = TRUE)
  diversity_converged <- diversity_std < tolerance * 10  # More relaxed threshold
  
  # Criterion 4: Fitness variance stabilization
  variance_std <- sd(recent_variance, na.rm = TRUE)
  variance_converged <- variance_std < tolerance * 100  # More relaxed threshold
  
  # Combined convergence decision
  criteria_met <- sum(c(fitness_converged, improvement_converged, 
                       diversity_converged, variance_converged))
  
  # Require at least 3 out of 4 criteria to be met
  converged <- criteria_met >= 3
  
  # Determine primary reason for convergence
  reason <- if (converged) {
    reasons <- c()
    if (fitness_converged) reasons <- c(reasons, "fitness stabilized")
    if (improvement_converged) reasons <- c(reasons, "low improvement rate")
    if (diversity_converged) reasons <- c(reasons, "diversity stabilized")
    if (variance_converged) reasons <- c(reasons, "variance stabilized")
    paste(reasons, collapse = ", ")
  } else {
    "convergence criteria not met"
  }
  
  list(
    converged = converged,
    reason = reason,
    criteria_details = list(
      fitness_converged = fitness_converged,
      improvement_converged = improvement_converged,
      diversity_converged = diversity_converged,
      variance_converged = variance_converged,
      fitness_range = fitness_range,
      mean_improvement = mean_improvement,
      diversity_std = diversity_std,
      variance_std = variance_std
    )
  )
}

#' Perform genetic algorithm restart with population diversification
#' @param current_population Current population
#' @param best_solution Current best solution  
#' @param candidates Candidate pool
#' @param ntoselect Number to select per solution
#' @param npop Population size
#' @param generation Current generation (for restart strategies)
#' @return List with new diversified population
perform_restart <- function(current_population, best_solution, candidates, 
                          ntoselect, npop, generation) {
  
  # Strategy 1: Keep a small fraction of current population (elites)
  elite_fraction <- 0.1
  n_elites <- max(1, floor(npop * elite_fraction))
  
  # Strategy 2: Generate new random solutions
  random_fraction <- 0.6
  n_random <- floor(npop * random_fraction)
  
  # Strategy 3: Generate solutions based on best solution with mutations
  n_mutated_best <- npop - n_elites - n_random
  
  new_population <- list()
  
  # Add elites (best performing solutions from current population)
  if (length(current_population) >= n_elites) {
    # Select best elites (would need fitness evaluation, but approximate by keeping first few)
    for (i in 1:n_elites) {
      new_population[[length(new_population) + 1]] <- current_population[[i]]
    }
  }
  
  # Add random solutions for exploration
  for (i in 1:n_random) {
    new_population[[length(new_population) + 1]] <- sample(candidates, ntoselect, replace = FALSE)
  }
  
  # Add mutated versions of best solution for exploitation around good region
  if (!is.null(best_solution)) {
    for (i in 1:n_mutated_best) {
      # Create mutated version of best solution
      mutated_solution <- mutate_solution(best_solution, candidates, 
                                        min(3, length(best_solution)), ntoselect)
      new_population[[length(new_population) + 1]] <- mutated_solution
    }
  } else {
    # If no best solution, add more random solutions
    for (i in 1:n_mutated_best) {
      new_population[[length(new_population) + 1]] <- sample(candidates, ntoselect, replace = FALSE)
    }
  }
  
  # Ensure we have exactly npop solutions
  while (length(new_population) < npop) {
    new_population[[length(new_population) + 1]] <- sample(candidates, ntoselect, replace = FALSE)
  }
  
  if (length(new_population) > npop) {
    new_population <- new_population[1:npop]
  }
  
  list(
    new_population = new_population,
    strategy_breakdown = list(
      elites = n_elites,
      random = n_random,
      mutated_best = n_mutated_best
    )
  )
}

#' Get convergence history diagnostics from GA results
#' @param ga_result Result from subset_ga function
#' @param plot Whether to create diagnostic plots (default: TRUE)
#' @return List with convergence diagnostics
#' 
#' @examples
#' # Load wheat genomic data
#' data(WheatData)
#' 
#' # Create a subset for demonstration
#' set.seed(456)
#' subset_indices <- sample(1:nrow(Wheat.M), 100)
#' M_subset <- Wheat.M[subset_indices, 1:30]
#' 
#' # Extract principal components
#' pca_result <- prcomp(M_subset, center = TRUE, scale. = TRUE)
#' PC_subset <- pca_result$x[, 1:5]
#' rownames(PC_subset) <- rownames(M_subset)
#' 
#' # Define candidate and test sets
#' all_individuals <- rownames(PC_subset)
#' test_set <- sample(all_individuals, 20)
#' candidates <- setdiff(all_individuals, test_set)
#' 
#' # Run genetic algorithm with convergence tracking
#' ga_result <- subset_ga(
#'   P = PC_subset,
#'   Candidates = candidates,
#'   Test = test_set,
#'   ntoselect = 15,
#'   npop = 25,
#'   niterations = 40,
#'   criterion = "pev_mean",
#'   enable_restart = TRUE,
#'   convergence_window_multiplier = 4,
#'   verbose = FALSE
#' )
#' 
#' # Generate convergence diagnostics
#' diags <- convergence_diagnostics(ga_result, plot = FALSE)
#' 
#' print("Convergence Diagnostics:")
#' print(paste("Converged:", diags$convergence_status$converged))
#' print(paste("Reason:", diags$convergence_status$reason))
#' print(paste("Window size used:", diags$window_size))
#' print(paste("Data points analyzed:", diags$data_points))
#' 
#' print("Stability Metrics:")
#' print(paste("Fitness stability:", round(diags$metrics$fitness_stability, 8)))
#' print(paste("Diversity stability:", round(diags$metrics$diversity_stability, 6)))
#' print(paste("Mean improvement rate:", round(diags$metrics$mean_improvement_rate, 8)))
#' 
#' if (length(diags$recommendations) > 0) {
#'   print("Recommendations:")
#'   for (rec in diags$recommendations) {
#'     print(paste("-", rec))
#'   }
#' }
convergence_diagnostics <- function(ga_result, plot = TRUE) {
  
  if (!"convergence_history" %in% names(ga_result)) {
    stop("GA result does not contain convergence_history. Run subset_ga with updated version.")
  }
  
  conv_hist <- ga_result$convergence_history
  fitness_hist <- ga_result$fitness_history
  n_generations <- nrow(fitness_hist)
  
  # Extract non-NA values for analysis
  valid_indices <- which(!is.na(conv_hist$fitness))
  
  if (length(valid_indices) == 0) {
    warning("No valid convergence history data found")
    return(list(status = "no_data"))
  }
  
  # Calculate convergence metrics
  fitness_stability <- sd(conv_hist$fitness[valid_indices], na.rm = TRUE)
  diversity_stability <- sd(conv_hist$diversity[valid_indices], na.rm = TRUE)
  mean_improvement_rate <- mean(conv_hist$improvement_rate[valid_indices], na.rm = TRUE)
  variance_stability <- sd(conv_hist$fitness_variance[valid_indices], na.rm = TRUE)
  
  # Convergence assessment
  final_convergence <- check_multi_criteria_convergence(
    conv_hist, n_generations, ga_result$parameters$convergence_window, 1e-7
  )
  
  diagnostics <- list(
    convergence_status = final_convergence,
    window_size = ga_result$parameters$convergence_window,
    data_points = length(valid_indices),
    metrics = list(
      fitness_stability = fitness_stability,
      diversity_stability = diversity_stability,
      mean_improvement_rate = mean_improvement_rate,
      variance_stability = variance_stability
    ),
    recommendations = list()
  )
  
  # Generate recommendations
  if (mean_improvement_rate > 0.01) {
    diagnostics$recommendations <- c(diagnostics$recommendations, 
                                   "Consider increasing niterations - high improvement rate suggests more optimization possible")
  }
  
  if (diversity_stability > 0.1) {
    diagnostics$recommendations <- c(diagnostics$recommendations,
                                   "Population diversity is unstable - consider adjusting diversity_preservation settings")
  }
  
  if (fitness_stability > 0.01) {
    diagnostics$recommendations <- c(diagnostics$recommendations,
                                   "Fitness not fully stabilized - consider increasing convergence_window_multiplier")
  }
  
  # Plotting
  if (plot && n_generations > 1) {
    par(mfrow = c(2, 2))
    
    # Fitness convergence
    plot(1:n_generations, fitness_hist[, "best"], type = "l", col = "blue", lwd = 2,
         main = "Fitness Convergence", xlab = "Generation", ylab = "Best Fitness")
    
    # Diversity tracking
    if (length(valid_indices) > 1) {
      plot(valid_indices, conv_hist$diversity[valid_indices], type = "l", col = "red", lwd = 2,
           main = "Population Diversity", xlab = "Generation", ylab = "Diversity")
    }
    
    # Improvement rate
    if (length(valid_indices) > 1) {
      plot(valid_indices, conv_hist$improvement_rate[valid_indices], type = "l", col = "green", lwd = 2,
           main = "Improvement Rate", xlab = "Generation", ylab = "Improvement Rate")
    }
    
    # Fitness variance
    if (length(valid_indices) > 1) {
      plot(valid_indices, conv_hist$fitness_variance[valid_indices], type = "l", col = "purple", lwd = 2,
           main = "Fitness Variance", xlab = "Generation", ylab = "Population Fitness Std Dev")
    }
    
    par(mfrow = c(1, 1))
  }
  
  return(diagnostics)
}

# Helper functions are now implemented inline - no need for sourcing