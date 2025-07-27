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
  
  # Initialize internal helper functions
  source_ga_operators()
  source_evaluation_functions()
  
  # Ridge regression function for internal use
  ridge_regression <- function(y, x, lambda = lambda) {
    n <- nrow(x)
    p <- ncol(x)
    mindim <- min(p, n)
    rownames(x) <- NULL
    
    if (p >= n) {
      # More predictors than observations - use eigendecomposition
      eigen_decomp <- eigen(crossprod(x) / (n - 1) + lambda * diag(p), symmetric = TRUE)
      eigenvectors <- eigen_decomp$vectors
      eigenvalues <- eigen_decomp$values
      
      # Compute coefficients
      beta <- eigenvectors %*% diag(1 / eigenvalues) %*% t(eigenvectors) %*% t(x) %*% y
      return(beta)
    } else {
      # Standard case - more observations than predictors
      XtX <- crossprod(x)
      Xty <- crossprod(x, y)
      
      # Use Cholesky decomposition for numerical stability
      tryCatch({
        L <- chol(XtX + lambda * diag(p))
        beta <- backsolve(L, forwardsolve(t(L), Xty))
        return(beta)
      }, error = function(e) {
        # Fallback to SVD if Cholesky fails
        svd_result <- svd(x)
        d <- svd_result$d
        d_reg <- d^2 + lambda
        beta <- svd_result$v %*% diag(d / d_reg) %*% t(svd_result$u) %*% y
        return(beta)
      })
    }
  }
  
  # Main GA execution logic will be inserted here
  # This is a placeholder - the actual implementation follows the original logic
  # but with cleaner organization and the new function names
  
  # Return results in standardized format
  list(
    best_solution = NULL,  # Best training set found
    best_fitness = NULL,   # Best fitness value
    fitness_history = NULL, # Fitness evolution over generations
    population_stats = NULL, # Population diversity and convergence stats
    parameters = list(     # Algorithm parameters used
      npop = npop,
      niterations = niterations,
      criterion = criterion,
      selection_method = selection_method,
      adaptive_mutation = adaptive_mutation
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

# Helper function to source operator functions (to be implemented)
source_ga_operators <- function() {
  # This will load the operator functions when they're created
  invisible(NULL)
}

# Helper function to source evaluation functions (to be implemented)  
source_evaluation_functions <- function() {
  # This will load the evaluation functions when they're created
  invisible(NULL)
}