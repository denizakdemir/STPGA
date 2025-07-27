# Demonstration of Genetic Algorithm Features in STPGA
# This script shows the main algorithm functionality

library(STPGA)

# Set seed for reproducibility
set.seed(123)

# Generate synthetic data for demonstration
n <- 100  # Number of individuals
p <- 20   # Number of markers

# Create synthetic genotype matrix
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
rownames(X) <- paste0("Ind", 1:n)
colnames(X) <- paste0("Marker", 1:p)

# Define training candidates and test set
test_indices <- sample(1:n, 30)
train_candidates <- rownames(X)[-test_indices]
test_set <- rownames(X)[test_indices]

cat("STPGA Genetic Algorithm Demo\n")
cat("============================\n\n")

# Basic usage with default adaptive parameters
cat("Running genetic algorithm for subset selection...\n")
result <- GenAlgForSubsetSelection(
  P = X,
  Candidates = train_candidates,
  Test = test_set,
  ntoselect = 15,
  npop = 30,
  nelite = 5,
  niterations = 50,
  errorstat = "pev_mean_normalized",
  verbose = TRUE
)

# Extract best solution
best_solution <- result[[1]]  # First element is best solution
best_fitness <- min(result[["Best criterion values over iterarions"]])

cat(sprintf("\nBest solution found: %d individuals selected\n", length(best_solution)))
cat(sprintf("Best fitness achieved: %.6f\n", best_fitness))

# Show convergence
convergence <- result[["Best criterion values over iterarions"]]
cat(sprintf("Convergence: %.6f -> %.6f (%.1f%% improvement)\n", 
            convergence[1], tail(convergence, 1),
            100 * (convergence[1] - tail(convergence, 1)) / convergence[1]))

# Display algorithm parameters used
if (!is.null(result$GAParameters)) {
  cat("\nAlgorithm configuration:\n")
  cat(sprintf("- Selection method: %s\n", result$GAParameters$selectionMethod))
  cat(sprintf("- Crossover method: %s\n", result$GAParameters$crossoverMethod))
  cat(sprintf("- Adaptive mutation: %s\n", result$GAParameters$adaptiveMutation))
  cat(sprintf("- Diversity preservation: %s\n", result$GAParameters$diversityPreservation))
  if (!is.na(result$GAParameters$finalDiversity)) {
    cat(sprintf("- Final population diversity: %.3f\n", result$GAParameters$finalDiversity))
  }
}

# Plot convergence if running interactively
if (interactive()) {
  par(mfrow = c(2, 2))
  
  # Convergence curve
  plot(convergence, type = "l", col = "blue", lwd = 2,
       xlab = "Generation", ylab = "Best Fitness", 
       main = "Algorithm Convergence")
  points(length(convergence), tail(convergence, 1), col = "red", pch = 19, cex = 1.5)
  
  # Show diversity over time if available
  if (!is.null(result$ConvergenceHistory) && !all(is.na(result$ConvergenceHistory$Diversity))) {
    plot(result$ConvergenceHistory$Diversity, type = "l", col = "green", lwd = 2,
         xlab = "Generation", ylab = "Population Diversity", 
         main = "Population Diversity")
    abline(h = 0.3, col = "red", lty = 2)
    legend("topright", legend = c("Actual", "Target"), col = c("green", "red"), lty = c(1, 2))
  }
  
  # Show mutation rate adaptation if available
  if (!is.null(result$ConvergenceHistory) && !all(is.na(result$ConvergenceHistory$MutationRate))) {
    plot(result$ConvergenceHistory$MutationRate, type = "l", col = "purple", lwd = 2,
         xlab = "Generation", ylab = "Mutation Rate", 
         main = "Adaptive Mutation Rate")
  }
  
  # Show fitness statistics
  if (!is.null(result$ConvergenceHistory)) {
    plot(result$ConvergenceHistory$Best, type = "l", col = "blue", lwd = 2,
         ylim = range(c(result$ConvergenceHistory$Best, result$ConvergenceHistory$Mean), na.rm = TRUE),
         xlab = "Generation", ylab = "Fitness", main = "Population Statistics")
    if (!all(is.na(result$ConvergenceHistory$Mean))) {
      lines(result$ConvergenceHistory$Mean, col = "orange", lwd = 2)
    }
    legend("topright", legend = c("Best", "Mean"), col = c("blue", "orange"), lty = 1, lwd = 2)
  }
  
  par(mfrow = c(1, 1))
}

# Demonstrate modern criteria functions
cat("\nTesting modern criteria functions on best solution:\n")
cat("=================================================\n")

# Evaluate best solution with different criteria
criteria_results <- list(
  a_optimality = a_optimality(best_solution, test_set, X),
  d_optimality = d_optimality(best_solution, test_set, X),
  e_optimality = e_optimality(best_solution, test_set, X),
  pev_mean = pev_mean(best_solution, test_set, X),
  pev_mean_norm = pev_mean(best_solution, test_set, X, normalized = TRUE),
  pev_max = pev_max(best_solution, test_set, X),
  cd_mean = cd_mean(best_solution, test_set, X)
)

# Display results
for (i in seq_along(criteria_results)) {
  criterion_name <- names(criteria_results)[i]
  criterion_value <- criteria_results[[i]]
  cat(sprintf("%-15s: %12.8f\n", criterion_name, criterion_value))
}

# Compare with random selection
cat("\nComparison with random selection:\n")
random_selection <- sample(train_candidates, length(best_solution))
random_pev <- pev_mean(random_selection, test_set, X, normalized = TRUE)
optimized_pev <- criteria_results$pev_mean_norm

cat(sprintf("Random selection PEV:    %.8f\n", random_pev))
cat(sprintf("Optimized selection PEV: %.8f\n", optimized_pev))
cat(sprintf("Improvement:             %.2f%%\n", 
            100 * (random_pev - optimized_pev) / random_pev))

cat("\nDemo completed successfully!\n")

# Show different parameter configurations
cat("\nAlternative configurations:\n")
cat("============================\n")
cat("# Different optimization criteria:\n")
cat("result <- GenAlgForSubsetSelection(..., errorstat = 'a_optimality')  # Minimize average variance\n")
cat("result <- GenAlgForSubsetSelection(..., errorstat = 'd_optimality')  # Minimize generalized variance\n")
cat("result <- GenAlgForSubsetSelection(..., errorstat = 'pev_mean')      # Minimize prediction error\n")
cat("result <- GenAlgForSubsetSelection(..., errorstat = 'cd_mean')       # Minimize Cook's distance\n\n")

cat("# For faster but less thorough search:\n")
cat("result <- GenAlgForSubsetSelection(..., selectionMethod = 'elite', crossoverMethod = 'traditional')\n\n")
cat("# For maximum diversity preservation:\n") 
cat("result <- GenAlgForSubsetSelection(..., diversityMethod = 'both', sharingRadius = 0.2)\n\n")
cat("# For detailed monitoring:\n")
cat("result <- GenAlgForSubsetSelection(..., verbose = TRUE, plotiters = TRUE)\n\n")

cat("# Multi-objective optimization (Pareto optimal solutions):\n")
cat("mo_result <- GenAlgForSubsetSelectionMO(\n")
cat("  Pcs = X, candidates = train_candidates, Test = test_set,\n")
cat("  selectionstats = c('a_optimality', 'pev_mean'),\n")
cat("  selectionstatstypes = c(TRUE, TRUE)  # Both minimize\n")
cat(")\n")