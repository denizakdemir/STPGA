# Enhanced Genetic Algorithm Demo - STPGA v7.0.0
# This script demonstrates the new genetic algorithm features

library(STPGA)

# Load and prepare data
data("WheatData")
set.seed(123)

# Create manageable subset
subset_indices <- sample(1:nrow(Wheat.M), 150)
M_subset <- Wheat.M[subset_indices, 1:50]  # Use first 50 markers
PC_subset <- prcomp(M_subset, center = TRUE, scale. = TRUE)$x[, 1:5]
rownames(PC_subset) <- rownames(M_subset)

# Define candidates and test set
all_individuals <- rownames(M_subset)
test_set <- sample(all_individuals, 25)
candidates <- setdiff(all_individuals, test_set)

cat("Enhanced GA Demo Setup:\n")
cat("- Total individuals:", length(all_individuals), "\n")
cat("- Candidates:", length(candidates), "\n")
cat("- Test set:", length(test_set), "\n")
cat("- Features:", ncol(PC_subset), "\n\n")

# Demo 1: Multi-criteria convergence detection
cat("=== Demo 1: Enhanced Convergence Detection ===\n")

# Traditional convergence (minimal window)
ga_traditional <- subset_ga(
  P = PC_subset,
  Candidates = candidates,
  Test = test_set,
  ntoselect = 20,
  npop = 30,
  niterations = 60,
  criterion = "pev_mean",
  enable_restart = FALSE,
  convergence_window_multiplier = 1,
  tolconv = 1e-6,
  verbose = FALSE
)

# Enhanced convergence (multi-criteria with larger window)
ga_enhanced <- subset_ga(
  P = PC_subset,
  Candidates = candidates,
  Test = test_set,
  ntoselect = 20,
  npop = 30,
  niterations = 60,
  criterion = "pev_mean",
  enable_restart = FALSE,
  convergence_window_multiplier = 4,
  tolconv = 1e-6,
  verbose = FALSE
)

cat("Traditional convergence generation:", ga_traditional$convergence_generation, "\n")
cat("Enhanced convergence generation:", ga_enhanced$convergence_generation, "\n")
cat("Traditional best fitness:", sprintf("%.6f", ga_traditional$best_fitness), "\n")
cat("Enhanced best fitness:", sprintf("%.6f", ga_enhanced$best_fitness), "\n\n")

# Demo 2: Restart mechanisms
cat("=== Demo 2: Restart Mechanisms ===\n")

# Without restart
ga_no_restart <- subset_ga(
  P = PC_subset,
  Candidates = candidates,
  Test = test_set,
  ntoselect = 20,
  npop = 25,
  niterations = 50,
  criterion = "pev_mean",
  enable_restart = FALSE,
  tolconv = 1e-4,
  minitbefstop = 15,
  verbose = FALSE
)

# With restart
ga_with_restart <- subset_ga(
  P = PC_subset,
  Candidates = candidates,
  Test = test_set,
  ntoselect = 20,
  npop = 25,
  niterations = 50,
  criterion = "pev_mean",
  enable_restart = TRUE,
  restart_threshold = 0.6,
  max_restarts = 2,
  tolconv = 1e-4,
  minitbefstop = 15,
  verbose = FALSE
)

cat("No restart - Best fitness:", sprintf("%.6f", ga_no_restart$best_fitness), "\n")
cat("With restart - Best fitness:", sprintf("%.6f", ga_with_restart$best_fitness), "\n")
cat("Restart count:", ga_with_restart$restart_count, "\n")

if (ga_with_restart$restart_count > 0) {
  cat("Restart history:\n")
  for (i in 1:ga_with_restart$restart_count) {
    restart_info <- ga_with_restart$restart_history[[i]]
    cat(sprintf("  Restart %d: Generation %d, Fitness %.6f\n",
                i, restart_info$generation, restart_info$best_fitness))
  }
}
cat("\n")

# Demo 3: Enhanced selection methods
cat("=== Demo 3: Selection Methods Comparison ===\n")

selection_methods <- c("tournament", "rank", "elite", "hybrid")
selection_results <- list()

for (method in selection_methods) {
  cat(paste("Testing", method, "selection...\n"))
  
  # Set method-specific parameters
  extra_params <- list()
  if (method == "rank") {
    extra_params$selection_pressure <- 1.4
  }
  
  ga_result <- do.call(subset_ga, c(
    list(
      P = PC_subset,
      Candidates = candidates,
      Test = test_set,
      ntoselect = 20,
      npop = 25,
      niterations = 30,
      criterion = "pev_mean",
      selection_method = method,
      enable_restart = FALSE,
      verbose = FALSE
    ),
    extra_params
  ))
  
  selection_results[[method]] <- ga_result
}

cat("\nSelection Method Comparison:\n")
cat("Method      | Best Fitness | Convergence Gen\n")
cat("------------|--------------|---------------\n")
for (method in selection_methods) {
  result <- selection_results[[method]]
  cat(sprintf("%-11s | %12.6f | %14d\n",
              method, result$best_fitness, result$convergence_generation))
}
cat("\n")

# Demo 4: Convergence diagnostics
cat("=== Demo 4: Convergence Diagnostics ===\n")

# Run GA with detailed tracking
ga_diagnostic <- subset_ga(
  P = PC_subset,
  Candidates = candidates,
  Test = test_set,
  ntoselect = 20,
  npop = 30,
  niterations = 40,
  criterion = "pev_mean",
  selection_method = "rank",
  selection_pressure = 1.3,
  enable_restart = TRUE,
  restart_threshold = 0.7,
  convergence_window_multiplier = 4,
  verbose = FALSE
)

# Generate diagnostics
diags <- convergence_diagnostics(ga_diagnostic, plot = FALSE)

cat("Convergence Status:\n")
cat("- Converged:", diags$convergence_status$converged, "\n")
cat("- Reason:", diags$convergence_status$reason, "\n")
cat("- Window size:", diags$window_size, "\n")

cat("\nMetrics:\n")
cat("- Fitness stability:", sprintf("%.8f", diags$metrics$fitness_stability), "\n")
cat("- Diversity stability:", sprintf("%.6f", diags$metrics$diversity_stability), "\n")
cat("- Mean improvement rate:", sprintf("%.8f", diags$metrics$mean_improvement_rate), "\n")

if (length(diags$recommendations) > 0) {
  cat("\nRecommendations:\n")
  for (rec in diags$recommendations) {
    cat(paste("-", rec, "\n"))
  }
}
cat("\n")

# Demo 5: Parameter validation
cat("=== Demo 5: Parameter Validation ===\n")

# Test parameter validation
validation_result <- validate_selection_parameters(
  selection_method = "rank",
  tournament_size = 100,  # Too large - will be corrected
  selection_pressure = 2.5,  # Out of range - will be corrected  
  population_size = 25
)

cat("Parameter Validation Results:\n")
cat("- Selection method:", validation_result$selection_method, "\n")
cat("- Tournament size (corrected):", validation_result$tournament_size, "\n")
cat("- Selection pressure (corrected):", validation_result$selection_pressure, "\n")

if (length(validation_result$warnings) > 0) {
  cat("\nWarnings/Corrections:\n")
  for (warning in validation_result$warnings) {
    cat(paste("-", warning, "\n"))
  }
}

cat("\n=== Demo Complete ===\n")
cat("This demo showcased the key enhancements in STPGA v7.0.0:\n")
cat("1. Multi-criteria convergence detection\n")
cat("2. Restart mechanisms for escaping local optima\n")
cat("3. Enhanced selection methods with rank-based selection\n")
cat("4. Comprehensive convergence diagnostics\n")
cat("5. Automatic parameter validation and correction\n\n")
cat("For more detailed examples, see the 'Advanced Genetic Algorithm Features' vignette.\n")