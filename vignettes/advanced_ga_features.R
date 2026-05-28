## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6
)

## ----load-package-------------------------------------------------------------
library(STPGA)

# Load wheat genomic data
data("WheatData")

# Create a manageable subset for demonstration
set.seed(123)
n_total <- nrow(Wheat.M)
subset_indices <- sample(1:n_total, 200)  # Use 200 individuals for faster computation

# Extract subset
M_subset <- Wheat.M[subset_indices, 1:100]  # Use first 100 markers
K_subset <- Wheat.K[subset_indices, subset_indices]

# Create principal components for non-mixed model criteria
pca_result <- prcomp(M_subset, center = TRUE, scale. = TRUE)
PC_subset <- pca_result$x[, 1:5]
rownames(PC_subset) <- rownames(M_subset)

# Define training candidates and test set
all_individuals <- rownames(M_subset)
test_set <- sample(all_individuals, 30)
candidates <- setdiff(all_individuals, test_set)

cat("Dataset Summary:\n")
cat("- Total individuals:", length(all_individuals), "\n")
cat("- Candidates:", length(candidates), "\n")
cat("- Test set:", length(test_set), "\n")
cat("- Markers:", ncol(M_subset), "\n")
cat("- Principal components:", ncol(PC_subset), "\n")

## ----convergence-comparison---------------------------------------------------
# Parameters for comparison
n_select <- 25
npop <- 40
niterations <- 100

cat("Comparing Convergence Detection Methods\n")
cat("=====================================\n")

# Simulate traditional convergence (disable restart and enhanced features for comparison)
set.seed(456)
ga_traditional <- subset_ga(
  P = PC_subset,
  Candidates = candidates,
  Test = test_set,
  ntoselect = n_select,
  npop = npop,
  niterations = niterations,
  criterion = "pev_mean",
  enable_restart = FALSE,
  tolconv = 1e-6,
  niterreg = 5,
  convergence_window_multiplier = 1,  # Minimal window
  verbose = FALSE
)

# Enhanced convergence detection
set.seed(456)  # Same seed for fair comparison
ga_enhanced <- subset_ga(
  P = PC_subset,
  Candidates = candidates,
  Test = test_set,
  ntoselect = n_select,
  npop = npop,
  niterations = niterations,
  criterion = "pev_mean",
  enable_restart = FALSE,  # Disable restart to see pure convergence detection
  tolconv = 1e-6,
  convergence_window_multiplier = 4,
  verbose = FALSE
)

cat("Traditional Approach:\n")
cat("- Convergence generation:", ga_traditional$convergence_generation, "\n")
cat("- Best fitness:", sprintf("%.6f", ga_traditional$best_fitness), "\n")
cat("- Window size:", ga_traditional$parameters$convergence_window, "\n")

cat("\nEnhanced Approach:\n")
cat("- Convergence generation:", ga_enhanced$convergence_generation, "\n") 
cat("- Best fitness:", sprintf("%.6f", ga_enhanced$best_fitness), "\n")
cat("- Window size:", ga_enhanced$parameters$convergence_window, "\n")

## ----plot-convergence-comparison, fig.width=10, fig.height=5------------------
# Plot convergence behavior
par(mfrow = c(1, 2))

# Plot fitness convergence
plot(1:nrow(ga_traditional$fitness_history), ga_traditional$fitness_history[, "best"],
     type = "l", col = "blue", lwd = 2,
     xlab = "Generation", ylab = "Best Fitness",
     main = "Convergence Detection Comparison")
lines(1:nrow(ga_enhanced$fitness_history), ga_enhanced$fitness_history[, "best"],
      col = "red", lwd = 2)

# Mark convergence points
abline(v = ga_traditional$convergence_generation, col = "blue", lty = 2)
abline(v = ga_enhanced$convergence_generation, col = "red", lty = 2)

legend("topright", 
       legend = c("Traditional", "Enhanced", "Conv. Traditional", "Conv. Enhanced"),
       col = c("blue", "red", "blue", "red"),
       lty = c(1, 1, 2, 2), lwd = c(2, 2, 1, 1))

# Plot population diversity comparison
trad_diversity <- sapply(ga_traditional$population_stats, function(x) x$diversity)
enh_diversity <- sapply(ga_enhanced$population_stats, function(x) x$diversity)

plot(1:length(trad_diversity), trad_diversity,
     type = "l", col = "blue", lwd = 2,
     xlab = "Generation", ylab = "Population Diversity",
     main = "Population Diversity Evolution")
lines(1:length(enh_diversity), enh_diversity, col = "red", lwd = 2)

abline(v = ga_traditional$convergence_generation, col = "blue", lty = 2)
abline(v = ga_enhanced$convergence_generation, col = "red", lty = 2)

legend("topright", 
       legend = c("Traditional", "Enhanced"),
       col = c("blue", "red"), lwd = 2)

par(mfrow = c(1, 1))

## ----restart-demonstration----------------------------------------------------
cat("Restart Mechanism Demonstration\n")
cat("==============================\n")

# GA without restart (may get stuck in local optima)
set.seed(789)
ga_no_restart <- subset_ga(
  P = PC_subset,
  Candidates = candidates,
  Test = test_set,
  ntoselect = n_select,
  npop = 30,
  niterations = 80,
  criterion = "pev_mean",
  enable_restart = FALSE,
  tolconv = 1e-4,  # More relaxed for quicker convergence
  minitbefstop = 20,
  verbose = FALSE
)

# GA with restart enabled
set.seed(789)  # Same seed for comparison
ga_with_restart <- subset_ga(
  P = PC_subset,
  Candidates = candidates,
  Test = test_set,
  ntoselect = n_select,
  npop = 30,
  niterations = 80,
  criterion = "pev_mean",
  enable_restart = TRUE,
  restart_threshold = 0.5,  # Restart if converged before 50% of iterations
  max_restarts = 2,
  tolconv = 1e-4,
  minitbefstop = 20,
  verbose = FALSE
)

cat("Without Restart:\n")
cat("- Best fitness:", sprintf("%.6f", ga_no_restart$best_fitness), "\n")
cat("- Convergence generation:", ga_no_restart$convergence_generation, "\n")
cat("- Total generations:", ga_no_restart$total_generations, "\n")

cat("\nWith Restart:\n")
cat("- Best fitness:", sprintf("%.6f", ga_with_restart$best_fitness), "\n")
cat("- Convergence generation:", ga_with_restart$convergence_generation, "\n")
cat("- Total generations:", ga_with_restart$total_generations, "\n")
cat("- Restart count:", ga_with_restart$restart_count, "\n")

if (ga_with_restart$restart_count > 0) {
  cat("\nRestart History:\n")
  for (i in 1:ga_with_restart$restart_count) {
    restart_info <- ga_with_restart$restart_history[[i]]
    cat(sprintf("  Restart %d: Generation %d, Fitness %.6f\n",
                i, restart_info$generation, restart_info$best_fitness))
  }
}

improvement <- (ga_no_restart$best_fitness - ga_with_restart$best_fitness) / ga_no_restart$best_fitness * 100
cat("\nImprovement with restart:", sprintf("%.2f%%", improvement), "\n")

## ----plot-restart-effect, fig.width=10, fig.height=6--------------------------
# Plot the effect of restart mechanisms
par(mfrow = c(1, 2))

# Fitness evolution comparison
plot(1:nrow(ga_no_restart$fitness_history), ga_no_restart$fitness_history[, "best"],
     type = "l", col = "blue", lwd = 2,
     xlab = "Generation", ylab = "Best Fitness",
     main = "Effect of Restart Mechanism")
lines(1:nrow(ga_with_restart$fitness_history), ga_with_restart$fitness_history[, "best"],
      col = "red", lwd = 2)

# Mark restart points
if (ga_with_restart$restart_count > 0) {
  for (i in 1:ga_with_restart$restart_count) {
    restart_gen <- ga_with_restart$restart_history[[i]]$generation
    abline(v = restart_gen, col = "red", lty = 3, lwd = 2)
  }
}

legend("topright", 
       legend = c("No Restart", "With Restart", "Restart Points"),
       col = c("blue", "red", "red"),
       lty = c(1, 1, 3), lwd = 2)

# Population diversity comparison
no_restart_div <- sapply(ga_no_restart$population_stats, function(x) x$diversity)
with_restart_div <- sapply(ga_with_restart$population_stats, function(x) x$diversity)

plot(1:length(no_restart_div), no_restart_div,
     type = "l", col = "blue", lwd = 2,
     xlab = "Generation", ylab = "Population Diversity",
     main = "Diversity Maintenance with Restarts")
lines(1:length(with_restart_div), with_restart_div, col = "red", lwd = 2)

# Mark restart points for diversity plot
if (ga_with_restart$restart_count > 0) {
  for (i in 1:ga_with_restart$restart_count) {
    restart_gen <- ga_with_restart$restart_history[[i]]$generation
    abline(v = restart_gen, col = "red", lty = 3, lwd = 2)
  }
}

legend("topright", 
       legend = c("No Restart", "With Restart"),
       col = c("blue", "red"), lwd = 2)

par(mfrow = c(1, 1))

## ----selection-methods-comparison---------------------------------------------
cat("Selection Methods Comparison\n")
cat("===========================\n")

# Compare different selection methods
selection_methods <- c("tournament", "elite", "rank", "hybrid")
results <- list()

set.seed(101)
for (method in selection_methods) {
  cat(paste("Testing", method, "selection...\n"))
  
  # Set method-specific parameters
  extra_params <- list()
  if (method == "rank") {
    extra_params$selection_pressure <- 1.4
  }
  
  # Run GA with current selection method
  ga_result <- do.call(subset_ga, c(
    list(
      P = PC_subset,
      Candidates = candidates,
      Test = test_set,
      ntoselect = n_select,
      npop = 25,
      niterations = 40,
      criterion = "pev_mean",
      selection_method = method,
      enable_restart = FALSE,
      verbose = FALSE
    ),
    extra_params
  ))
  
  results[[method]] <- ga_result
}

# Compare results
cat("\nSelection Method Results:\n")
cat("Method      | Best Fitness | Convergence Gen | Final Diversity\n")
cat("------------|--------------|-----------------|----------------\n")

for (method in selection_methods) {
  result <- results[[method]]
  final_div <- result$population_stats[[length(result$population_stats)]]$diversity
  cat(sprintf("%-11s | %12.6f | %15d | %14.4f\n",
              method, result$best_fitness, result$convergence_generation, final_div))
}

