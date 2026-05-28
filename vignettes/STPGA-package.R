## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----load-package-------------------------------------------------------------
library(STPGA)

# Load the wheat dataset
data("WheatData")

# Explore the data structure
cat("Dataset components:\n")
cat("- Wheat.Y: Phenotype data (", nrow(Wheat.Y), " observations)\n")
cat("- Wheat.M: Marker data (", nrow(Wheat.M), " individuals x ", ncol(Wheat.M), " markers)\n")
cat("- Wheat.K: Kinship matrix (", nrow(Wheat.K), " x ", ncol(Wheat.K), ")\n")

## ----explore-data-------------------------------------------------------------
# Look at the phenotype data
head(Wheat.Y)

# Check marker data dimensions and first few markers
cat("Marker matrix preview:\n")
Wheat.M[1:5, 1:5]

# Examine kinship matrix structure
cat("Kinship matrix diagonal (first 10):\n")
diag(Wheat.K)[1:10]
cat("Kinship matrix off-diagonal (first 5x5):\n")
Wheat.K[1:5, 1:5]

## ----setup-demo---------------------------------------------------------------
# Define candidate pool and test sets
all_individuals <- rownames(Wheat.M)
n_total <- length(all_individuals)

# Create training candidates and test set
set.seed(123)
test_set <- sample(all_individuals, 40)
candidates <- setdiff(all_individuals, test_set)

cat("Total individuals:", n_total, "\n")
cat("Test set size:", length(test_set), "\n")
cat("Candidate pool size:", length(candidates), "\n")

# Extract first 5 principal components for non-mixed model criteria
# This reduces computational complexity for high-dimensional marker data
cat("\nExtracting first 5 principal components for non-mixed model criteria...\n")
pca_result <- prcomp(Wheat.M, center = TRUE, scale. = TRUE)
Wheat.PC5 <- pca_result$x[, 1:5]
rownames(Wheat.PC5) <- rownames(Wheat.M)

cat("PC matrix dimensions:", nrow(Wheat.PC5), "x", ncol(Wheat.PC5), "\n")
cat("Variance explained by first 5 PCs:", sprintf("%.2f%%", 100 * sum(pca_result$sdev[1:5]^2) / sum(pca_result$sdev^2)), "\n")

## ----classical-criteria-------------------------------------------------------
# Select training sets of different sizes
train_small <- sample(candidates, 30)
train_medium <- sample(candidates, 60)
train_large <- sample(candidates, 100)

# Compute A-optimality using first 5 PCs (minimizes average prediction variance)
aopt_small <- a_optimality(train_small, test_set, Wheat.PC5)
aopt_medium <- a_optimality(train_medium, test_set, Wheat.PC5)
aopt_large <- a_optimality(train_large, test_set, Wheat.PC5)

cat("A-optimality results:\n")
cat("Small training (n=30):", sprintf("%.2e", aopt_small), "\n")
cat("Medium training (n=60):", sprintf("%.2e", aopt_medium), "\n")
cat("Large training (n=100):", sprintf("%.2e", aopt_large), "\n")
cat("Improvement (large vs small):", sprintf("%.1fx", aopt_small/aopt_large), "\n")

# D-optimality using first 5 PCs (minimizes generalized variance)
dopt_small <- d_optimality(train_small, test_set, Wheat.PC5)
dopt_medium <- d_optimality(train_medium, test_set, Wheat.PC5)
dopt_large <- d_optimality(train_large, test_set, Wheat.PC5)

cat("\nD-optimality results:\n")
cat("Small training (n=30):", sprintf("%.2f", dopt_small), "\n")
cat("Medium training (n=60):", sprintf("%.2f", dopt_medium), "\n")
cat("Large training (n=100):", sprintf("%.2f", dopt_large), "\n")

# E-optimality using first 5 PCs (minimizes maximum variance)
eopt_small <- e_optimality(train_small, test_set, Wheat.PC5)
eopt_medium <- e_optimality(train_medium, test_set, Wheat.PC5)
eopt_large <- e_optimality(train_large, test_set, Wheat.PC5)

cat("\nE-optimality results:\n")
cat("Small training (n=30):", sprintf("%.2f", eopt_small), "\n")
cat("Medium training (n=60):", sprintf("%.2f", eopt_medium), "\n")
cat("Large training (n=100):", sprintf("%.2f", eopt_large), "\n")

## ----pev-criteria-------------------------------------------------------------
# Mean PEV using first 5 PCs
pev_mean_small <- pev_mean(train_small, test_set, Wheat.PC5)
pev_mean_large <- pev_mean(train_large, test_set, Wheat.PC5)

# Normalized PEV using first 5 PCs (scale-invariant)
pev_norm_small <- pev_mean(train_small, test_set, Wheat.PC5, normalized = TRUE)
pev_norm_large <- pev_mean(train_large, test_set, Wheat.PC5, normalized = TRUE)

# Maximum PEV using first 5 PCs (worst-case scenario)
pev_max_small <- pev_max(train_small, test_set, Wheat.PC5)
pev_max_large <- pev_max(train_large, test_set, Wheat.PC5)

cat("PEV Results:\n")
cat("Mean PEV (small):", sprintf("%.2e", pev_mean_small), "\n")
cat("Mean PEV (large):", sprintf("%.2e", pev_mean_large), "\n")
cat("Improvement:", sprintf("%.1fx", pev_mean_small/pev_mean_large), "\n")

cat("\nNormalized PEV:\n")
cat("Small training:", sprintf("%.6f", pev_norm_small), "\n")
cat("Large training:", sprintf("%.6f", pev_norm_large), "\n")

cat("\nMax PEV:\n")
cat("Small training:", sprintf("%.2e", pev_max_small), "\n")
cat("Large training:", sprintf("%.2e", pev_max_large), "\n")
cat("Max/Mean ratio (small):", sprintf("%.2f", pev_max_small/pev_mean_small), "\n")
cat("Max/Mean ratio (large):", sprintf("%.2f", pev_max_large/pev_mean_large), "\n")

## ----cook-distance------------------------------------------------------------
# Mean Coefficient of Determination (R²) using first 5 PCs
cd_mean_small <- cd_mean(train_small, test_set, Wheat.PC5)
cd_mean_large <- cd_mean(train_large, test_set, Wheat.PC5)

# Normalized Coefficient of Determination using first 5 PCs
cd_norm_small <- cd_mean(train_small, test_set, Wheat.PC5, normalized = TRUE)
cd_norm_large <- cd_mean(train_large, test_set, Wheat.PC5, normalized = TRUE)

cat("Coefficient of Determination Results:\n")
cat("Mean CD (small):", sprintf("%.2e", cd_mean_small), "\n")
cat("Mean CD (large):", sprintf("%.2e", cd_mean_large), "\n")

cat("\nNormalized R²:\n")
cat("Small training:", sprintf("%.6f", cd_norm_small), "\n")
cat("Large training:", sprintf("%.6f", cd_norm_large), "\n")

## ----mixed-models-------------------------------------------------------------
# Select smaller sets for computational efficiency with mixed models
train_mm <- sample(candidates, 50)
test_mm <- sample(test_set, 20)

# Extract relevant kinship submatrices
K_train <- Wheat.K[train_mm, train_mm]
K_test <- Wheat.K[test_mm, test_mm]
K_test_train <- Wheat.K[test_mm, train_mm]

# Mixed model PEV with different heritability scenarios
h2_high <- list(Vg = 0.8, Ve = 0.2)  # High heritability (80%)
h2_low <- list(Vg = 0.3, Ve = 0.7)   # Low heritability (30%)

# Use smaller dataset for computational efficiency
# For mixed models, use subset of individuals but all relationships
all_mm_individuals <- c(train_mm, test_mm)
P_subset <- Wheat.M[all_mm_individuals, 1:100]  # Use only first 100 markers
K_subset <- Wheat.K[all_mm_individuals, all_mm_individuals]

# Compute mixed model criteria with smaller dataset
pev_mm_high <- pev_mean_mm(train_mm, test_mm, P_subset, K_subset, 
                           Vg = h2_high$Vg, Ve = h2_high$Ve)
pev_mm_low <- pev_mean_mm(train_mm, test_mm, P_subset, K_subset, 
                          Vg = h2_low$Vg, Ve = h2_low$Ve)

cd_mm_high <- cd_mean_mm(train_mm, test_mm, P_subset, K_subset, 
                         Vg = h2_high$Vg, Ve = h2_high$Ve)
cd_mm_low <- cd_mean_mm(train_mm, test_mm, P_subset, K_subset, 
                        Vg = h2_low$Vg, Ve = h2_low$Ve)

cat("Mixed Model Results:\n")
cat("PEV (h² = 0.8):", sprintf("%.6f", pev_mm_high), "\n")
cat("PEV (h² = 0.3):", sprintf("%.6f", pev_mm_low), "\n")
cat("Heritability effect:", sprintf("%.2fx", pev_mm_low/pev_mm_high), "\n")

cat("\nCoefficient of Determination (mixed model):\n")
cat("CD (h² = 0.8):", sprintf("%.6f", cd_mm_high), "\n")
cat("CD (h² = 0.3):", sprintf("%.6f", cd_mm_low), "\n")

## ----unified-interface--------------------------------------------------------
# Define a medium-sized training set
train_unified <- sample(candidates, 60)

# Test all modern criteria
modern_criteria <- c("a_optimality", "d_optimality", "e_optimality", 
                     "pev_mean", "pev_max", "cd_mean")

cat("Modern Criteria Results:\n")
cat("========================\n")
for (crit in modern_criteria) {
  value <- criterion(train_unified, test_set, Wheat.PC5, criterion = crit)
  cat(sprintf("%-15s: %12.6e\n", crit, value))
}

# Test normalized versions
cat("\nNormalized Criteria:\n")
cat("===================\n")
pev_norm <- criterion(train_unified, test_set, Wheat.PC5, criterion = "pev_mean_normalized")
cd_norm <- criterion(train_unified, test_set, Wheat.PC5, criterion = "cd_mean_normalized")
cat(sprintf("%-20s: %12.6f\n", "pev_mean_normalized", pev_norm))
cat(sprintf("%-20s: %12.6f\n", "cd_mean_normalized", cd_norm))

# Test mixed model criteria (using subset for speed)
cat("\nMixed Model Criteria (100 markers):\n")
cat("===================================\n")
P_unified_subset <- Wheat.M[c(train_unified, test_set), 1:100]
K_unified_subset <- Wheat.K[c(train_unified, test_set), c(train_unified, test_set)]
pev_mm <- criterion(train_unified, test_set, P_unified_subset, criterion = "pev_mean_mm", K = K_unified_subset)
cd_mm <- criterion(train_unified, test_set, P_unified_subset, criterion = "cd_mean_mm", K = K_unified_subset)
cat(sprintf("%-15s: %12.6e\n", "pev_mean_mm", pev_mm))
cat(sprintf("%-15s: %12.6e\n", "cd_mean_mm", cd_mm))

# Demonstrate legacy compatibility
cat("\nLegacy Compatibility:\n")
cat("====================\n")
aopt_modern <- criterion(train_unified, test_set, Wheat.PC5, criterion = "a_optimality")
aopt_legacy <- criterion(train_unified, test_set, Wheat.PC5, criterion = "AOPT")
cat("Modern a_optimality:", sprintf("%.6e", aopt_modern), "\n")
cat("Legacy AOPT:        ", sprintf("%.6e", aopt_legacy), "\n")
cat("Values match:", aopt_modern == aopt_legacy, "\n")

## ----single-objective-ga------------------------------------------------------
# Define optimization parameters
n_select <- 40  # Number of individuals to select
n_candidates <- length(candidates)

cat("Enhanced Genetic Algorithm Optimization\n")
cat("======================================\n")
cat("Candidate pool size:", n_candidates, "\n")
cat("Selection target:", n_select, "\n")
cat("Test set size:", length(test_set), "\n")

# Demonstrate the genetic algorithm with advanced features
set.seed(456)

# Basic genetic algorithm run
cat("\nRunning genetic algorithm with default settings...\n")
ga_basic <- subset_ga(
  P = Wheat.PC5,
  Candidates = candidates,
  Test = test_set,
  ntoselect = n_select,
  npop = 30,
  niterations = 50,
  criterion = "pev_mean",
  verbose = FALSE
)

cat("Basic GA Results:\n")
cat("Best fitness:", sprintf("%.8f", ga_basic$best_fitness), "\n")
cat("Convergence generation:", ga_basic$convergence_generation, "\n")
cat("Selection method:", ga_basic$parameters$selection_method, "\n")

# Enhanced genetic algorithm with rank-based selection and restart
cat("\nRunning enhanced GA with rank-based selection and restart...\n")
ga_enhanced <- subset_ga(
  P = Wheat.PC5,
  Candidates = candidates,
  Test = test_set,
  ntoselect = n_select,
  npop = 30,
  niterations = 50,
  criterion = "pev_mean",
  selection_method = "rank",
  selection_pressure = 1.4,
  enable_restart = TRUE,
  restart_threshold = 0.6,
  max_restarts = 1,
  convergence_window_multiplier = 3,
  verbose = FALSE
)

cat("Enhanced GA Results:\n")
cat("Best fitness:", sprintf("%.8f", ga_enhanced$best_fitness), "\n")
cat("Convergence generation:", ga_enhanced$convergence_generation, "\n")
cat("Total generations:", ga_enhanced$total_generations, "\n")
cat("Restart count:", ga_enhanced$restart_count, "\n")
cat("Selection method:", ga_enhanced$parameters$selection_method, "\n")
cat("Selection pressure:", ga_enhanced$parameters$selection_pressure, "\n")

# Compare with random selection
random_selection <- sample(candidates, n_select)
random_fitness <- pev_mean(random_selection, test_set, Wheat.PC5, normalized = TRUE)

cat("\nComparison Results:\n")
cat("Random selection PEV:", sprintf("%.8f", random_fitness), "\n")
cat("Basic GA PEV:       ", sprintf("%.8f", ga_basic$best_fitness), "\n")
cat("Enhanced GA PEV:    ", sprintf("%.8f", ga_enhanced$best_fitness), "\n")

basic_improvement <- (random_fitness - ga_basic$best_fitness) / random_fitness * 100
enhanced_improvement <- (random_fitness - ga_enhanced$best_fitness) / random_fitness * 100

cat("Basic GA improvement:", sprintf("%.2f%%", basic_improvement), "\n")
cat("Enhanced GA improvement:", sprintf("%.2f%%", enhanced_improvement), "\n")

