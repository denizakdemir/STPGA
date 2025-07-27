#' STPGA Optimization Criteria Examples
#' 
#' This script demonstrates all the modern optimization criteria functions
#' in the STPGA package with practical examples.

library(STPGA)

# Set seed for reproducibility
set.seed(12345)

cat("=== STPGA Optimization Criteria Examples ===\n\n")

# ============================================================================
# 1. Basic Setup: Create Example Data
# ============================================================================

cat("1. Creating example genomic data...\n")

# Simulate genomic prediction matrix
n_individuals <- 100
n_markers <- 500
P <- matrix(rnorm(n_individuals * n_markers), n_individuals, n_markers)
rownames(P) <- paste0("ind", 1:n_individuals)
colnames(P) <- paste0("snp", 1:n_markers)

# Define training and test sets
train_large <- paste0("ind", 1:60)  # Large training set
train_small <- paste0("ind", 1:30)  # Small training set
test <- paste0("ind", 71:90)        # Test set
validation <- paste0("ind", 91:100) # Validation set

cat(sprintf("- Created prediction matrix: %d individuals x %d markers\n", 
            nrow(P), ncol(P)))
cat(sprintf("- Large training set: %d individuals\n", length(train_large)))
cat(sprintf("- Small training set: %d individuals\n", length(train_small)))
cat(sprintf("- Test set: %d individuals\n", length(test)))

# ============================================================================
# 2. Classical Optimality Criteria
# ============================================================================

cat("\n2. Classical Optimality Criteria\n")
cat("=================================\n")

# A-optimality (minimizes average prediction variance)
aopt_large <- a_optimality(train_large, test, P)
aopt_small <- a_optimality(train_small, test, P)

cat(sprintf("A-optimality (large training): %.6f\n", aopt_large))
cat(sprintf("A-optimality (small training): %.6f\n", aopt_small))
cat(sprintf("Improvement ratio: %.2f\n", aopt_small / aopt_large))

# D-optimality (minimizes generalized variance)
dopt_large <- d_optimality(train_large, test, P)
dopt_small <- d_optimality(train_small, test, P)

cat(sprintf("\nD-optimality (large training): %.6f\n", dopt_large))
cat(sprintf("D-optimality (small training): %.6f\n", dopt_small))
cat(sprintf("Improvement: %.6f\n", dopt_large - dopt_small))

# E-optimality (minimizes maximum variance)
eopt_large <- e_optimality(train_large, test, P)
eopt_small <- e_optimality(train_small, test, P)

cat(sprintf("\nE-optimality (large training): %.6f\n", eopt_large))
cat(sprintf("E-optimality (small training): %.6f\n", eopt_small))
cat(sprintf("Improvement: %.6f\n", eopt_large - eopt_small))

# ============================================================================
# 3. Prediction Error Variance (PEV) Criteria
# ============================================================================

cat("\n3. Prediction Error Variance Criteria\n")
cat("=====================================\n")

# Mean PEV (unnormalized)
pev_mean_large <- pev_mean(train_large, test, P, normalized = FALSE)
pev_mean_small <- pev_mean(train_small, test, P, normalized = FALSE)

cat(sprintf("Mean PEV (large training): %.8f\n", pev_mean_large))
cat(sprintf("Mean PEV (small training): %.8f\n", pev_mean_small))
cat(sprintf("Relative improvement: %.2f%%\n", 
            (pev_mean_small - pev_mean_large) / pev_mean_small * 100))

# Mean PEV (normalized)
pev_norm_large <- pev_mean(train_large, test, P, normalized = TRUE)
pev_norm_small <- pev_mean(train_small, test, P, normalized = TRUE)

cat(sprintf("\nNormalized Mean PEV (large): %.8f\n", pev_norm_large))
cat(sprintf("Normalized Mean PEV (small): %.8f\n", pev_norm_small))

# Maximum PEV
pev_max_large <- pev_max(train_large, test, P)
pev_max_small <- pev_max(train_small, test, P)

cat(sprintf("\nMax PEV (large training): %.8f\n", pev_max_large))
cat(sprintf("Max PEV (small training): %.8f\n", pev_max_small))
cat(sprintf("Max/Mean ratio (large): %.2f\n", pev_max_large / pev_mean_large))
cat(sprintf("Max/Mean ratio (small): %.2f\n", pev_max_small / pev_mean_small))

# ============================================================================
# 4. Cook's Distance Criteria
# ============================================================================

cat("\n4. Cook's Distance Criteria\n")
cat("===========================\n")

# Mean Cook's Distance
cd_mean_large <- cd_mean(train_large, test, P)
cd_mean_small <- cd_mean(train_small, test, P)

cat(sprintf("Mean Cook's Distance (large): %.6f\n", cd_mean_large))
cat(sprintf("Mean Cook's Distance (small): %.6f\n", cd_mean_small))

# Normalized Cook's Distance
cd_norm_large <- cd_mean(train_large, test, P, normalized = TRUE)
cd_norm_small <- cd_mean(train_small, test, P, normalized = TRUE)

cat(sprintf("Normalized CD (large): %.6f\n", cd_norm_large))
cat(sprintf("Normalized CD (small): %.6f\n", cd_norm_small))

# ============================================================================
# 5. Mixed Model Criteria (Genomic Selection)
# ============================================================================

cat("\n5. Mixed Model Criteria (Genomic Selection)\n")
cat("==========================================\n")

# Create kinship matrix (simulated genomic relationships)
K <- matrix(0.05, n_individuals, n_individuals)  # Base relationship
diag(K) <- 1.0  # Self-relationship

# Add some family structure
family1 <- 1:20
family2 <- 21:40
family3 <- 41:60

K[family1, family1] <- 0.4
K[family2, family2] <- 0.4  
K[family3, family3] <- 0.4
diag(K) <- 1.0

rownames(K) <- colnames(K) <- rownames(P)

cat("Created kinship matrix with family structure\n")

# Mixed model PEV with different heritabilities
h2_high <- list(Vg = 0.8, Ve = 0.2)  # High heritability
h2_low <- list(Vg = 0.2, Ve = 0.8)   # Low heritability

pev_mm_high <- pev_mean_mm(train_large, test, P, K, Vg = h2_high$Vg, Ve = h2_high$Ve)
pev_mm_low <- pev_mean_mm(train_large, test, P, K, Vg = h2_low$Vg, Ve = h2_low$Ve)

cat(sprintf("Mixed model PEV (h² = 0.8): %.8f\n", pev_mm_high))
cat(sprintf("Mixed model PEV (h² = 0.2): %.8f\n", pev_mm_low))
cat(sprintf("Heritability effect: %.2fx\n", pev_mm_low / pev_mm_high))

# Mixed model Cook's Distance
cd_mm_high <- cd_mean_mm(train_large, test, P, K, Vg = h2_high$Vg, Ve = h2_high$Ve)
cd_mm_low <- cd_mean_mm(train_large, test, P, K, Vg = h2_low$Vg, Ve = h2_low$Ve)

cat(sprintf("\nMixed model CD (h² = 0.8): %.8f\n", cd_mm_high))
cat(sprintf("Mixed model CD (h² = 0.2): %.8f\n", cd_mm_low))

# ============================================================================
# 6. Unified Criterion Interface
# ============================================================================

cat("\n6. Unified Criterion Interface\n")
cat("==============================\n")

# Compare all criteria using the unified interface
criteria_names <- c("a_optimality", "d_optimality", "e_optimality", 
                   "pev_mean", "pev_max", "cd_mean")

cat("Comparing all criteria for large training set:\n")
for (crit in criteria_names) {
  value <- criterion(train_large, test, P, criterion = crit)
  cat(sprintf("  %-15s: %12.8f\n", crit, value))
}

cat("\nMixed model criteria:\n")
mm_criteria <- c("pev_mean_mm", "cd_mean_mm")
for (crit in mm_criteria) {
  value <- criterion(train_large, test, P, criterion = crit, K = K)
  cat(sprintf("  %-15s: %12.8f\n", crit, value))
}

# ============================================================================
# 7. Legacy Compatibility
# ============================================================================

cat("\n7. Legacy Compatibility\n")
cat("=======================\n")

# Test that legacy names work
legacy_pairs <- list(
  c("AOPT", "a_optimality"),
  c("DOPT", "d_optimality"),
  c("PEVMEAN", "pev_mean"),
  c("PEVMEAN2", "pev_mean_normalized"),
  c("CDMEAN", "cd_mean")
)

cat("Verifying legacy compatibility:\n")
for (pair in legacy_pairs) {
  legacy_val <- criterion(train_large, test, P, criterion = pair[1])
  modern_val <- criterion(train_large, test, P, criterion = pair[2])
  if (pair[2] == "pev_mean_normalized") {
    modern_val <- criterion(train_large, test, P, criterion = "pev_mean", normalized = TRUE)
  }
  match <- abs(legacy_val - modern_val) < 1e-10
  cat(sprintf("  %s == %s: %s\n", pair[1], pair[2], if(match) "✓" else "✗"))
}

# ============================================================================
# 8. Practical Training Set Optimization
# ============================================================================

cat("\n8. Practical Training Set Optimization Example\n")
cat("==============================================\n")

# Simulate optimizing training set selection
candidates <- paste0("ind", 1:70)  # Pool of candidates
n_select <- 30
n_random_sets <- 20

cat(sprintf("Optimizing selection of %d individuals from %d candidates\n", 
            n_select, length(candidates)))

# Generate random training sets
random_results <- data.frame(
  set_id = 1:n_random_sets,
  a_opt = numeric(n_random_sets),
  pev_mean = numeric(n_random_sets),
  cd_mean = numeric(n_random_sets)
)

for (i in 1:n_random_sets) {
  train_random <- sample(candidates, n_select)
  random_results$a_opt[i] <- a_optimality(train_random, test, P)
  random_results$pev_mean[i] <- pev_mean(train_random, test, P)
  random_results$cd_mean[i] <- cd_mean(train_random, test, P)
}

# Show summary statistics
cat("\nRandom training set optimization results:\n")
cat(sprintf("A-optimality:    mean=%.6f, sd=%.6f, range=[%.6f, %.6f]\n",
            mean(random_results$a_opt), sd(random_results$a_opt),
            min(random_results$a_opt), max(random_results$a_opt)))
cat(sprintf("PEV mean:        mean=%.8f, sd=%.8f, range=[%.8f, %.8f]\n",
            mean(random_results$pev_mean), sd(random_results$pev_mean),
            min(random_results$pev_mean), max(random_results$pev_mean)))
cat(sprintf("Cook's Distance: mean=%.6f, sd=%.6f, range=[%.6f, %.6f]\n",
            mean(random_results$cd_mean), sd(random_results$cd_mean),
            min(random_results$cd_mean), max(random_results$cd_mean)))

# Find best set for each criterion
best_aopt_idx <- which.min(random_results$a_opt)
best_pev_idx <- which.min(random_results$pev_mean)

cat(sprintf("\nBest A-optimality: Set #%d (value: %.6f)\n", 
            best_aopt_idx, random_results$a_opt[best_aopt_idx]))
cat(sprintf("Best PEV:          Set #%d (value: %.8f)\n", 
            best_pev_idx, random_results$pev_mean[best_pev_idx]))

# ============================================================================
# 9. Ridge Regularization Effects
# ============================================================================

cat("\n9. Ridge Regularization Effects\n")
cat("===============================\n")

# Test different lambda values
lambda_values <- c(1e-8, 1e-6, 1e-4, 1e-2, 1e-1)
cat("Testing ridge regularization effects:\n")
cat("Lambda\t\tA-optimality\tPEV Mean\n")
cat("------\t\t------------\t--------\n")

for (lambda in lambda_values) {
  aopt_ridge <- a_optimality(train_large, test, P, lambda = lambda)
  pev_ridge <- pev_mean(train_large, test, P, lambda = lambda)
  cat(sprintf("%.0e\t\t%.6f\t%.8f\n", lambda, aopt_ridge, pev_ridge))
}

cat("\n=== Example Complete ===\n")
cat("All STPGA optimization criteria have been demonstrated!\n")
cat("See individual function documentation for more details.\n")