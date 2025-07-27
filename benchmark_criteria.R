# Benchmark current criteria performance
library(STPGA)

# Create test data
set.seed(123)
n_individuals <- 1000
n_markers <- 500
P_full <- matrix(rnorm(n_individuals * n_markers), n_individuals, n_markers)
rownames(P_full) <- paste0("ind", 1:n_individuals)

# Extract first 5 principal components for non-mixed model criteria
# This reduces computational complexity for high-dimensional marker data
cat("Extracting first 5 principal components...\n")
pca_result <- prcomp(P_full, center = TRUE, scale. = TRUE)
P <- pca_result$x[, 1:5]
rownames(P) <- rownames(P_full)

cat("Original matrix:", nrow(P_full), "x", ncol(P_full), "\n")
cat("PCA matrix:", nrow(P), "x", ncol(P), "\n")
cat("Variance explained by first 5 PCs:", sprintf("%.2f%%", 100 * sum(pca_result$sdev[1:5]^2) / sum(pca_result$sdev^2)), "\n\n")

# Create kinship matrix
K <- matrix(0.1, n_individuals, n_individuals) + diag(0.9, n_individuals)
rownames(K) <- colnames(K) <- rownames(P)

# Test different training set sizes
train_sizes <- c(50, 100, 200, 500)

benchmark_results <- list()

for (size in train_sizes) {
  cat("Testing training set size:", size, "\n")
  
  train <- paste0("ind", 1:size)
  test <- paste0("ind", (size+1):(size+100))
  
  # Benchmark different criteria
  times <- list()
  
  # A-optimality (using first 5 PCs)
  times$a_opt <- system.time({
    a_val <- a_optimality(train, test, P)
  })[3]
  
  # PEV mean (using first 5 PCs)
  times$pev_mean <- system.time({
    pev_val <- pev_mean(train, test, P)
  })[3]
  
  # Mixed model PEV (most expensive - uses full marker data)
  times$pev_mm <- system.time({
    pev_mm_val <- pev_mean_mm(train, test, P_full, K)
  })[3]
  
  # CD mixed model (uses full marker data)
  times$cd_mm <- system.time({
    cd_mm_val <- cd_mean_mm(train, test, P_full, K)
  })[3]
  
  benchmark_results[[as.character(size)]] <- times
  
  cat("Times - A-opt (PC5):", times$a_opt, "PEV (PC5):", times$pev_mean, 
      "PEV-MM (full):", times$pev_mm, "CD-MM (full):", times$cd_mm, "\n\n")
}

# Print summary
cat("=== PERFORMANCE SUMMARY ===\n")
cat("Note: Non-mixed model criteria use first 5 PCs for efficiency\n")
cat("Mixed model criteria use full marker matrix as required\n\n")
for (size in names(benchmark_results)) {
  times <- benchmark_results[[size]]
  cat("Size", size, "- Total time:", sum(unlist(times)), "seconds\n")
  cat("  Speedup from PCA (non-mixed):", sprintf("%.1fx faster", ncol(P_full)/ncol(P)), "\n")
}