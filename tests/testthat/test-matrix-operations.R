test_that("compute_amatrix function works correctly", {
  # Create test marker matrix
  set.seed(123)
  n_individuals <- 20
  n_markers <- 50
  
  # Create marker matrix coded as -1, 0, 1
  M <- matrix(sample(c(-1, 0, 1), n_individuals * n_markers, replace = TRUE), 
              n_individuals, n_markers)
  rownames(M) <- paste0("ind", 1:n_individuals)
  colnames(M) <- paste0("marker", 1:n_markers)
  
  # Test basic functionality
  A <- compute_amatrix(M)
  
  expect_type(A, "double")
  expect_equal(dim(A), c(n_individuals, n_individuals))
  expect_equal(rownames(A), rownames(M))
  expect_equal(colnames(A), rownames(M))
  
  # Test symmetry
  expect_true(all(abs(A - t(A)) < 1e-10))
  
  # Test diagonal elements (should be around 1)
  diag_vals <- diag(A)
  expect_true(all(diag_vals > 0))
  expect_true(mean(abs(diag_vals - 1)) < 0.5)  # Should be close to 1
  
  # Test with different number of pieces
  A_pieces5 <- compute_amatrix(M, pieces = 5)
  expect_equal(dim(A_pieces5), dim(A))
  
  # Results should be similar (but not identical due to chunking)
  correlation <- cor(as.vector(A), as.vector(A_pieces5))
  expect_true(correlation > 0.95)
  
  # Test with parallel processing (mc.cores = 1 for testing)
  A_parallel <- compute_amatrix(M, pieces = 3, mc.cores = 1)
  expect_equal(dim(A_parallel), dim(A))
  
  # Test error conditions
  expect_error(compute_amatrix("not_a_matrix"))
  expect_error(compute_amatrix(M, pieces = 0))
  expect_error(compute_amatrix(M, pieces = n_markers + 1))
})

test_that("genomic_relationship_matrix function works correctly", {
  # Create test marker matrix
  set.seed(123)
  n_individuals <- 15
  n_markers <- 30
  
  M <- matrix(sample(c(-1, 0, 1), n_individuals * n_markers, replace = TRUE),
              n_individuals, n_markers)
  rownames(M) <- paste0("ind", 1:n_individuals)
  
  # Test VanRaden method
  G_vanraden <- genomic_relationship_matrix(M, method = "vanraden")
  
  expect_type(G_vanraden, "double")
  expect_equal(dim(G_vanraden), c(n_individuals, n_individuals))
  expect_equal(rownames(G_vanraden), rownames(M))
  expect_equal(colnames(G_vanraden), rownames(M))
  expect_true(all(abs(G_vanraden - t(G_vanraden)) < 1e-10))  # Symmetry
  
  # Test pieces method
  G_pieces <- genomic_relationship_matrix(M, method = "pieces", pieces = 5)
  expect_equal(dim(G_pieces), dim(G_vanraden))
  
  # Test standard method
  G_standard <- genomic_relationship_matrix(M, method = "standard")
  expect_equal(dim(G_standard), dim(G_vanraden))
  
  # Test that methods give similar results
  cor_van_pieces <- cor(as.vector(G_vanraden), as.vector(G_pieces))
  expect_true(cor_van_pieces > 0.9)
  
  # Test with parallel processing
  G_parallel <- genomic_relationship_matrix(M, method = "pieces", 
                                           pieces = 3, mc.cores = 1)
  expect_equal(dim(G_parallel), dim(G_vanraden))
  
  # Test invalid method
  expect_error(genomic_relationship_matrix(M, method = "invalid"))
})

test_that("safe_matrix_inverse function works correctly", {
  set.seed(123)
  
  # Test with well-conditioned matrix
  n <- 10
  X <- matrix(rnorm(n^2), n, n)
  X <- X %*% t(X) + diag(n)  # Make positive definite
  
  # Test Cholesky method
  inv_chol <- safe_matrix_inverse(X, method = "cholesky")
  expect_equal(dim(inv_chol), c(n, n))
  
  # Check that X * inv(X) ≈ I
  product <- X %*% inv_chol
  identity_approx <- abs(product - diag(n))
  expect_true(max(identity_approx) < 1e-10)
  
  # Test SVD method
  inv_svd <- safe_matrix_inverse(X, method = "svd")
  expect_equal(dim(inv_svd), c(n, n))
  
  # Test eigen method
  inv_eigen <- safe_matrix_inverse(X, method = "eigen")
  expect_equal(dim(inv_eigen), c(n, n))
  
  # Results should be similar
  expect_true(max(abs(inv_chol - inv_svd)) < 1e-8)
  expect_true(max(abs(inv_chol - inv_eigen)) < 1e-8)
  
  # Test with ridge regularization
  inv_ridge <- safe_matrix_inverse(X, lambda = 0.1)
  expect_equal(dim(inv_ridge), c(n, n))
  
  # Test with ill-conditioned matrix
  X_ill <- matrix(1, 5, 5)  # Rank 1 matrix
  X_ill <- X_ill + 1e-10 * diag(5)  # Slightly regularized
  
  expect_warning(inv_ill_chol <- safe_matrix_inverse(X_ill, method = "cholesky"))
  inv_ill_svd <- safe_matrix_inverse(X_ill, method = "svd")
  
  expect_equal(dim(inv_ill_svd), c(5, 5))
  expect_true(all(is.finite(inv_ill_svd)))
  
  # Test invalid method
  expect_error(safe_matrix_inverse(X, method = "invalid"))
})

test_that("matrix_stability_check function works correctly", {
  set.seed(123)
  
  # Test with well-conditioned symmetric matrix
  n <- 8
  X_good <- matrix(rnorm(n^2), n, n)
  X_good <- X_good %*% t(X_good) + diag(n)  # Positive definite
  
  stability_good <- matrix_stability_check(X_good)
  
  expect_type(stability_good, "list")
  expect_true("condition_number" %in% names(stability_good))
  expect_true("is_positive_definite" %in% names(stability_good))
  expect_true("is_well_conditioned" %in% names(stability_good))
  expect_true("rank_deficient" %in% names(stability_good))
  expect_true("numerical_rank" %in% names(stability_good))
  
  expect_true(stability_good$condition_number > 1)
  expect_true(stability_good$is_positive_definite)
  expect_true(stability_good$is_well_conditioned)
  expect_false(stability_good$rank_deficient)
  expect_equal(stability_good$numerical_rank, n)
  
  # Test with ill-conditioned matrix
  X_ill <- matrix(1, 4, 4)  # Rank 1
  X_ill <- X_ill + 1e-15 * diag(4)  # Tiny regularization
  
  stability_ill <- matrix_stability_check(X_ill)
  expect_true(stability_ill$condition_number > 1e10)
  expect_false(stability_ill$is_well_conditioned)
  expect_true(stability_ill$rank_deficient)
  expect_true(stability_ill$numerical_rank < 4)
  
  # Test with non-symmetric matrix
  X_nonsym <- matrix(rnorm(12), 3, 4)
  stability_nonsym <- matrix_stability_check(X_nonsym)
  expect_true(is.na(stability_nonsym$is_positive_definite))
  expect_true(is.finite(stability_nonsym$condition_number))
  
  # Test singular matrix
  X_singular <- matrix(0, 3, 3)
  stability_singular <- matrix_stability_check(X_singular)
  expect_equal(stability_singular$numerical_rank, 0)
  expect_true(stability_singular$rank_deficient)
})

test_that("ridge_regression_cv function works correctly", {
  set.seed(123)
  
  # Create test data
  n <- 50
  p <- 8
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- rnorm(p)
  y <- X %*% beta_true + rnorm(n, sd = 0.5)
  
  # Test basic functionality
  ridge_result <- ridge_regression_cv(y, X, cv_folds = 3)
  
  expect_type(ridge_result, "list")
  expect_true("lambda_optimal" %in% names(ridge_result))
  expect_true("beta" %in% names(ridge_result))
  expect_true("lambda_seq" %in% names(ridge_result))
  expect_true("cv_errors" %in% names(ridge_result))
  
  expect_type(ridge_result$lambda_optimal, "double")
  expect_length(ridge_result$lambda_optimal, 1)
  expect_true(ridge_result$lambda_optimal > 0)
  
  expect_type(ridge_result$beta, "double")
  expect_length(ridge_result$beta, p)
  
  # Test with custom lambda sequence
  lambda_seq_custom <- c(0.001, 0.01, 0.1, 1)
  ridge_custom <- ridge_regression_cv(y, X, lambda_seq = lambda_seq_custom, cv_folds = 3)
  expect_true(ridge_custom$lambda_optimal %in% lambda_seq_custom)
  
  # Test with AIC criterion
  ridge_aic <- ridge_regression_cv(y, X, criterion = "aic", cv_folds = 1)
  expect_true(is.finite(ridge_aic$lambda_optimal))
  
  # Test with BIC criterion
  ridge_bic <- ridge_regression_cv(y, X, criterion = "bic", cv_folds = 1)
  expect_true(is.finite(ridge_bic$lambda_optimal))
  
  # Test prediction accuracy
  y_pred <- X %*% ridge_result$beta
  correlation <- cor(y, y_pred)
  expect_true(correlation > 0.5)  # Should have reasonable prediction
  
  # Test with high-dimensional data (p > n)
  n_small <- 20
  p_large <- 30
  X_hd <- matrix(rnorm(n_small * p_large), n_small, p_large)
  y_hd <- rnorm(n_small)
  
  ridge_hd <- ridge_regression_cv(y_hd, X_hd, cv_folds = 3)
  expect_length(ridge_hd$beta, p_large)
  expect_true(is.finite(ridge_hd$lambda_optimal))
  
  # Test invalid criterion
  expect_error(ridge_regression_cv(y, X, criterion = "invalid"))
})

test_that("legacy Amat.pieces wrapper works", {
  # Create test data
  set.seed(123)
  n <- 15
  p <- 25
  M <- matrix(sample(c(-1, 0, 1), n * p, replace = TRUE), n, p)
  rownames(M) <- paste0("ind", 1:n)
  
  # Test that legacy function matches new function
  A_new <- compute_amatrix(M, pieces = 5)
  A_legacy <- Amat.pieces(M, pieces = 5)
  
  expect_equal(A_new, A_legacy)
  
  # Test with different parameters
  A_legacy_par <- Amat.pieces(M, pieces = 3, mc.cores = 1)
  expect_equal(dim(A_legacy_par), c(n, n))
})

test_that("numerical stability and edge cases", {
  # Test with very small matrix
  M_small <- matrix(c(-1, 1, 0, 1), 2, 2)
  rownames(M_small) <- c("ind1", "ind2")
  
  A_small <- compute_amatrix(M_small)
  expect_equal(dim(A_small), c(2, 2))
  expect_true(all(is.finite(A_small)))
  
  # Test with constant markers (zero variance)
  M_constant <- matrix(1, 10, 5)  # All markers are constant
  A_constant <- compute_amatrix(M_constant)
  expect_equal(dim(A_constant), c(10, 10))
  expect_true(all(is.finite(A_constant)))
  
  # Test genomic relationship with extreme values
  M_extreme <- matrix(c(rep(-1, 50), rep(1, 50)), 10, 10)
  G_extreme <- genomic_relationship_matrix(M_extreme, method = "vanraden")
  expect_true(all(is.finite(G_extreme)))
  
  # Test matrix inversion with near-singular matrix
  X_near_singular <- diag(5) + 1e-12
  inv_near_singular <- safe_matrix_inverse(X_near_singular, lambda = 1e-10)
  expect_true(all(is.finite(inv_near_singular)))
  
  # Test ridge regression with perfect collinearity
  X_collinear <- cbind(1:20, 2 * (1:20), 3 * (1:20))  # Perfectly collinear
  y_collinear <- rnorm(20)
  
  ridge_collinear <- ridge_regression_cv(y_collinear, X_collinear, cv_folds = 3)
  expect_true(is.finite(ridge_collinear$lambda_optimal))
  expect_true(ridge_collinear$lambda_optimal > 1e-8)  # Should use regularization
})

test_that("parallel processing works correctly", {
  # Note: Testing with mc.cores = 1 to avoid issues in test environment
  set.seed(123)
  n <- 20
  p <- 40
  M <- matrix(sample(c(-1, 0, 1), n * p, replace = TRUE), n, p)
  rownames(M) <- paste0("ind", 1:n)
  
  # Test that parallel and sequential give same results
  A_sequential <- compute_amatrix(M, pieces = 4, mc.cores = 1)
  A_parallel <- compute_amatrix(M, pieces = 4, mc.cores = 1)  # Keep as 1 for testing
  
  expect_equal(A_sequential, A_parallel)
  
  # Test genomic relationship matrix with parallel
  G_sequential <- genomic_relationship_matrix(M, method = "pieces", pieces = 3, mc.cores = 1)
  G_parallel <- genomic_relationship_matrix(M, method = "pieces", pieces = 3, mc.cores = 1)
  
  expect_equal(G_sequential, G_parallel)
})

test_that("memory efficiency and large matrices", {
  # Test that function can handle moderately large matrices
  set.seed(123)
  n <- 50
  p <- 100
  
  # Don't actually create huge matrices in tests, but test the logic
  # Create smaller version to test chunking logic
  M <- matrix(sample(c(-1, 0, 1), n * p, replace = TRUE), n, p)
  rownames(M) <- paste0("ind", 1:n)
  
  # Test with many pieces (should still work)
  A_many_pieces <- compute_amatrix(M, pieces = 20)
  expect_equal(dim(A_many_pieces), c(n, n))
  expect_true(all(is.finite(A_many_pieces)))
  
  # Test that chunking doesn't break symmetry
  expect_true(all(abs(A_many_pieces - t(A_many_pieces)) < 1e-10))
})