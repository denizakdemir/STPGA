# Test safe matrix operations
library(STPGA)

test_that("safe_matrix_inverse works correctly", {
  # Well-conditioned matrix
  X <- matrix(c(2, 1, 1, 2), nrow = 2)
  X_inv <- safe_matrix_inverse(X)
  
  # Check if inverse is correct
  expect_true(all(abs(X %*% X_inv - diag(2)) < 1e-10))
  
  # Positive definite matrix (Cholesky should work)
  A <- matrix(rnorm(25), 5, 5)
  A <- t(A) %*% A + diag(5)  # Make positive definite
  A_inv <- safe_matrix_inverse(A)
  expect_true(all(abs(A %*% A_inv - diag(5)) < 1e-10))
  
  # Test ridge regularization
  X_reg <- safe_matrix_inverse(X, lambda = 0.1)
  expect_true(is.matrix(X_reg))
  expect_equal(nrow(X_reg), ncol(X_reg))
  
  # Test error handling
  expect_error(safe_matrix_inverse("not_matrix"), "must be a matrix")
  expect_error(safe_matrix_inverse(matrix(1:6, 2, 3)), "square")
  
  # Test with NA values
  X_na <- X
  X_na[1, 1] <- NA
  expect_error(safe_matrix_inverse(X_na), "non-finite")
})

test_that("safe_matrix_inverse handles rank deficient matrices", {
  # Create rank deficient matrix
  X <- matrix(c(1, 2, 2, 4), nrow = 2)  # rank 1
  
  # Should warn about rank deficiency and use pseudo-inverse
  expect_warning(X_inv <- safe_matrix_inverse(X), "rank deficient")
  
  # Check pseudo-inverse property: A * A+ * A = A
  expect_true(all(abs(X %*% X_inv %*% X - X) < 1e-10))
})

test_that("compute_prediction_core works correctly", {
  # Create test data
  set.seed(123)
  P <- matrix(rnorm(100), nrow = 20, ncol = 5)
  rownames(P) <- paste0("ind", 1:20)
  
  p_train <- P[1:15, ]
  p_test <- P[16:20, ]
  
  # Test basic functionality
  core <- compute_prediction_core(p_train, p_test)
  
  expect_true(is.list(core))
  expect_true("inv_xtx_reg" %in% names(core))
  expect_true("pev_matrix" %in% names(core))
  expect_equal(nrow(core$pev_matrix), nrow(p_test))
  expect_equal(ncol(core$pev_matrix), nrow(p_test))
  
  # Test with training set only
  core_train <- compute_prediction_core(p_train)
  expect_equal(nrow(core_train$pev_matrix), nrow(p_train))
  
  # Test with contrast matrix
  C <- matrix(c(1, -1, 0, 0, 0), nrow = 1)
  core_contrast <- compute_prediction_core(p_train, p_test, C = C)
  expect_true(!is.null(core_contrast$contrast_core))
  
  # Test dimension validation
  p_test_bad <- matrix(rnorm(15), nrow = 5, ncol = 3)  # Wrong number of columns
  expect_error(compute_prediction_core(p_train, p_test_bad), "same number of columns")
  
  # Test minimum size requirements
  p_tiny <- matrix(1, nrow = 1, ncol = 5)
  expect_error(compute_prediction_core(p_tiny), "at least 2 individuals")
})

test_that("compute_prediction_core handles adaptive ridge", {
  # Create ill-conditioned data
  set.seed(456)
  P <- matrix(rnorm(50), nrow = 10, ncol = 5)
  P[, 5] <- P[, 1] + 1e-12 * rnorm(10)  # Near collinearity
  
  # Should use adaptive ridge
  expect_message(core <- compute_prediction_core(P, adaptive_ridge = TRUE), 
                 "adaptive ridge")
  
  # Should not use adaptive ridge when disabled
  expect_silent(core_no_adapt <- compute_prediction_core(P, adaptive_ridge = FALSE))
})

test_that("h2_to_variances works correctly", {
  # Test basic conversion
  result <- h2_to_variances(0.8)
  expect_equal(result$Vg, 0.8)
  expect_equal(result$Ve, 0.2)
  expect_equal(result$Vg + result$Ve, 1)
  
  # Test with different total variance
  result2 <- h2_to_variances(0.6, total_var = 2)
  expect_equal(result2$Vg, 1.2)
  expect_equal(result2$Ve, 0.8)
  expect_equal(result2$Vg + result2$Ve, 2)
  
  # Test edge cases
  result_zero <- h2_to_variances(0)
  expect_equal(result_zero$Vg, 0)
  expect_equal(result_zero$Ve, 1)
  
  result_one <- h2_to_variances(1)
  expect_equal(result_one$Vg, 1)
  expect_equal(result_one$Ve, 0)
})