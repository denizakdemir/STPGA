# Test validation functions
library(STPGA)

test_that("validate_matrix_params works correctly", {
  # Create test data
  P <- matrix(rnorm(50), nrow = 10, ncol = 5)
  rownames(P) <- paste0("ind", 1:10)
  colnames(P) <- paste0("var", 1:5)
  
  train <- c("ind1", "ind2", "ind3")
  test <- c("ind4", "ind5")
  
  # Test valid input
  expect_silent(validate_matrix_params(P, train = train, test = test))
  
  # Test invalid P matrix
  expect_error(validate_matrix_params("not_a_matrix"), "P must be a matrix")
  
  P_bad <- P
  P_bad[1, 1] <- NA
  expect_error(validate_matrix_params(P_bad), "non-finite values")
  
  # Test missing row names
  P_no_names <- P
  rownames(P_no_names) <- NULL
  expect_error(validate_matrix_params(P_no_names), "row names")
  
  # Test invalid lambda
  expect_error(validate_matrix_params(P, lambda = -1), "non-negative")
  expect_error(validate_matrix_params(P, lambda = c(1, 2)), "scalar")
  
  # Test invalid training set
  expect_error(validate_matrix_params(P, train = c("missing")), "not in P")
  expect_error(validate_matrix_params(P, train = c("ind1")), "at least 2")
  expect_error(validate_matrix_params(P, train = c("ind1", "ind1")), "duplicate")
  
  # Test overlap between train and test
  expect_error(validate_matrix_params(P, train = train, test = c("ind1", "ind6")), "overlap")
})

test_that("validate_kinship_matrix works correctly", {
  # Create valid kinship matrix
  K <- matrix(c(1.0, 0.5, 0.5, 1.0), nrow = 2)
  rownames(K) <- colnames(K) <- c("ind1", "ind2")
  
  # Test valid matrix
  expect_silent(validate_kinship_matrix(K, c("ind1", "ind2")))
  
  # Test non-square matrix
  K_bad <- K[1:2, 1, drop = FALSE]  # Keep as matrix
  expect_error(validate_kinship_matrix(K_bad, c("ind1")), "square")
  
  # Test asymmetric matrix
  K_asym <- K
  K_asym[1, 2] <- 0.3
  expect_error(validate_kinship_matrix(K_asym, c("ind1", "ind2")), "symmetric")
  
  # Test negative definite matrix
  K_neg <- matrix(c(1, 2, 2, 1), nrow = 2)
  rownames(K_neg) <- colnames(K_neg) <- c("ind1", "ind2")
  expect_error(validate_kinship_matrix(K_neg, c("ind1", "ind2")), "positive semi-definite")
  
  # Test missing individuals
  expect_error(validate_kinship_matrix(K, c("ind1", "ind3")), "missing individuals")
})

test_that("check_memory_usage works correctly", {
  # Should not warn for small matrices
  expect_silent(check_memory_usage(100, 50, "matrix_mult"))
  
  # Should warn for large matrices (reduce threshold for testing)
  expect_warning(check_memory_usage(5000, 5000, "matrix_mult"), "Large computation")
})

test_that("validate_variance_components works correctly", {
  # Valid components
  expect_silent(validate_variance_components(0.8, 0.2))
  expect_silent(validate_variance_components(1, 0))
  expect_silent(validate_variance_components(0, 1))
  
  # Invalid components
  expect_error(validate_variance_components(-1, 1), "non-negative")
  expect_error(validate_variance_components(1, -1), "non-negative")
  expect_error(validate_variance_components(0, 0), "At least one")
  
  # Cross-validation with heritability
  expect_silent(validate_variance_components(0.8, 0.2, h2 = 0.8))
  expect_error(validate_variance_components(0.8, 0.2, h2 = 0.5), "Inconsistent")
})

test_that("get_adaptive_ridge works correctly", {
  # Well-conditioned matrix
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  X <- scale(X)  # Improve conditioning
  lambda <- get_adaptive_ridge(X)
  expect_true(lambda >= 1e-10)
  expect_true(lambda <= 1e-2)
  
  # Ill-conditioned matrix
  X_bad <- cbind(X, X[, 1] + 1e-10 * rnorm(20))  # Near collinearity
  lambda_bad <- get_adaptive_ridge(X_bad)
  expect_true(lambda_bad > lambda)
})

test_that("check_matrix_dimensions works correctly", {
  A <- matrix(1:6, nrow = 2, ncol = 3)
  B <- matrix(1:6, nrow = 3, ncol = 2)
  C <- matrix(1:4, nrow = 2, ncol = 2)
  
  # Compatible matrices for multiplication
  expect_silent(check_matrix_dimensions(A = A, B = B))  # Can multiply A %*% B
  
  # Same dimensions for addition/subtraction
  D <- matrix(1:6, nrow = 2, ncol = 3)  # Same dimensions as A
  expect_silent(check_matrix_dimensions(A = A, D = D))
  
  # Single matrix
  expect_silent(check_matrix_dimensions(A = A))
  
  # No error for compatible dimensions
  expect_silent(check_matrix_dimensions(matrix(1:4, 2, 2), matrix(1:4, 2, 2)))
})