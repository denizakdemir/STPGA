# Test criteria functions with validation
library(STPGA)

test_that("Classical optimality criteria work correctly", {
  # Create test data
  set.seed(789)
  P <- matrix(rnorm(100), nrow = 20, ncol = 5)
  rownames(P) <- paste0("ind", 1:20)
  colnames(P) <- paste0("var", 1:5)
  
  train <- paste0("ind", 1:10)
  test <- paste0("ind", 11:15)
  
  # Test A-optimality
  aopt <- a_optimality(train, test, P)
  expect_true(is.numeric(aopt))
  expect_true(length(aopt) == 1)
  expect_true(aopt > 0)
  
  # Test D-optimality
  dopt <- d_optimality(train, test, P)
  expect_true(is.numeric(dopt))
  expect_true(length(dopt) == 1)
  
  # Test E-optimality
  eopt <- e_optimality(train, test, P)
  expect_true(is.numeric(eopt))
  expect_true(length(eopt) == 1)
  
  # Test G-optimality
  gopt <- g_optimality(train, test, P)
  expect_true(is.numeric(gopt))
  expect_true(length(gopt) == 1)
  expect_true(gopt > 0)
  
  # Test I-optimality
  iopt <- i_optimality(train, test, P)
  expect_true(is.numeric(iopt))
  expect_true(length(iopt) == 1)
  expect_true(iopt > 0)
})

test_that("PEV criteria work correctly", {
  set.seed(101)
  P <- matrix(rnorm(80), nrow = 16, ncol = 5)
  rownames(P) <- paste0("ind", 1:16)
  
  train <- paste0("ind", 1:10)
  test <- paste0("ind", 11:14)
  
  # Test mean PEV
  pev_mean_val <- pev_mean(train, test, P)
  expect_true(is.numeric(pev_mean_val))
  expect_true(pev_mean_val > 0)
  
  # Test normalized PEV
  pev_norm <- pev_mean(train, test, P, normalized = TRUE)
  expect_true(pev_norm >= 0)
  expect_true(pev_norm <= 1)
  
  # Test max PEV
  pev_max_val <- pev_max(train, test, P)
  expect_true(pev_max_val >= pev_mean_val)  # Max should be >= mean
  
  # Test coefficient of determination
  cd_val <- cd_mean(train, test, P)
  expect_true(is.numeric(cd_val))
  expect_true(cd_val >= 0)
})

test_that("Mixed model criteria work correctly", {
  set.seed(202)
  n <- 12
  P <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  rownames(P) <- paste0("ind", 1:n)
  
  # Create kinship matrix
  K <- matrix(0.1, nrow = n, ncol = n)
  diag(K) <- 1
  rownames(K) <- colnames(K) <- rownames(P)
  
  train <- paste0("ind", 1:8)
  test <- paste0("ind", 9:10)
  
  # Test mixed model PEV
  pev_mm <- pev_mean_mm(train, test, P, K, Vg = 0.8, Ve = 0.2)
  expect_true(is.numeric(pev_mm))
  expect_true(pev_mm >= 0)  # Should be non-negative after fix
  
  # Test mixed model CD
  cd_mm <- cd_mean_mm(train, test, P, K, Vg = 0.8, Ve = 0.2)
  expect_true(is.numeric(cd_mm))
  expect_true(cd_mm >= 0)
  expect_true(cd_mm <= 1)
  
  # Test heritability-based wrappers
  pev_h2 <- pev_mean_mm_h2(train, test, P, K, h2 = 0.8)
  expect_true(is.numeric(pev_h2))
  
  cd_h2 <- cd_mean_mm_h2(train, test, P, K, h2 = 0.8)
  expect_true(is.numeric(cd_h2))
  expect_true(cd_h2 >= 0)
  expect_true(cd_h2 <= 1)
})

test_that("Criteria input validation works", {
  P <- matrix(rnorm(50), nrow = 10, ncol = 5)
  rownames(P) <- paste0("ind", 1:10)
  
  train <- c("ind1", "ind2", "ind3")
  test <- c("ind4", "ind5")
  
  # Valid call should work
  expect_silent(a_optimality(train, test, P))
  
  # Invalid training set should fail
  expect_error(a_optimality(c("missing"), test, P), "not in P")
  
  # Overlapping train/test should fail
  expect_error(a_optimality(c("ind1", "ind2"), c("ind2", "ind3"), P), "overlap")
  
  # Invalid lambda should fail
  expect_error(a_optimality(train, test, P, lambda = -1), "non-negative")
})

test_that("Unified criterion function works", {
  set.seed(303)
  P <- matrix(rnorm(60), nrow = 12, ncol = 5)
  rownames(P) <- paste0("ind", 1:12)
  
  train <- paste0("ind", 1:8)
  test <- paste0("ind", 9:10)
  
  # Test modern criterion names
  aopt_unified <- criterion(train, test, P, criterion = "a_optimality")
  aopt_direct <- a_optimality(train, test, P)
  expect_equal(aopt_unified, aopt_direct)
  
  # Test legacy names
  aopt_legacy <- criterion(train, test, P, criterion = "AOPT")
  expect_equal(aopt_legacy, aopt_direct)
  
  # Test normalized versions
  pev_norm <- criterion(train, test, P, criterion = "pev_mean_normalized")
  expect_true(pev_norm >= 0)
  expect_true(pev_norm <= 1)
  
  # Test unknown criterion
  expect_error(criterion(train, test, P, criterion = "unknown"), "Unknown criterion")
})

test_that("Legacy compatibility functions work", {
  set.seed(404)
  P <- matrix(rnorm(40), nrow = 8, ncol = 5)
  rownames(P) <- paste0("ind", 1:8)
  
  train <- paste0("ind", 1:5)
  test <- paste0("ind", 6:7)
  
  # Test legacy influence measure
  influence <- influence_measure_legacy(train, test, P)
  expect_true(is.numeric(influence))
  expect_true(influence > 0)
  
  # Should be different from proper CD
  cd_proper <- cd_mean(train, test, P)
  expect_true(abs(influence - cd_proper) > 1e-10)  # Should be different
})