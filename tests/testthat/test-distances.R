test_that("distance_to_ideal function works correctly", {
  # Create test data
  set.seed(123)
  X <- matrix(c(
    1, 2, 3,
    4, 5, 6,
    0, 1, 2,
    3, 4, 5
  ), nrow = 4, byrow = TRUE)
  
  # Test euclidean distance
  dist_euclidean <- distance_to_ideal(X, method = "euclidean")
  expect_type(dist_euclidean, "double")
  expect_length(dist_euclidean, 4)
  expect_true(all(dist_euclidean >= 0))
  
  # Test manhattan distance
  dist_manhattan <- distance_to_ideal(X, method = "manhattan")
  expect_type(dist_manhattan, "double")
  expect_length(dist_manhattan, 4)
  expect_true(all(dist_manhattan >= 0))
  
  # Test chebyshev distance
  dist_chebyshev <- distance_to_ideal(X, method = "chebyshev")
  expect_type(dist_chebyshev, "double")
  expect_length(dist_chebyshev, 4)
  expect_true(all(dist_chebyshev >= 0))
  
  # Test with zero ranges (constant columns)
  X_constant <- cbind(X, rep(5, 4))
  expect_warning(distance_to_ideal(X_constant, handle_zeros = "warning"))
  
  # Test error handling for zero ranges
  expect_error(distance_to_ideal(X_constant, handle_zeros = "error"))
  
  # Test ignore option for zero ranges
  dist_ignore <- distance_to_ideal(X_constant, handle_zeros = "ignore")
  expect_type(dist_ignore, "double")
  
  # Test empty matrix
  expect_equal(distance_to_ideal(matrix(numeric(0), 0, 0)), numeric(0))
  
  # Test data.frame input
  X_df <- as.data.frame(X)
  dist_df <- distance_to_ideal(X_df)
  expect_equal(dist_df, dist_euclidean)
  
  # Test invalid method
  expect_error(distance_to_ideal(X, method = "invalid"))
})

test_that("distance matrix functions work correctly", {
  # Create test distance matrix
  set.seed(123)
  n <- 20
  individuals <- paste0("ind", 1:n)
  dist_matrix <- matrix(runif(n * n, 0, 10), n, n)
  diag(dist_matrix) <- 0
  dist_matrix[lower.tri(dist_matrix)] <- t(dist_matrix)[lower.tri(dist_matrix)]
  rownames(dist_matrix) <- colnames(dist_matrix) <- individuals
  
  train <- individuals[1:10]
  test <- individuals[11:15]
  
  # Test max distance from training to test
  max_dist <- distance_train_to_test_max(train, test, dist_matrix)
  expect_type(max_dist, "double")
  expect_length(max_dist, 1)
  expect_true(max_dist >= 0)
  
  # Test mean distance from training to test
  mean_dist <- distance_train_to_test_mean(train, test, dist_matrix)
  expect_type(mean_dist, "double")
  expect_length(mean_dist, 1)
  expect_true(mean_dist >= 0)
  expect_true(mean_dist <= max_dist)  # Mean should be <= max
  
  # Test internal distances
  neg_min_internal <- distance_internal_min(train, distance_matrix = dist_matrix)
  neg_mean_internal <- distance_internal_mean(train, distance_matrix = dist_matrix)
  
  expect_type(neg_min_internal, "double")
  expect_type(neg_mean_internal, "double")
  expect_true(neg_min_internal <= 0)  # Should be negative
  expect_true(neg_mean_internal <= 0)  # Should be negative
  
  # Test with empty sets
  expect_equal(distance_train_to_test_max(character(0), test, dist_matrix), Inf)
  expect_equal(distance_train_to_test_mean(train, character(0), dist_matrix), Inf)
  expect_equal(distance_internal_min(character(0)), 0)
  expect_equal(distance_internal_min(train[1]), 0)  # Single individual
  
  # Test warning for missing individuals
  invalid_train <- c(train[1:5], "invalid_ind")
  expect_warning(distance_train_to_test_max(invalid_train, test, dist_matrix))
})

test_that("distance_criterion unified function works", {
  # Create test data
  set.seed(123)
  n <- 15
  individuals <- paste0("ind", 1:n)
  dist_matrix <- matrix(runif(n * n, 0, 5), n, n)
  diag(dist_matrix) <- 0
  dist_matrix[lower.tri(dist_matrix)] <- t(dist_matrix)[lower.tri(dist_matrix)]
  rownames(dist_matrix) <- colnames(dist_matrix) <- individuals
  
  train <- individuals[1:8]
  test <- individuals[9:12]
  
  # Test all criterion types
  max_to_test <- distance_criterion(train, test, dist_matrix, "max_to_test")
  mean_to_test <- distance_criterion(train, test, dist_matrix, "mean_to_test")
  neg_min_internal <- distance_criterion(train, test, dist_matrix, "neg_min_internal")
  neg_mean_internal <- distance_criterion(train, test, dist_matrix, "neg_mean_internal")
  
  expect_type(max_to_test, "double")
  expect_type(mean_to_test, "double")
  expect_type(neg_min_internal, "double")
  expect_type(neg_mean_internal, "double")
  
  # Test that they match individual functions
  expect_equal(max_to_test, distance_train_to_test_max(train, test, dist_matrix))
  expect_equal(mean_to_test, distance_train_to_test_mean(train, test, dist_matrix))
  expect_equal(neg_min_internal, distance_internal_min(train, test, dist_matrix))
  expect_equal(neg_mean_internal, distance_internal_mean(train, test, dist_matrix))
  
  # Test invalid criterion
  expect_error(distance_criterion(train, test, dist_matrix, "invalid"))
})

test_that("population_distances function works correctly", {
  set.seed(123)
  candidates <- paste0("item", 1:30)
  
  # Create test population
  population <- list()
  for (i in 1:6) {
    population[[i]] <- sample(candidates, 10)
  }
  
  # Test hamming distance
  dist_hamming <- population_distances(population, "hamming")
  expect_type(dist_hamming, "double")
  expect_equal(dim(dist_hamming), c(6, 6))
  expect_true(all(diag(dist_hamming) == 0))  # Diagonal should be 0
  expect_true(all(dist_hamming >= 0 & dist_hamming <= 1))  # Should be in [0,1]
  
  # Test jaccard distance
  dist_jaccard <- population_distances(population, "jaccard")
  expect_equal(dim(dist_jaccard), c(6, 6))
  expect_true(all(diag(dist_jaccard) == 0))
  expect_true(all(dist_jaccard >= 0 & dist_jaccard <= 1))
  
  # Test sorensen distance
  dist_sorensen <- population_distances(population, "sorensen")
  expect_equal(dim(dist_sorensen), c(6, 6))
  expect_true(all(diag(dist_sorensen) == 0))
  expect_true(all(dist_sorensen >= 0 & dist_sorensen <= 1))
  
  # Test symmetry
  expect_true(all(abs(dist_hamming - t(dist_hamming)) < 1e-10))
  expect_true(all(abs(dist_jaccard - t(dist_jaccard)) < 1e-10))
  expect_true(all(abs(dist_sorensen - t(dist_sorensen)) < 1e-10))
  
  # Test with identical solutions
  identical_pop <- list(population[[1]], population[[1]])
  dist_identical <- population_distances(identical_pop, "jaccard")
  expect_equal(dist_identical[1, 2], 0)
  expect_equal(dist_identical[2, 1], 0)
  
  # Test with single solution
  single_pop <- list(population[[1]])
  dist_single <- population_distances(single_pop, "jaccard")
  expect_equal(dim(dist_single), c(1, 1))
  expect_equal(dist_single[1, 1], 0)
  
  # Test empty population
  empty_pop <- list()
  dist_empty <- population_distances(empty_pop, "jaccard")
  expect_equal(dim(dist_empty), c(0, 0))
})

test_that("diversity_summary function works correctly", {
  set.seed(123)
  candidates <- paste0("item", 1:25)
  
  # Create diverse population
  population_diverse <- list()
  for (i in 1:8) {
    population_diverse[[i]] <- sample(candidates, 8)
  }
  
  # Test diversity summary
  summary_diverse <- diversity_summary(population_diverse, "jaccard")
  
  expect_type(summary_diverse, "list")
  expect_true("mean_distance" %in% names(summary_diverse))
  expect_true("min_distance" %in% names(summary_diverse))
  expect_true("max_distance" %in% names(summary_diverse))
  expect_true("diversity_index" %in% names(summary_diverse))
  expect_true("distance_matrix" %in% names(summary_diverse))
  
  expect_true(summary_diverse$mean_distance >= 0)
  expect_true(summary_diverse$min_distance >= 0)
  expect_true(summary_diverse$max_distance >= summary_diverse$mean_distance)
  expect_true(summary_diverse$mean_distance >= summary_diverse$min_distance)
  expect_equal(summary_diverse$diversity_index, summary_diverse$mean_distance)
  
  # Test with identical population (low diversity)
  base_solution <- sample(candidates, 8)
  population_identical <- rep(list(base_solution), 5)
  summary_identical <- diversity_summary(population_identical)
  
  expect_equal(summary_identical$mean_distance, 0)
  expect_equal(summary_identical$min_distance, 0)
  expect_equal(summary_identical$max_distance, 0)
  expect_equal(summary_identical$diversity_index, 0)
  
  # Test with single solution
  population_single <- list(sample(candidates, 8))
  summary_single <- diversity_summary(population_single)
  
  expect_equal(summary_single$mean_distance, 0)
  expect_equal(summary_single$diversity_index, 0)
  
  # Test empty population
  population_empty <- list()
  summary_empty <- diversity_summary(population_empty)
  
  expect_equal(summary_empty$mean_distance, 0)
  expect_equal(summary_empty$diversity_index, 0)
})

test_that("legacy wrapper functions work correctly", {
  # Test that legacy functions match new implementations
  set.seed(123)
  X <- matrix(rnorm(20), 5, 4)
  
  # Test disttoideal wrapper
  dist_new <- distance_to_ideal(X)
  dist_legacy <- disttoideal(X)
  expect_equal(dist_new, dist_legacy)
  
  # Test distance matrix wrappers
  n <- 10
  individuals <- paste0("ind", 1:n)
  dist_matrix <- matrix(runif(n * n, 0, 5), n, n)
  diag(dist_matrix) <- 0
  rownames(dist_matrix) <- colnames(dist_matrix) <- individuals
  
  train <- individuals[1:5]
  test <- individuals[6:8]
  
  # Test dist_to_test wrappers
  expect_equal(dist_to_test(train, test, dist_matrix),
               distance_train_to_test_max(train, test, dist_matrix))
  expect_equal(dist_to_test2(train, test, dist_matrix),
               distance_train_to_test_mean(train, test, dist_matrix))
  
  # Test internal distance wrappers
  expect_equal(neg_dist_in_train(train, test, dist_matrix),
               distance_internal_min(train, test, dist_matrix))
  expect_equal(neg_dist_in_train2(train, test, dist_matrix),
               distance_internal_mean(train, test, dist_matrix))
  
  # Test unified wrapper
  expect_equal(unified_distance_criterion(train, test, dist_matrix, "max_to_test"),
               distance_criterion(train, test, dist_matrix, "max_to_test"))
  
  # Test population distances wrapper
  population <- list(
    sample(paste0("item", 1:20), 8),
    sample(paste0("item", 1:20), 8)
  )
  expect_equal(calculate_population_distances(population, "jaccard"),
               population_distances(population, "jaccard"))
})

test_that("input validation works correctly", {
  # Test invalid matrix input for distance_to_ideal
  expect_error(distance_to_ideal("not_a_matrix"))
  expect_error(distance_to_ideal(list(1, 2, 3)))
  
  # Test invalid method
  X <- matrix(rnorm(12), 3, 4)
  expect_error(distance_to_ideal(X, method = "invalid_method"))
  expect_error(population_distances(list(), "invalid_method"))
  
  # Test empty inputs for distance matrix functions
  empty_matrix <- matrix(numeric(0), 0, 0)
  expect_error(distance_train_to_test_max(character(0), "test", empty_matrix))
  
  # Test mismatched names in distance matrix
  dist_matrix <- matrix(1:9, 3, 3)
  rownames(dist_matrix) <- c("a", "b", "c")
  colnames(dist_matrix) <- c("x", "y", "z")
  expect_warning(distance_train_to_test_max(c("a"), c("missing"), dist_matrix))
})

test_that("numerical properties are maintained", {
  set.seed(123)
  
  # Test that distance functions return reasonable values
  X <- matrix(runif(50, 0, 10), 10, 5)
  distances <- distance_to_ideal(X)
  
  expect_true(all(is.finite(distances)))
  expect_true(all(distances >= 0))
  expect_true(all(distances <= sqrt(5)))  # Max possible distance for normalized [0,1]^5
  
  # Test triangle inequality for population distances
  candidates <- paste0("item", 1:15)
  pop1 <- sample(candidates, 6)
  pop2 <- sample(candidates, 6)
  pop3 <- sample(candidates, 6)
  
  population <- list(pop1, pop2, pop3)
  dist_matrix <- population_distances(population, "jaccard")
  
  # Triangle inequality: d(i,k) <= d(i,j) + d(j,k)
  for (i in 1:3) {
    for (j in 1:3) {
      for (k in 1:3) {
        expect_true(dist_matrix[i, k] <= dist_matrix[i, j] + dist_matrix[j, k] + 1e-10)
      }
    }
  }
})
