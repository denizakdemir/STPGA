test_that("evaluate_population function works correctly", {
  # Create test data
  set.seed(123)
  n <- 30
  p <- 8
  P <- matrix(rnorm(n * p), n, p)
  rownames(P) <- paste0("ind", 1:n)
  
  # Create test population
  population <- list()
  for (i in 1:5) {
    population[[i]] <- sample(rownames(P), 10)
  }
  
  test_set <- paste0("ind", 21:25)
  
  # Test basic functionality
  fitness <- evaluate_population(population, P, test_set, "pev_mean_normalized")
  
  expect_type(fitness, "double")
  expect_length(fitness, 5)
  expect_true(all(is.finite(fitness)))
  
  # Test different criteria
  fitness_aopt <- evaluate_population(population, P, test_set, "a_optimality")
  fitness_dopt <- evaluate_population(population, P, test_set, "d_optimality")
  fitness_cd <- evaluate_population(population, P, test_set, "cd_mean")
  
  expect_type(fitness_aopt, "double")
  expect_type(fitness_dopt, "double")
  expect_type(fitness_cd, "double")
  expect_length(fitness_aopt, 5)
  expect_length(fitness_dopt, 5)
  expect_length(fitness_cd, 5)
  
  # Test without test set
  fitness_single <- evaluate_population(population, P, NULL, "a_optimality")
  expect_type(fitness_single, "double")
  expect_length(fitness_single, 5)
  
  # Test with caching
  fitness_cached1 <- evaluate_population(population, P, test_set, "pev_mean_normalized", use_cache = TRUE)
  fitness_cached2 <- evaluate_population(population, P, test_set, "pev_mean_normalized", use_cache = TRUE)
  expect_equal(fitness_cached1, fitness_cached2)
  
  # Test without caching
  fitness_no_cache <- evaluate_population(population, P, test_set, "pev_mean_normalized", use_cache = FALSE)
  expect_equal(fitness, fitness_no_cache)
  
  # Test with mixed models
  K <- matrix(0.1, n, n)
  diag(K) <- 1
  rownames(K) <- colnames(K) <- rownames(P)
  
  fitness_mm <- evaluate_population(population, P, test_set, "pev_mean_mm", K = K)
  expect_type(fitness_mm, "double")
  expect_length(fitness_mm, 5)
  expect_true(all(is.finite(fitness_mm)))
})

test_that("evaluate_population_smart function works correctly", {
  # Create test data
  set.seed(123)
  n <- 40
  p <- 6
  P <- matrix(rnorm(n * p), n, p)
  rownames(P) <- paste0("ind", 1:n)
  
  # Create larger population for batch testing
  population <- list()
  for (i in 1:15) {
    population[[i]] <- sample(rownames(P), 12)
  }
  
  test_set <- paste0("ind", 31:35)
  
  # Test basic functionality
  fitness_smart <- evaluate_population_smart(population, P, test_set, "pev_mean")
  
  expect_type(fitness_smart, "double")
  expect_length(fitness_smart, 15)
  expect_true(all(is.finite(fitness_smart)))
  
  # Compare with regular evaluation
  fitness_regular <- evaluate_population(population, P, test_set, "pev_mean", use_cache = FALSE)
  
  # Results should be very similar
  correlation <- cor(fitness_smart, fitness_regular)
  expect_true(correlation > 0.99)
  
  # Test with custom batch size
  fitness_batch <- evaluate_population_smart(population, P, test_set, "pev_mean", 
                                             batch_size = 5)
  expect_length(fitness_batch, 15)
  
  # Test with different criteria
  fitness_smart_cd <- evaluate_population_smart(population, P, test_set, "cd_mean")
  expect_type(fitness_smart_cd, "double")
  expect_length(fitness_smart_cd, 15)
})

test_that("create_evaluation_cache works correctly", {
  # Test cache creation and basic operations
  cache <- create_evaluation_cache()
  
  expect_type(cache, "list")
  expect_true("get" %in% names(cache))
  expect_true("store" %in% names(cache))
  expect_true("clear" %in% names(cache))
  expect_true("size" %in% names(cache))
  
  # Test storing and retrieving
  solution <- c("ind1", "ind2", "ind3")
  test_set <- c("ind4", "ind5")
  criterion <- "pev_mean"
  lambda <- 1e-6
  fitness_value <- 0.5
  
  # Initially should be empty
  expect_null(cache$get(solution, test_set, criterion, lambda))
  expect_equal(cache$size(), 0)
  
  # Store value
  cache$store(solution, test_set, criterion, lambda, fitness_value)
  expect_equal(cache$size(), 1)
  
  # Retrieve value
  retrieved <- cache$get(solution, test_set, criterion, lambda)
  expect_equal(retrieved, fitness_value)
  
  # Store another value
  solution2 <- c("ind6", "ind7", "ind8")
  cache$store(solution2, test_set, criterion, lambda, 0.8)
  expect_equal(cache$size(), 2)
  
  # Clear cache
  cache$clear()
  expect_equal(cache$size(), 0)
  expect_null(cache$get(solution, test_set, criterion, lambda))
  
  # Test with NULL test set
  cache$store(solution, NULL, criterion, lambda, 0.3)
  retrieved_null <- cache$get(solution, NULL, criterion, lambda)
  expect_equal(retrieved_null, 0.3)
})

test_that("transform_fitness function works correctly", {
  set.seed(123)
  
  # Test with mixed positive/negative fitness values
  fitness <- c(-2, -1, 0, 1, 2, 3)
  
  # Test linear transformation
  fitness_linear <- transform_fitness(fitness, "linear")
  expect_type(fitness_linear, "double")
  expect_length(fitness_linear, 6)
  expect_true(all(fitness_linear > 0))  # Should be all positive
  expect_equal(min(fitness_linear), 1)  # Minimum should be 1
  
  # Test rank transformation
  fitness_rank <- transform_fitness(fitness, "rank")
  expect_type(fitness_rank, "double")
  expect_length(fitness_rank, 6)
  expect_equal(sort(fitness_rank), 1:6)  # Should be ranks 1-6
  
  # Test tournament transformation
  fitness_tournament <- transform_fitness(fitness, "tournament", scaling_factor = 3)
  expect_type(fitness_tournament, "double")
  expect_length(fitness_tournament, 6)
  expect_true(all(fitness_tournament >= 1))
  expect_true(max(fitness_tournament) <= 3)
  
  # Test with all positive values
  fitness_positive <- c(1, 2, 3, 4, 5)
  fitness_pos_linear <- transform_fitness(fitness_positive, "linear")
  expect_equal(fitness_pos_linear, fitness_positive)  # Should be unchanged
  
  # Test with constant values
  fitness_constant <- rep(5, 4)
  fitness_const_tournament <- transform_fitness(fitness_constant, "tournament")
  expect_true(all(fitness_const_tournament == 1))  # Should all be 1
  
  # Test with NA values
  fitness_na <- c(1, 2, NA, 4, 5)
  fitness_na_rank <- transform_fitness(fitness_na, "rank")
  expect_true(is.na(fitness_na_rank[3]))  # NA should remain NA
  
  # Test invalid method
  expect_error(transform_fitness(fitness, "invalid_method"))
})

test_that("compute_population_stats function works correctly", {
  set.seed(123)
  candidates <- paste0("ind", 1:30)
  
  # Create test population
  population <- list()
  for (i in 1:8) {
    population[[i]] <- sample(candidates, 10)
  }
  
  # Create fitness values with some invalid ones
  fitness <- c(0.5, 0.8, 1.2, Inf, 0.9, 1.1, NA, 0.7)
  
  # Test basic functionality
  stats <- compute_population_stats(population, fitness)
  
  expect_type(stats, "list")
  expect_true("best_fitness" %in% names(stats))
  expect_true("worst_fitness" %in% names(stats))
  expect_true("mean_fitness" %in% names(stats))
  expect_true("fitness_std" %in% names(stats))
  expect_true("diversity" %in% names(stats))
  expect_true("convergence_metric" %in% names(stats))
  expect_true("valid_solutions" %in% names(stats))
  expect_true("invalid_solutions" %in% names(stats))
  
  # Check values
  expect_equal(stats$best_fitness, 0.5)  # Minimum of valid values
  expect_equal(stats$worst_fitness, 1.2)  # Maximum of valid values
  expect_equal(stats$valid_solutions, 6)  # Should exclude Inf and NA
  expect_equal(stats$invalid_solutions, 2)  # Inf and NA
  
  expect_true(stats$mean_fitness > 0)
  expect_true(stats$fitness_std >= 0)
  expect_true(stats$diversity >= 0 && stats$diversity <= 1)
  expect_true(stats$convergence_metric >= 0)
  
  # Test with all valid fitness
  fitness_valid <- c(0.5, 0.8, 1.2, 0.9, 1.1, 0.7, 1.0, 0.6)
  stats_valid <- compute_population_stats(population, fitness_valid)
  
  expect_equal(stats_valid$valid_solutions, 8)
  expect_equal(stats_valid$invalid_solutions, 0)
  
  # Test with all invalid fitness
  fitness_invalid <- c(Inf, NA, Inf, NA)
  population_small <- population[1:4]
  stats_invalid <- compute_population_stats(population_small, fitness_invalid)
  
  expect_equal(stats_invalid$best_fitness, Inf)
  expect_equal(stats_invalid$worst_fitness, Inf)
  expect_equal(stats_invalid$mean_fitness, Inf)
  expect_equal(stats_invalid$valid_solutions, 0)
  expect_equal(stats_invalid$invalid_solutions, 4)
  
  # Test convergence metric calculation
  fitness_converged <- rep(1.0, 5)  # Very low variance
  population_converged <- population[1:5]
  stats_converged <- compute_population_stats(population_converged, fitness_converged)
  expect_equal(stats_converged$convergence_metric, 0)  # Should be 0 for no variance
})

test_that("batch processing functions work correctly", {
  # Test determine_optimal_batch_size
  batch_size_sequential <- determine_optimal_batch_size(100, 1)
  batch_size_parallel <- determine_optimal_batch_size(100, 4)
  
  expect_type(batch_size_sequential, "double")
  expect_type(batch_size_parallel, "double")
  expect_true(batch_size_sequential <= 50)
  expect_true(batch_size_parallel >= 10)
  
  # Test group_similar_solutions
  set.seed(123)
  candidates <- paste0("item", 1:50)
  population <- list()
  for (i in 1:12) {
    size <- sample(8:15, 1)  # Variable sizes
    population[[i]] <- sample(candidates, size)
  }
  
  groups <- group_similar_solutions(population, target_batch_size = 4)
  
  expect_type(groups, "list")
  expect_true(length(groups) >= 1)
  
  # Check that all solutions are included
  all_indices <- unlist(lapply(groups, function(g) g$indices))
  expect_equal(sort(all_indices), 1:12)
  
  # Test with small population
  small_pop <- population[1:3]
  groups_small <- group_similar_solutions(small_pop, target_batch_size = 5)
  expect_length(groups_small, 1)
  expect_equal(groups_small[[1]]$indices, 1:3)
})

test_that("legacy wrapper functions work correctly", {
  # Create test data
  set.seed(123)
  n <- 25
  p <- 5
  P <- matrix(rnorm(n * p), n, p)
  rownames(P) <- paste0("ind", 1:n)
  
  population <- list()
  for (i in 1:4) {
    population[[i]] <- sample(rownames(P), 8)
  }
  
  test_set <- paste0("ind", 20:23)
  
  # Test evaluate_population_optimized wrapper
  fitness_new <- evaluate_population(population, P, test_set, "pev_mean")
  fitness_legacy <- evaluate_population_optimized(population, P, test_set, "pev_mean")
  
  expect_equal(fitness_new, fitness_legacy)
  
  # Test evaluate_population_legacy wrapper
  fitness_legacy2 <- evaluate_population_legacy(population, "pev_mean", test_set, P)
  expect_equal(fitness_new, fitness_legacy2)
  
  # Test with additional parameters
  fitness_legacy_params <- evaluate_population_optimized(population, P, test_set, 
                                                        "pev_mean", lambda = 1e-5, 
                                                        mc.cores = 1)
  expect_type(fitness_legacy_params, "double")
  expect_length(fitness_legacy_params, 4)
})

test_that("error handling works correctly", {
  # Create test data
  set.seed(123)
  n <- 20
  p <- 4
  P <- matrix(rnorm(n * p), n, p)
  rownames(P) <- paste0("ind", 1:n)
  
  population <- list(
    sample(rownames(P), 6),
    sample(rownames(P), 6)
  )
  
  test_set <- paste0("ind", 15:18)
  
  # Test with invalid criterion
  expect_warning(fitness_invalid <- evaluate_population(population, P, test_set, "INVALID"))
  expect_true(all(is.infinite(fitness_invalid)))
  
  # Test with missing individuals in population
  invalid_population <- list(c("ind1", "invalid_ind", "ind3"))
  expect_warning(fitness_missing <- evaluate_population(invalid_population, P, test_set, "a_optimality"))
  expect_true(is.infinite(fitness_missing[1]))
  
  # Test transform_fitness with invalid method
  fitness <- c(1, 2, 3)
  expect_error(transform_fitness(fitness, "invalid"))
  
  # Test cache with invalid inputs
  cache <- create_evaluation_cache()
  expect_null(cache$get(character(0), test_set, "a_optimality", 1e-6))
})

test_that("parallel processing works correctly", {
  # Note: Testing with mc.cores = 1 to avoid issues in test environment
  set.seed(123)
  n <- 30
  p <- 6
  P <- matrix(rnorm(n * p), n, p)
  rownames(P) <- paste0("ind", 1:n)
  
  # Create larger population for parallel testing
  population <- list()
  for (i in 1:20) {
    population[[i]] <- sample(rownames(P), 10)
  }
  
  test_set <- paste0("ind", 25:28)
  
  # Test sequential vs parallel (both with mc.cores = 1 for testing)
  fitness_sequential <- evaluate_population(population, P, test_set, "pev_mean", mc.cores = 1)
  fitness_parallel <- evaluate_population(population, P, test_set, "pev_mean", mc.cores = 1)
  
  expect_equal(fitness_sequential, fitness_parallel)
  
  # Test smart evaluation with parallel processing
  fitness_smart_seq <- evaluate_population_smart(population, P, test_set, "pev_mean", mc.cores = 1)
  fitness_smart_par <- evaluate_population_smart(population, P, test_set, "pev_mean", mc.cores = 1)
  
  expect_equal(fitness_smart_seq, fitness_smart_par)
})

test_that("numerical stability and edge cases", {
  # Test with ill-conditioned prediction matrix
  set.seed(123)
  n <- 20
  p <- 18  # Close to n
  P <- matrix(rnorm(n * p), n, p)
  # Add some collinearity
  P[, p] <- P[, 1] + rnorm(n, sd = 0.01)
  rownames(P) <- paste0("ind", 1:n)
  
  population <- list(
    sample(rownames(P), 12),
    sample(rownames(P), 12)
  )
  
  test_set <- paste0("ind", 18:20)
  
  # Should handle ill-conditioned matrices gracefully
  fitness_ill <- evaluate_population(population, P, test_set, "a_optimality", lambda = 1e-3)
  expect_true(all(is.finite(fitness_ill)))
  
  # Test with very small training sets
  small_population <- list(rownames(P)[1:2])
  fitness_small <- evaluate_population(small_population, P, test_set, "pev_mean")
  expect_true(is.finite(fitness_small[1]))
  
  # Test with identical solutions in population
  identical_population <- rep(list(rownames(P)[1:8]), 3)
  fitness_identical <- evaluate_population(identical_population, P, test_set, "a_optimality")
  expect_equal(fitness_identical[1], fitness_identical[2])
  expect_equal(fitness_identical[2], fitness_identical[3])
  
  # Test diversity calculation with identical solutions
  diversity_identical <- calculate_diversity(identical_population)
  expect_equal(diversity_identical, 0)
})