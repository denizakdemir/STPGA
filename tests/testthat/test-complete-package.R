# Comprehensive test for all completed implementations
library(STPGA)

test_that("Complete GA pipeline works end-to-end", {
  # Create test data
  set.seed(123)
  n_ind <- 50
  n_markers <- 20
  n_select <- 10
  
  # Prediction matrix
  P <- matrix(rnorm(n_ind * n_markers), nrow = n_ind, ncol = n_markers)
  rownames(P) <- paste0("ind", 1:n_ind)
  colnames(P) <- paste0("marker", 1:n_markers)
  
  # Candidates and test sets
  candidates <- rownames(P)[1:40]
  test_set <- rownames(P)[41:50]
  
  # Test single-objective GA (suppress messages from adaptive ridge)
  expect_message({
    ga_result <- subset_ga(
      P = P, 
      Candidates = candidates, 
      Test = test_set, 
      ntoselect = n_select,
      npop = 20,
      niterations = 10,
      criterion = "pev_mean",
      verbose = FALSE
    )
  }, "adaptive ridge|Using adaptive")
  
  # Check results structure
  expect_true(is.list(ga_result))
  expect_true("best_solution" %in% names(ga_result))
  expect_true("best_fitness" %in% names(ga_result))
  expect_true("fitness_history" %in% names(ga_result))
  
  # Check solution validity
  expect_equal(length(ga_result$best_solution), n_select)
  expect_true(all(ga_result$best_solution %in% candidates))
  expect_true(is.numeric(ga_result$best_fitness))
  expect_true(is.finite(ga_result$best_fitness))
})

test_that("Multi-objective optimization works", {
  # Create test data
  set.seed(456)
  n_ind <- 30
  n_markers <- 15
  
  P <- matrix(rnorm(n_ind * n_markers), nrow = n_ind, ncol = n_markers)
  rownames(P) <- paste0("ind", 1:n_ind)
  
  candidates <- rownames(P)[1:25]
  test_set <- rownames(P)[26:30]
  
  # Test multi-objective GA (may have messages)
  suppressMessages({
    mo_result <- subset_ga_multiobjective(
      Pcs = P,
      candidates = candidates,
      test = test_set,
      ntoselect = 8,
      criteria = c("pev_mean", "cd_mean"),
      criteria_types = c(TRUE, FALSE),  # minimize PEV, maximize CD
      plot_directions = c(1, -1),
      npop = 15,
      niterations = 5,
      plot_iterations = FALSE
    )
  })
  
  # Check results
  expect_true(is.list(mo_result))
  expect_true("pareto_solutions" %in% names(mo_result))
  expect_true("pareto_fitness" %in% names(mo_result))
  expect_true(length(mo_result$pareto_solutions) > 0)
})

test_that("Distance functions work correctly", {
  # Test distance calculations
  set.seed(789)
  
  # Create distance matrix
  n <- 20
  dist_matrix <- matrix(runif(n*n, 0, 10), n, n)
  dist_matrix <- (dist_matrix + t(dist_matrix)) / 2  # Make symmetric
  diag(dist_matrix) <- 0
  rownames(dist_matrix) <- colnames(dist_matrix) <- paste0("ind", 1:n)
  
  train <- paste0("ind", 1:10)
  test <- paste0("ind", 11:15)
  
  # Test distance functions
  max_dist <- distance_train_to_test_max(train, test, dist_matrix)
  mean_dist <- distance_train_to_test_mean(train, test, dist_matrix)
  min_internal <- distance_internal_min(train, NULL, dist_matrix)
  mean_internal <- distance_internal_mean(train, NULL, dist_matrix)
  
  expect_true(is.numeric(max_dist) && is.finite(max_dist))
  expect_true(is.numeric(mean_dist) && is.finite(mean_dist))
  expect_true(is.numeric(min_internal) && is.finite(min_internal))
  expect_true(is.numeric(mean_internal) && is.finite(mean_internal))
  
  expect_true(max_dist >= mean_dist)  # Max should be >= mean
})

test_that("Matrix operations work correctly", {
  # Test genomic relationship matrix computation
  set.seed(101)
  n_ind <- 20
  n_markers <- 50
  
  # Create marker matrix with -1, 0, 1 coding
  M <- matrix(sample(c(-1, 0, 1), n_ind * n_markers, replace = TRUE), 
              nrow = n_ind, ncol = n_markers)
  rownames(M) <- paste0("ind", 1:n_ind)
  
  # Test A matrix computation
  expect_silent({
    A_vanraden <- genomic_relationship_matrix(M, method = "vanraden")
    A_pieces <- genomic_relationship_matrix(M, method = "pieces", pieces = 5)
    A_standard <- genomic_relationship_matrix(M, method = "standard")
  })
  
  # Check properties
  expect_equal(dim(A_vanraden), c(n_ind, n_ind))
  expect_true(isSymmetric(A_vanraden))
  expect_equal(rownames(A_vanraden), rownames(M))
  
  # All methods should give similar results
  expect_true(all(abs(A_vanraden - A_pieces) < 1e-10))
})

test_that("Evaluation functions work with different criteria", {
  # Test population evaluation
  set.seed(202)
  n_ind <- 25
  n_markers <- 15
  
  P <- matrix(rnorm(n_ind * n_markers), nrow = n_ind, ncol = n_markers)
  rownames(P) <- paste0("ind", 1:n_ind)
  
  # Create population of solutions
  population <- list(
    paste0("ind", 1:8),
    paste0("ind", 3:10),
    paste0("ind", 5:12),
    paste0("ind", 7:14)
  )
  
  test_set <- paste0("ind", 20:25)
  
  # Test different criteria
  criteria_to_test <- c("a_optimality", "d_optimality", "pev_mean", "cd_mean")
  
  for (criterion in criteria_to_test) {
    suppressMessages({
      fitness <- evaluate_population(population, P, test_set, criterion)
    })
    
    expect_equal(length(fitness), length(population))
    expect_true(all(is.finite(fitness)))
  }
})

test_that("GA operators work correctly", {
  # Test crossover
  parent1 <- c("ind1", "ind2", "ind3", "ind4", "ind5")
  parent2 <- c("ind3", "ind4", "ind6", "ind7", "ind8")
  candidates <- paste0("ind", 1:20)
  
  offspring <- crossover(parent1, parent2, 5, candidates, "uniform")
  expect_equal(length(offspring), 5)
  expect_true(all(offspring %in% candidates))
  
  # Test mutation
  solution <- paste0("ind", 1:5)
  mutated <- mutate_solution(solution, candidates, 2, 5)
  expect_equal(length(mutated), 5)
  expect_true(all(mutated %in% candidates))
  
  # Test population diversity
  pop <- list(
    paste0("ind", 1:5),
    paste0("ind", 2:6),
    paste0("ind", 10:14)
  )
  diversity <- calculate_diversity(pop)
  expect_true(is.numeric(diversity))
  expect_true(diversity >= 0 && diversity <= 1)
})

test_that("Validation functions prevent errors", {
  # Test parameter validation
  P <- matrix(rnorm(100), nrow = 20, ncol = 5)
  rownames(P) <- paste0("ind", 1:20)
  
  # Valid case should work
  expect_silent(validate_matrix_params(P, train = paste0("ind", 1:5), test = paste0("ind", 6:8)))
  
  # Invalid cases should fail
  expect_error(validate_matrix_params("not_a_matrix"), "matrix")
  expect_error(validate_matrix_params(P, train = "missing_ind"), "not in P")
  expect_error(validate_matrix_params(P, lambda = -1), "non-negative")
})

test_that("Performance functions work", {
  # Test progress bar creation
  pb <- create_progress_bar(10, "Test", verbose = FALSE)
  expect_true(is.null(pb) || is.list(pb))
  
  # Test chunked operations
  set.seed(303)
  X <- matrix(rnorm(200), nrow = 20, ncol = 10)
  
  # Test chunked crossprod
  result1 <- crossprod(X)
  result2 <- chunked_crossprod(X, chunk_size = 5)
  
  expect_true(all(abs(result1 - result2) < 1e-10))
})

test_that("Naming standards work", {
  # Test name standardization
  expect_equal(standardize_name("camelCase"), "camel_case")
  expect_equal(standardize_name("PascalCase"), "pascal_case")
  expect_equal(standardize_name("already_snake"), "already_snake")
  
  # Test parameter checking
  param_check <- check_parameter_naming(c("validName", "invalid-name", "TRUE"))
  expect_true(is.list(param_check))
  expect_true(length(param_check$warnings) > 0)
})