# Test multi-criteria convergence detection
library(STPGA)

test_that("Multi-criteria convergence detection works", {
  # Test the convergence function directly
  set.seed(123)
  window_size <- 10
  tolerance <- 1e-6
  
  # Create mock convergence history with stable patterns
  convergence_history <- list(
    fitness = c(rep(10.0, 10)),  # Very stable fitness
    diversity = c(rep(0.5, 10)),   # Stable diversity  
    improvement_rate = c(rep(1e-8, 10)),  # Very low improvement
    fitness_variance = c(rep(0.001, 10))  # Very stable variance
  )
  
  result <- check_multi_criteria_convergence(
    convergence_history, 10, window_size, tolerance
  )
  
  expect_true(is.list(result))
  expect_true(result$converged)
  expect_true(is.character(result$reason))
  expect_true(length(result$criteria_details) > 0)
})

test_that("Multi-criteria convergence with real GA run", {
  # Create small test data
  set.seed(456)
  n_ind <- 20
  n_markers <- 8
  n_select <- 5
  
  P <- matrix(rnorm(n_ind * n_markers), nrow = n_ind, ncol = n_markers)
  rownames(P) <- paste0("ind", 1:n_ind)
  colnames(P) <- paste0("marker", 1:n_markers)
  
  candidates <- rownames(P)[1:15]
  test_set <- rownames(P)[16:20]
  
  # Run GA with small parameters to ensure quick convergence
  suppressMessages({
    ga_result <- subset_ga(
      P = P,
      Candidates = candidates,
      Test = test_set,
      ntoselect = n_select,
      npop = 10,
      niterations = 20,
      criterion = "pev_mean",
      verbose = FALSE,
      tolconv = 1e-3,  # More relaxed tolerance for quick convergence
      minitbefstop = 5
    )
  })
  
  # Check that convergence history is included
  expect_true("convergence_history" %in% names(ga_result))
  expect_true(is.list(ga_result$convergence_history))
  expect_true("fitness" %in% names(ga_result$convergence_history))
  expect_true("diversity" %in% names(ga_result$convergence_history))
  expect_true("improvement_rate" %in% names(ga_result$convergence_history))
  expect_true("fitness_variance" %in% names(ga_result$convergence_history))
  
  # Check that convergence window parameter is included
  expect_true("convergence_window" %in% names(ga_result$parameters))
  expect_true(ga_result$parameters$convergence_window >= 5)
})

test_that("Convergence history tracking with configurable windows", {
  # Test different window sizes
  set.seed(789)
  n_ind <- 15
  n_markers <- 6
  
  P <- matrix(rnorm(n_ind * n_markers), nrow = n_ind, ncol = n_markers)
  rownames(P) <- paste0("ind", 1:n_ind)
  
  candidates <- rownames(P)[1:12]
  test_set <- rownames(P)[13:15]
  
  # Test with different convergence window multipliers
  suppressMessages({
    ga_small_window <- subset_ga(
      P = P, Candidates = candidates, Test = test_set, ntoselect = 4,
      npop = 8, niterations = 15, convergence_window_multiplier = 2,
      niterreg = 3, verbose = FALSE
    )
    
    ga_large_window <- subset_ga(
      P = P, Candidates = candidates, Test = test_set, ntoselect = 4,
      npop = 8, niterations = 15, convergence_window_multiplier = 6,
      niterreg = 3, verbose = FALSE
    )
  })
  
  expect_equal(ga_small_window$parameters$convergence_window, 10)  # max(3 * 2, 10) = 10
  expect_equal(ga_large_window$parameters$convergence_window, 18) # max(3 * 6, 10) = 18
})

test_that("Convergence diagnostics function works", {
  # Create test GA result
  set.seed(101)
  n_ind <- 12
  n_markers <- 5
  
  P <- matrix(rnorm(n_ind * n_markers), nrow = n_ind, ncol = n_markers)
  rownames(P) <- paste0("ind", 1:n_ind)
  
  candidates <- rownames(P)[1:10]
  test_set <- rownames(P)[11:12]
  
  suppressMessages({
    ga_result <- subset_ga(
      P = P, Candidates = candidates, Test = test_set, ntoselect = 3,
      npop = 6, niterations = 10, verbose = FALSE
    )
  })
  
  # Test diagnostics function
  diags <- convergence_diagnostics(ga_result, plot = FALSE)
  
  expect_true(is.list(diags))
  expect_true("convergence_status" %in% names(diags))
  expect_true("metrics" %in% names(diags))
  expect_true("recommendations" %in% names(diags))
  expect_true(is.numeric(diags$window_size))
  expect_true(diags$data_points > 0)
})

test_that("Restart mechanism works correctly", {
  # Create test data that should converge quickly
  set.seed(202)
  n_ind <- 15
  n_markers <- 5
  
  P <- matrix(rnorm(n_ind * n_markers), nrow = n_ind, ncol = n_markers)
  rownames(P) <- paste0("ind", 1:n_ind)
  
  candidates <- rownames(P)[1:12]
  test_set <- rownames(P)[13:15]
  
  # Test restart mechanism with very low restart threshold to force restarts
  suppressMessages({
    ga_with_restart <- subset_ga(
      P = P, Candidates = candidates, Test = test_set, ntoselect = 3,
      npop = 8, niterations = 50, enable_restart = TRUE,
      restart_threshold = 0.3, max_restarts = 1,  # Force restart early
      minitbefstop = 5, tolconv = 1e-2, verbose = FALSE
    )
  })
  
  # Check restart tracking
  expect_true("restart_history" %in% names(ga_with_restart))
  expect_true("restart_count" %in% names(ga_with_restart))
  expect_true("total_generations" %in% names(ga_with_restart))
  expect_true(is.numeric(ga_with_restart$restart_count))
  expect_true(ga_with_restart$restart_count >= 0)
  
  # Test restart disabled
  suppressMessages({
    ga_no_restart <- subset_ga(
      P = P, Candidates = candidates, Test = test_set, ntoselect = 3,
      npop = 8, niterations = 20, enable_restart = FALSE,
      verbose = FALSE
    )
  })
  
  expect_equal(ga_no_restart$restart_count, 0)
  expect_equal(length(ga_no_restart$restart_history), 0)
})

test_that("Restart population generation works", {
  # Test the restart function directly
  set.seed(303)
  candidates <- paste0("ind", 1:20)
  current_pop <- list(
    paste0("ind", 1:5),
    paste0("ind", 3:7),
    paste0("ind", 8:12)
  )
  best_solution <- paste0("ind", 1:5)
  
  restart_result <- perform_restart(current_pop, best_solution, candidates, 5, 10, 20)
  
  expect_true(is.list(restart_result))
  expect_true("new_population" %in% names(restart_result))
  expect_true("strategy_breakdown" %in% names(restart_result))
  expect_equal(length(restart_result$new_population), 10)
  
  # Check that all solutions are valid
  for (solution in restart_result$new_population) {
    expect_equal(length(solution), 5)
    expect_true(all(solution %in% candidates))
    expect_equal(length(unique(solution)), 5)  # No duplicates
  }
})

test_that("Convergence function handles edge cases", {
  # Test with insufficient data
  short_history <- list(
    fitness = c(10.0, 11.0),
    diversity = c(0.5, 0.6),
    improvement_rate = c(0.1, 0.05),
    fitness_variance = c(0.2, 0.3)
  )
  
  result <- check_multi_criteria_convergence(
    short_history, 2, 10, 1e-6
  )
  
  expect_false(result$converged)
  expect_equal(result$reason, "Insufficient data")
  
  # Test with non-converged data
  volatile_history <- list(
    fitness = c(1:10),  # Continuously improving
    diversity = runif(10, 0, 1),  # Random diversity
    improvement_rate = runif(10, 0.1, 0.5),  # High improvement rates
    fitness_variance = runif(10, 1, 5)  # High variance
  )
  
  result2 <- check_multi_criteria_convergence(
    volatile_history, 10, 10, 1e-6
  )
  
  expect_false(result2$converged)
  expect_equal(result2$reason, "convergence criteria not met")
})