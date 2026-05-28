# Test selection mechanisms and validation
library(STPGA)

test_that("Rank-based selection works correctly", {
  # Create test population and fitness
  set.seed(123)
  population <- list(
    paste0("ind", 1:5),
    paste0("ind", 2:6),
    paste0("ind", 3:7),
    paste0("ind", 4:8),
    paste0("ind", 5:9)
  )
  fitness <- c(10, 5, 15, 8, 12)  # fitness[2] is best
  
  # Test linear rank selection
  selected <- rank_selection(population, fitness, 3, selection_pressure = 1.5, "linear")
  
  expect_equal(length(selected), 3)
  expect_true(all(selected %in% 1:5))
  
  # Test exponential rank selection
  selected_exp <- rank_selection(population, fitness, 2, selection_pressure = 1.2, "exponential")
  
  expect_equal(length(selected_exp), 2)
  expect_true(all(selected_exp %in% 1:5))
})

test_that("Selection parameter validation works", {
  # Test valid parameters
  result1 <- validate_selection_parameters("tournament", 3, 1.5, 20)
  expect_equal(result1$selection_method, "tournament")
  expect_equal(result1$tournament_size, 3)
  expect_equal(length(result1$warnings), 0)
  
  # Test invalid selection method
  suppressWarnings({
    result2 <- validate_selection_parameters("invalid_method", 3, 1.5, 20)
  })
  expect_equal(result2$selection_method, "tournament")
  expect_true(length(result2$warnings) > 0)
  
  # Test invalid tournament size
  result3 <- validate_selection_parameters("tournament", 25, 1.5, 20)
  expect_equal(result3$tournament_size, 20)  # Adjusted to population size
  expect_true(any(grepl("Tournament size", result3$warnings)))
  
  # Test invalid selection pressure
  result4 <- validate_selection_parameters("rank", 3, 3.0, 20)
  expect_equal(result4$selection_pressure, 2.0)  # Adjusted to max
  expect_true(any(grepl("Selection pressure", result4$warnings)))
})

test_that("GA with rank-based selection works", {
  # Create test data
  set.seed(456)
  n_ind <- 20
  n_markers <- 8
  
  P <- matrix(rnorm(n_ind * n_markers), nrow = n_ind, ncol = n_markers)
  rownames(P) <- paste0("ind", 1:n_ind)
  
  candidates <- rownames(P)[1:15]
  test_set <- rownames(P)[16:20]
  
  # Test rank selection
  suppressMessages({
    ga_rank <- subset_ga(
      P = P, Candidates = candidates, Test = test_set, ntoselect = 5,
      npop = 10, niterations = 15, selection_method = "rank",
      selection_pressure = 1.3, verbose = FALSE
    )
  })
  
  expect_true("best_solution" %in% names(ga_rank))
  expect_equal(ga_rank$parameters$selection_method, "rank")
  expect_equal(ga_rank$parameters$selection_pressure, 1.3)
  
  # Test hybrid selection
  suppressMessages({
    ga_hybrid <- subset_ga(
      P = P, Candidates = candidates, Test = test_set, ntoselect = 5,
      npop = 12, niterations = 10, selection_method = "hybrid",
      nelite = 6, verbose = FALSE
    )
  })
  
  expect_equal(ga_hybrid$parameters$selection_method, "hybrid")
})

test_that("Selection pressure boundary conditions", {
  # Test boundary values
  set.seed(789)
  population <- list(
    paste0("ind", 1:3),
    paste0("ind", 2:4),
    paste0("ind", 3:5)
  )
  fitness <- c(1, 2, 3)
  
  # Test minimum selection pressure
  suppressWarnings({
    selected_min <- rank_selection(population, fitness, 2, selection_pressure = 1.0)
  })
  expect_equal(length(selected_min), 2)
  
  # Test maximum selection pressure
  suppressWarnings({
    selected_max <- rank_selection(population, fitness, 2, selection_pressure = 2.0)
  })
  expect_equal(length(selected_max), 2)
  
  # Test out-of-bounds values (should be adjusted with warning)
  suppressWarnings({
    selected_low <- rank_selection(population, fitness, 2, selection_pressure = 0.5)
    selected_high <- rank_selection(population, fitness, 2, selection_pressure = 3.0)
  })
  expect_equal(length(selected_low), 2)
  expect_equal(length(selected_high), 2)
})

test_that("Selection methods produce different behaviors", {
  # Create population with clear fitness differences
  set.seed(101)
  population <- lapply(1:10, function(i) paste0("ind", i:(i+3)))
  fitness <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)  # Clear gradient
  
  n_select <- 5
  n_trials <- 100
  
  # Test elite selection (should always select same best individuals)
  elite_selections <- replicate(n_trials, {
    elite_selection(population, fitness, n_select)
  })
  elite_variance <- apply(elite_selections, 1, var)
  expect_true(all(elite_variance == 0))  # Should be deterministic
  
  # Test tournament selection (should have some randomness)
  tournament_selections <- replicate(n_trials, {
    tournament_selection(population, fitness, 3, n_select)
  })
  tournament_variance <- apply(tournament_selections, 1, var)
  expect_true(any(tournament_variance > 0))  # Should have some variance
  
  # Test rank selection (should prefer good individuals but allow diversity)
  rank_selections <- replicate(n_trials, {
    rank_selection(population, fitness, n_select, 1.5)
  })
  rank_variance <- apply(rank_selections, 1, var)
  expect_true(any(rank_variance > 0))  # Should have some variance
})