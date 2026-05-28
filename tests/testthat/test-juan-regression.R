library(STPGA)

test_that("Juan regression: cd_mean GA runs when selecting all candidates", {
  set.seed(260528)
  n <- 60
  p <- 8
  P <- matrix(rnorm(n * p), nrow = n, ncol = p)
  rownames(P) <- paste0("ind", seq_len(n))

  Candidates <- rownames(P)[1:30]
  Test <- rownames(P)[31:60]

  res <- subset_ga(
    P = P,
    Candidates = Candidates,
    Test = Test,
    ntoselect = 30,
    criterion = "cd_mean",
    npop = 12,
    niterations = 8,
    minitbefstop = 20,
    verbose = FALSE
  )

  expect_true(is.list(res))
  expect_equal(length(res$best_solution), 30)
  expect_true(setequal(res$best_solution, Candidates))
  expect_true(is.finite(res$best_fitness))
})

test_that("Juan regression: distance criteria accept MDS coordinates and D matrix", {
  set.seed(260529)
  n <- 50
  coords <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  rownames(coords) <- paste0("ind", seq_len(n))
  D <- as.matrix(dist(coords))

  Candidates <- rownames(coords)[1:25]
  Test <- rownames(coords)[26:50]

  direct_max <- distance_train_to_test_max(Candidates, Test, D)
  coord_max <- distance_train_to_test_max(Candidates, Test, coords)
  direct_mean <- distance_train_to_test_mean(Candidates, Test, D)
  coord_mean <- distance_train_to_test_mean(Candidates, Test, coords)

  expect_equal(coord_max, direct_max)
  expect_equal(coord_mean, direct_mean)
  expect_lte(coord_mean, coord_max)

  res <- subset_ga(
    P = coords,
    Candidates = Candidates,
    Test = Test,
    ntoselect = 10,
    criterion = "max_to_test",
    npop = 10,
    niterations = 6,
    minitbefstop = 20,
    verbose = FALSE
  )

  expect_true(is.finite(res$best_fitness))
  expect_equal(length(res$best_solution), 10)
  expect_true(all(res$best_solution %in% Candidates))
})

test_that("Juan regression: subset_ga_single preserves supplied Candidates and Test", {
  set.seed(260530)
  n <- 45
  coords <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  rownames(coords) <- paste0("ind", seq_len(n))

  Candidates <- rownames(coords)[1:20]
  Test <- rownames(coords)[21:45]

  res <- subset_ga_single(
    P = coords,
    Candidates = Candidates,
    Test = Test,
    ntoselect = 10,
    criterion = "mean_to_test",
    npop = 10,
    niterations = 6,
    minitbefstop = 20,
    verbose = FALSE
  )

  expect_true(is.finite(res$best_fitness))
  expect_equal(length(res$best_solution), 10)
  expect_true(all(res$best_solution %in% Candidates))
})

test_that("Juan regression: no-test cd_mean returns immediately", {
  set.seed(260531)
  P <- matrix(rnorm(60 * 6), nrow = 60, ncol = 6)
  rownames(P) <- paste0("ind", seq_len(60))

  res <- subset_ga_single(
    P = P,
    ntoselect = 20,
    criterion = "cd_mean",
    npop = 20,
    niterations = 50,
    verbose = FALSE
  )

  expect_equal(res$convergence_generation, 0)
  expect_equal(res$parameters$niterations, 0)
  expect_equal(length(res$best_solution), 20)
  expect_true(is.finite(res$best_fitness))
})

test_that("selection helpers are available through the installed namespace", {
  expect_true(exists("validate_selection_parameters", envir = asNamespace("STPGA")))
  expect_true(exists("rank_selection", envir = asNamespace("STPGA")))
  expect_true("validate_selection_parameters" %in% getNamespaceExports("STPGA"))
  expect_true("rank_selection" %in% getNamespaceExports("STPGA"))
})
