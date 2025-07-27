test_that("subset_ga_multiobjective basic functionality works", {
  # Create test data
  set.seed(123)
  n <- 30
  p <- 8
  Pcs <- matrix(rnorm(n * p), n, p)
  rownames(Pcs) <- paste0("ind", 1:n)
  colnames(Pcs) <- paste0("PC", 1:p)
  
  candidates <- rownames(Pcs)
  test_set <- paste0("ind", 21:25)
  ntoselect <- 12
  
  criteria <- c("pev_mean", "cd_mean")
  criteria_types <- c(TRUE, TRUE)  # Both minimize
  plot_directions <- c(1, 1)
  
  # Test basic functionality (few iterations for speed)
  result <- subset_ga_multiobjective(
    Pcs = Pcs, 
    candidates = candidates,
    test = test_set,
    ntoselect = ntoselect,
    criteria = criteria,
    criteria_types = criteria_types,
    plot_directions = plot_directions,
    npop = 10,
    niterations = 5,
    plot_iterations = FALSE
  )
  
  expect_type(result, "list")
  expect_true("pareto_solutions" %in% names(result))
  expect_true("pareto_fitness" %in% names(result))
  expect_true("all_archive" %in% names(result))
  expect_true("convergence_history" %in% names(result))
  expect_true("generation_stats" %in% names(result))
  expect_true("parameters" %in% names(result))
  
  # Check Pareto solutions
  expect_type(result$pareto_solutions, "list")
  expect_true(length(result$pareto_solutions) >= 1)
  
  # Check Pareto fitness
  expect_type(result$pareto_fitness, "double")
  expect_equal(ncol(result$pareto_fitness), 2)  # Two criteria
  expect_equal(nrow(result$pareto_fitness), length(result$pareto_solutions))
  
  # Check each solution has correct size
  for (sol in result$pareto_solutions) {
    expect_length(sol, ntoselect)
    expect_true(all(sol %in% candidates))
  }
  
  # Check convergence history
  expect_type(result$convergence_history, "double")
  expect_equal(ncol(result$convergence_history), 2)
  expect_true(nrow(result$convergence_history) <= 5)  # Up to 5 iterations
  
  # Check parameters
  expect_equal(result$parameters$npop, 10)
  expect_equal(result$parameters$criteria, criteria)
  expect_equal(result$parameters$criteria_types, criteria_types)
})

test_that("subset_ga_multiobjective_single works without test set", {
  # Create test data
  set.seed(123)
  n <- 25
  p <- 6
  Pcs <- matrix(rnorm(n * p), n, p)
  rownames(Pcs) <- paste0("ind", 1:n)
  
  candidates <- rownames(Pcs)
  ntoselect <- 10
  
  criteria <- c("a_optimality", "d_optimality")
  criteria_types <- c(TRUE, TRUE)  # Both minimize
  plot_directions <- c(1, 1)
  
  # Test without test set
  result <- subset_ga_multiobjective_single(
    Pcs = Pcs,
    candidates = candidates,
    ntoselect = ntoselect,
    criteria = criteria,
    criteria_types = criteria_types,
    plot_directions = plot_directions,
    npop = 8,
    niterations = 3
  )
  
  expect_type(result, "list")
  expect_true("pareto_solutions" %in% names(result))
  expect_true("pareto_fitness" %in% names(result))
  
  # Check that solutions are valid
  for (sol in result$pareto_solutions) {
    expect_length(sol, ntoselect)
    expect_true(all(sol %in% candidates))
  }
})

test_that("determine_primary_matrix works correctly", {
  # Create test matrices
  Pcs <- matrix(rnorm(20), 4, 5)
  Dist <- matrix(runif(16), 4, 4)
  Kernel <- matrix(runif(16, 0.8, 1.2), 4, 4)
  
  # Test with Pcs priority
  primary1 <- determine_primary_matrix(Pcs, Dist, Kernel)
  expect_equal(primary1, Pcs)
  
  # Test with Kernel when no Pcs
  primary2 <- determine_primary_matrix(NULL, Dist, Kernel)
  expect_equal(primary2, Kernel)
  
  # Test with Dist when no Pcs or Kernel
  primary3 <- determine_primary_matrix(NULL, Dist, NULL)
  expect_equal(primary3, Dist)
  
  # Test error when nothing provided
  expect_error(determine_primary_matrix(NULL, NULL, NULL))
})

test_that("initialize_population works correctly", {
  candidates <- paste0("ind", 1:50)
  ntoselect <- 15
  npop <- 10
  
  set.seed(123)
  population <- initialize_population(candidates, ntoselect, npop)
  
  expect_type(population, "list")
  expect_length(population, npop)
  
  # Check each solution
  for (i in 1:npop) {
    expect_length(population[[i]], ntoselect)
    expect_true(all(population[[i]] %in% candidates))
    expect_equal(length(unique(population[[i]])), ntoselect)  # No duplicates
  }
})

test_that("evaluate_multiobjective_population works correctly", {
  # Create test data
  set.seed(123)
  n <- 20
  p <- 6
  P <- matrix(rnorm(n * p), n, p)
  rownames(P) <- paste0("ind", 1:n)
  
  # Create population
  population <- list()
  for (i in 1:5) {
    population[[i]] <- sample(rownames(P), 8)
  }
  
  test_set <- paste0("ind", 16:18)
  criteria <- c("a_optimality", "pev_mean", "cd_mean")
  
  # Test evaluation
  fitness_matrix <- evaluate_multiobjective_population(population, P, test_set, criteria, 
                                                      lambda = 1e-6, C = NULL, mc.cores = 1)
  
  expect_type(fitness_matrix, "double")
  expect_equal(dim(fitness_matrix), c(5, 3))  # 5 solutions, 3 criteria
  expect_equal(colnames(fitness_matrix), criteria)
  expect_true(all(is.finite(fitness_matrix)))
  
  # Test without test set
  fitness_single <- evaluate_multiobjective_population(population, P, NULL, criteria,
                                                      lambda = 1e-6, C = NULL, mc.cores = 1)
  expect_equal(dim(fitness_single), c(5, 3))
})

test_that("non_dominated_sorting works correctly", {
  # Create test fitness matrix
  # Solution 1: (1, 3) - dominates solution 3
  # Solution 2: (2, 2) - non-dominated
  # Solution 3: (3, 4) - dominated by solution 1
  # Solution 4: (1.5, 2.5) - non-dominated
  fitness_matrix <- matrix(c(
    1, 3,
    2, 2, 
    3, 4,
    1.5, 2.5
  ), nrow = 4, byrow = TRUE)
  
  criteria_types <- c(TRUE, TRUE)  # Both minimize
  
  ranks <- non_dominated_sorting(fitness_matrix, criteria_types)
  
  expect_type(ranks, "double")
  expect_length(ranks, 4)
  expect_true(all(ranks >= 1))
  
  # Solutions 1, 2, and 4 should be in front 1
  expect_equal(ranks[1], 1)  # (1, 3)
  expect_equal(ranks[2], 1)  # (2, 2)
  expect_equal(ranks[4], 1)  # (1.5, 2.5)
  
  # Solution 3 should be in front 2 (dominated)
  expect_equal(ranks[3], 2)  # (3, 4)
  
  # Test with maximization
  criteria_types_max <- c(FALSE, FALSE)  # Both maximize
  ranks_max <- non_dominated_sorting(fitness_matrix, criteria_types_max)
  
  # Now solution 3 (3, 4) should dominate solution 1 (1, 3)
  expect_equal(ranks_max[3], 1)  # (3, 4) is best for maximization
  expect_equal(ranks_max[1], 2)  # (1, 3) is dominated
})

test_that("check_dominance works correctly", {
  # Test minimization
  criteria_types <- c(TRUE, TRUE)  # Both minimize
  
  # Solution 1 dominates solution 2
  sol1 <- c(1, 2)
  sol2 <- c(2, 3)
  expect_equal(check_dominance(sol1, sol2, criteria_types), 1)
  expect_equal(check_dominance(sol2, sol1, criteria_types), -1)
  
  # Non-dominated solutions
  sol3 <- c(1, 3)
  sol4 <- c(2, 2)
  expect_equal(check_dominance(sol3, sol4, criteria_types), 0)
  expect_equal(check_dominance(sol4, sol3, criteria_types), 0)
  
  # Identical solutions
  expect_equal(check_dominance(sol1, sol1, criteria_types), 0)
  
  # Test maximization
  criteria_types_max <- c(FALSE, FALSE)  # Both maximize
  expect_equal(check_dominance(sol1, sol2, criteria_types_max), -1)  # sol2 dominates sol1
  expect_equal(check_dominance(sol2, sol1, criteria_types_max), 1)
  
  # Test mixed criteria
  criteria_mixed <- c(TRUE, FALSE)  # Minimize first, maximize second
  sol5 <- c(1, 1)  # Low first (good), low second (bad)
  sol6 <- c(2, 2)  # High first (bad), high second (good)
  expect_equal(check_dominance(sol5, sol6, criteria_mixed), 0)  # Non-dominated
})

test_that("find_non_dominated works correctly", {
  # Create test fitness matrix
  fitness_matrix <- matrix(c(
    1, 3,    # Non-dominated
    2, 2,    # Non-dominated
    3, 4,    # Dominated
    1.5, 2.5 # Non-dominated
  ), nrow = 4, byrow = TRUE)
  
  criteria_types <- c(TRUE, TRUE)  # Both minimize
  
  non_dom_indices <- find_non_dominated(fitness_matrix, criteria_types)
  
  expect_type(non_dom_indices, "integer")
  expect_true(all(non_dom_indices %in% 1:4))
  expect_true(1 %in% non_dom_indices)  # (1, 3)
  expect_true(2 %in% non_dom_indices)  # (2, 2)
  expect_true(4 %in% non_dom_indices)  # (1.5, 2.5)
  expect_false(3 %in% non_dom_indices) # (3, 4) is dominated
})

test_that("calculate_crowding_distance works correctly", {
  # Create test fitness matrix for a single front
  fitness_front <- matrix(c(
    1, 4,    # Boundary point
    2, 3,    # Middle point
    3, 2,    # Middle point
    4, 1     # Boundary point
  ), nrow = 4, byrow = TRUE)
  
  crowding_dist <- calculate_crowding_distance(fitness_front)
  
  expect_type(crowding_dist, "double")
  expect_length(crowding_dist, 4)
  expect_true(all(crowding_dist >= 0))
  
  # Boundary points should have infinite distance
  expect_equal(crowding_dist[1], Inf)  # Boundary in both objectives
  expect_equal(crowding_dist[4], Inf)  # Boundary in both objectives
  
  # Test with small front
  small_front <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  small_crowding <- calculate_crowding_distance(small_front)
  expect_true(all(small_crowding == Inf))  # All points are boundary points
  
  # Test with single point
  single_front <- matrix(c(1, 2), nrow = 1)
  single_crowding <- calculate_crowding_distance(single_front)
  expect_equal(single_crowding, Inf)
})

test_that("multiobjective_selection works correctly", {
  # Create test data
  fitness_matrix <- matrix(c(
    1, 3,    # Front 1
    2, 2,    # Front 1
    3, 4,    # Front 2
    1.5, 2.5,# Front 1
    4, 5     # Front 2
  ), nrow = 5, byrow = TRUE)
  
  criteria_types <- c(TRUE, TRUE)  # Both minimize
  ranks <- non_dominated_sorting(fitness_matrix, criteria_types)
  
  # Select top 3
  selected <- multiobjective_selection(fitness_matrix, ranks, n_select = 3)
  
  expect_type(selected, "integer")
  expect_length(selected, 3)
  expect_true(all(selected %in% 1:5))
  
  # All selected should be from front 1 (since front 1 has 3 solutions)
  selected_ranks <- ranks[selected]
  expect_true(all(selected_ranks == 1))
  
  # Test selecting more than front 1
  selected_large <- multiobjective_selection(fitness_matrix, ranks, n_select = 4)
  expect_length(selected_large, 4)
  
  # Should include all of front 1 and one from front 2
  selected_large_ranks <- ranks[selected_large]
  expect_equal(sum(selected_large_ranks == 1), 3)  # All front 1
  expect_equal(sum(selected_large_ranks == 2), 1)  # One from front 2
})

test_that("generate_multiobjective_offspring works correctly", {
  set.seed(123)
  candidates <- paste0("ind", 1:40)
  
  # Create parent population
  parents <- list()
  for (i in 1:6) {
    parents[[i]] <- sample(candidates, 12)
  }
  
  ntoselect <- 12
  n_offspring <- 10
  mutprob <- 0.4
  
  # Test offspring generation
  offspring <- generate_multiobjective_offspring(parents, candidates, ntoselect, 
                                                 n_offspring, mutprob, 
                                                 mutintensity = 1, mc.cores = 1)
  
  expect_type(offspring, "list")
  expect_length(offspring, n_offspring)
  
  # Check each offspring
  for (i in 1:n_offspring) {
    expect_length(offspring[[i]], ntoselect)
    expect_true(all(offspring[[i]] %in% candidates))
    expect_equal(length(unique(offspring[[i]])), ntoselect)  # No duplicates
  }
  
  # Test with high mutation probability
  offspring_high_mut <- generate_multiobjective_offspring(parents, candidates, ntoselect,
                                                          5, mutprob = 0.9,
                                                          mutintensity = 2, mc.cores = 1)
  expect_length(offspring_high_mut, 5)
})

test_that("update_pareto_archive works correctly", {
  set.seed(123)
  candidates <- paste0("ind", 1:30)
  
  # Create current archive
  archive <- list()
  archive[[1]] <- sample(candidates, 10)
  archive[[2]] <- sample(candidates, 10)
  
  archive_fitness <- matrix(c(
    1, 3,
    2, 2
  ), nrow = 2, byrow = TRUE)
  
  # Create new population
  population <- list()
  population[[1]] <- sample(candidates, 10)  # New solution
  population[[2]] <- archive[[1]]            # Duplicate of archive
  
  population_fitness <- matrix(c(
    0.5, 3.5,  # New non-dominated solution
    1, 3       # Duplicate
  ), nrow = 2, byrow = TRUE)
  
  criteria_types <- c(TRUE, TRUE)  # Both minimize
  
  # Test archive update
  result <- update_pareto_archive(population, population_fitness, archive, 
                                 archive_fitness, max_size = 10, criteria_types)
  
  expect_type(result, "list")
  expect_true("archive" %in% names(result))
  expect_true("fitness" %in% names(result))
  
  # Archive should contain non-dominated solutions
  expect_true(length(result$archive) >= 2)
  expect_equal(nrow(result$fitness), length(result$archive))
  
  # Test archive size limit
  large_population <- list()
  large_fitness <- matrix(0, 15, 2)
  for (i in 1:15) {
    large_population[[i]] <- sample(candidates, 10)
    large_fitness[i, ] <- c(runif(1, 0, 2), runif(1, 0, 2))
  }
  
  result_limited <- update_pareto_archive(large_population, large_fitness, list(),
                                         matrix(0, 0, 2), max_size = 5, criteria_types)
  
  expect_true(length(result_limited$archive) <= 5)
  expect_equal(nrow(result_limited$fitness), length(result_limited$archive))
})

test_that("legacy wrapper functions work correctly", {
  # Test GenAlgForSubsetSelectionMO wrapper
  set.seed(123)
  n <- 20
  p <- 6
  Pcs <- matrix(rnorm(n * p), n, p)
  rownames(Pcs) <- paste0("ind", 1:n)
  
  candidates <- rownames(Pcs)
  test_set <- paste0("ind", 16:18)
  ntoselect <- 8
  
  # Test legacy function
  result_legacy <- GenAlgForSubsetSelectionMO(
    Pcs = Pcs,
    candidates = candidates,
    Test = test_set,
    ntoselect = ntoselect,
    selectionstats = c("a_optimality", "pev_mean"),
    selectionstatstypes = c(TRUE, TRUE),
    plotdirections = c(1, 1),
    npopGA = 6,
    nitGA = 3,
    plotiters = FALSE
  )
  
  expect_type(result_legacy, "list")
  expect_true("pareto_solutions" %in% names(result_legacy))
  
  # Test GenAlgForSubsetSelectionMONoTest wrapper
  result_notest <- GenAlgForSubsetSelectionMONoTest(
    Pcs = Pcs,
    Candidates = candidates,
    ntoselect = ntoselect,
    selectionstats = c("a_optimality", "d_optimality"),
    selectionstatstypes = c(TRUE, TRUE),
    plotdirections = c(1, 1),
    npopGA = 5,
    nitGA = 2
  )
  
  expect_type(result_notest, "list")
  expect_true("pareto_solutions" %in% names(result_notest))
})

test_that("error handling works correctly", {
  # Test with no matrices provided
  expect_error(subset_ga_multiobjective(
    Pcs = NULL, Dist = NULL, Kernel = NULL,
    candidates = c("a", "b"), test = "c", ntoselect = 1,
    criteria = "a_optimality", criteria_types = TRUE, plot_directions = 1
  ))
  
  # Test with mismatched criteria and types lengths
  Pcs <- matrix(rnorm(20), 4, 5)
  rownames(Pcs) <- paste0("ind", 1:4)
  
  expect_error(subset_ga_multiobjective(
    Pcs = Pcs,
    candidates = rownames(Pcs),
    test = NULL,
    ntoselect = 2,
    criteria = c("a_optimality", "d_optimality"),
    criteria_types = c(TRUE),  # Wrong length
    plot_directions = c(1, 1)
  ))
  
  # Test non_dominated_sorting with empty matrix
  empty_fitness <- matrix(numeric(0), 0, 2)
  expect_length(non_dominated_sorting(empty_fitness, c(TRUE, TRUE)), 0)
  
  # Test invalid criterion in multiobjective evaluation
  population <- list(c("ind1", "ind2"))
  P <- matrix(rnorm(8), 4, 2)
  rownames(P) <- paste0("ind", 1:4)
  
  expect_warning(
    fitness_invalid <- evaluate_multiobjective_population(
      population, P, NULL, c("INVALID"), 1e-6, NULL, 1
    )
  )
  expect_true(all(is.infinite(fitness_invalid)))
})

test_that("numerical properties are maintained", {
  # Test that Pareto fronts maintain dominance relationships
  set.seed(123)
  fitness_matrix <- matrix(runif(20, 0, 10), 10, 2)
  criteria_types <- c(TRUE, TRUE)
  
  ranks <- non_dominated_sorting(fitness_matrix, criteria_types)
  
  # Check that no solution in front i dominates any solution in front j where i > j
  for (i in 1:(max(ranks)-1)) {
    for (j in (i+1):max(ranks)) {
      front_i_indices <- which(ranks == i)
      front_j_indices <- which(ranks == j)
      
      for (ii in front_i_indices) {
        for (jj in front_j_indices) {
          dominance <- check_dominance(fitness_matrix[ii, ], fitness_matrix[jj, ], criteria_types)
          expect_true(dominance != -1)  # Solution in better front should not be dominated
        }
      }
    }
  }
  
  # Test crowding distance properties
  if (sum(ranks == 1) > 2) {  # If front 1 has more than 2 solutions
    front1_indices <- which(ranks == 1)
    front1_fitness <- fitness_matrix[front1_indices, , drop = FALSE]
    crowding_dist <- calculate_crowding_distance(front1_fitness)
    
    # All distances should be non-negative
    expect_true(all(crowding_dist >= 0))
    
    # At least some boundary points should have infinite distance
    if (nrow(front1_fitness) > 2) {
      expect_true(any(is.infinite(crowding_dist)))
    }
  }
})

test_that("multiobjective GA converges appropriately", {
  # Test that the algorithm produces reasonable results
  set.seed(42)  # Fixed seed for reproducibility
  
  n <- 25
  p <- 5
  Pcs <- matrix(rnorm(n * p), n, p)
  rownames(Pcs) <- paste0("ind", 1:n)
  
  candidates <- rownames(Pcs)
  test_set <- paste0("ind", 21:23)
  ntoselect <- 10
  
  criteria <- c("a_optimality", "d_optimality")
  criteria_types <- c(TRUE, TRUE)
  plot_directions <- c(1, 1)
  
  # Run with more iterations to test convergence
  result <- subset_ga_multiobjective(
    Pcs = Pcs,
    candidates = candidates,
    test = test_set,
    ntoselect = ntoselect,
    criteria = criteria,
    criteria_types = criteria_types,
    plot_directions = plot_directions,
    npop = 15,
    niterations = 10,
    mutprob = 0.6,
    plot_iterations = FALSE
  )
  
  # Check that we found multiple Pareto solutions
  expect_true(length(result$pareto_solutions) >= 1)
  expect_true(nrow(result$pareto_fitness) >= 1)
  
  # Check that solutions are truly non-dominated
  if (nrow(result$pareto_fitness) > 1) {
    for (i in 1:(nrow(result$pareto_fitness)-1)) {
      for (j in (i+1):nrow(result$pareto_fitness)) {
        dominance1 <- check_dominance(result$pareto_fitness[i, ], 
                                     result$pareto_fitness[j, ], criteria_types)
        dominance2 <- check_dominance(result$pareto_fitness[j, ], 
                                     result$pareto_fitness[i, ], criteria_types)
        expect_equal(dominance1, 0)  # Should be non-dominated
        expect_equal(dominance2, 0)  # Should be non-dominated
      }
    }
  }
  
  # Check that all solutions are valid
  for (sol in result$pareto_solutions) {
    expect_length(sol, ntoselect)
    expect_true(all(sol %in% candidates))
    expect_equal(length(unique(sol)), ntoselect)
  }
  
  # Check convergence history
  expect_true(nrow(result$convergence_history) <= 10)
  expect_equal(ncol(result$convergence_history), 2)
  expect_true(all(is.finite(result$convergence_history)))
})