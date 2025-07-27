test_that("crossover function works correctly", {
  set.seed(123)
  candidates <- paste0("ind", 1:50)
  parent1 <- sample(candidates, 10)
  parent2 <- sample(candidates, 10)
  ntoselect <- 12
  
  # Test uniform crossover
  offspring_uniform <- crossover(parent1, parent2, ntoselect, candidates, "uniform")
  expect_length(offspring_uniform, ntoselect)
  expect_true(all(offspring_uniform %in% candidates))
  expect_equal(length(unique(offspring_uniform)), ntoselect)  # No duplicates
  
  # Test two-point crossover
  offspring_twopoint <- crossover(parent1, parent2, ntoselect, candidates, "twopoint")
  expect_length(offspring_twopoint, ntoselect)
  expect_true(all(offspring_twopoint %in% candidates))
  
  # Test proportional crossover
  offspring_prop <- crossover(parent1, parent2, ntoselect, candidates, "proportional")
  expect_length(offspring_prop, ntoselect)
  expect_true(all(offspring_prop %in% candidates))
  
  # Test with bias
  offspring_biased <- crossover(parent1, parent2, ntoselect, candidates, "uniform", bias = 0.8)
  expect_length(offspring_biased, ntoselect)
  expect_true(all(offspring_biased %in% candidates))
  
  # Test edge cases
  # Empty parents
  offspring_empty <- crossover(character(0), character(0), 5, candidates, "uniform")
  expect_length(offspring_empty, 5)
  expect_true(all(offspring_empty %in% candidates))
  
  # One parent empty
  offspring_half_empty <- crossover(parent1, character(0), 8, candidates, "uniform")
  expect_length(offspring_half_empty, 8)
})

test_that("generate_offspring function works correctly", {
  set.seed(123)
  candidates <- paste0("ind", 1:100)
  
  # Create elite population
  elites <- list()
  for (i in 1:5) {
    elites[[i]] <- sample(candidates, 15)
  }
  
  npop <- 20
  mutprob <- 0.5
  ntoselect <- 15
  
  # Test basic functionality
  offspring <- generate_offspring(elites, candidates, npop, mutprob, 
                                  ntoselect = ntoselect)
  
  expect_type(offspring, "list")
  expect_length(offspring, npop)
  
  # Check each offspring
  for (i in 1:npop) {
    expect_length(offspring[[i]], ntoselect)
    expect_true(all(offspring[[i]] %in% candidates))
    expect_equal(length(unique(offspring[[i]])), ntoselect)  # No duplicates
  }
  
  # Test with tabu memory
  tabu_memory <- list(list(sample(candidates, 15)))
  offspring_tabu <- generate_offspring(elites, candidates, 5, mutprob,
                                       tabu_memory = tabu_memory,
                                       ntoselect = ntoselect)
  expect_length(offspring_tabu, 5)
  
  # Test with different crossover methods
  offspring_twopoint <- generate_offspring(elites, candidates, 10, mutprob,
                                           crossover_method = "twopoint",
                                           ntoselect = ntoselect)
  expect_length(offspring_twopoint, 10)
  
  # Test parallel processing
  offspring_parallel <- generate_offspring(elites, candidates, 60, mutprob,
                                           mc.cores = 1, # Keep as 1 for testing
                                           parallel_threshold = 50,
                                           ntoselect = ntoselect)
  expect_length(offspring_parallel, 60)
})

test_that("mutate_solution function works correctly", {
  set.seed(123)
  candidates <- paste0("ind", 1:50)
  solution <- sample(candidates, 15)
  ntoselect <- length(solution)
  
  # Test light mutation
  mutated_light <- mutate_solution(solution, candidates, 0.1, ntoselect)
  expect_length(mutated_light, ntoselect)
  expect_true(all(mutated_light %in% candidates))
  
  # Test heavy mutation
  mutated_heavy <- mutate_solution(solution, candidates, 0.8, ntoselect)
  expect_length(mutated_heavy, ntoselect)
  expect_true(all(mutated_heavy %in% candidates))
  
  # Test with integer intensity
  mutated_int <- mutate_solution(solution, candidates, 3, ntoselect)
  expect_length(mutated_int, ntoselect)
  expect_true(all(mutated_int %in% candidates))
  
  # Check that some mutation occurred with high intensity
  mutated_max <- mutate_solution(solution, candidates, 1.0, ntoselect)
  expect_true(any(!mutated_max %in% solution))  # At least some change
})

test_that("calculate_diversity function works correctly", {
  set.seed(123)
  candidates <- paste0("ind", 1:50)
  
  # Create population with varying diversity
  population_diverse <- list()
  for (i in 1:10) {
    population_diverse[[i]] <- sample(candidates, 10)
  }
  
  diversity_diverse <- calculate_diversity(population_diverse)
  expect_type(diversity_diverse, "double")
  expect_length(diversity_diverse, 1)
  expect_true(diversity_diverse >= 0 && diversity_diverse <= 1)
  
  # Create identical population (low diversity)
  base_solution <- sample(candidates, 10)
  population_identical <- rep(list(base_solution), 5)
  diversity_low <- calculate_diversity(population_identical)
  expect_equal(diversity_low, 0)
  
  # Single solution
  population_single <- list(sample(candidates, 10))
  diversity_single <- calculate_diversity(population_single)
  expect_equal(diversity_single, 0)
  
  # Empty population
  population_empty <- list()
  diversity_empty <- calculate_diversity(population_empty)
  expect_equal(diversity_empty, 0)
})

test_that("adaptive_mutation_rate function works correctly", {
  base_rate <- 0.5
  max_generations <- 100
  
  # Test early generation, low diversity
  rate_early_low <- adaptive_mutation_rate(base_rate, 0.1, 10, max_generations)
  expect_type(rate_early_low, "double")
  expect_true(rate_early_low >= 0.01 && rate_early_low <= 0.95)
  expect_true(rate_early_low > base_rate)  # Should increase with low diversity
  
  # Test mid generation, high diversity
  rate_mid_high <- adaptive_mutation_rate(base_rate, 0.8, 50, max_generations)
  expect_true(rate_mid_high <= rate_early_low)  # Should be lower with high diversity
  
  # Test late generation
  rate_late <- adaptive_mutation_rate(base_rate, 0.5, 90, max_generations)
  expect_true(rate_late > base_rate)  # Should increase late in evolution
  
  # Test bounds
  rate_bounded <- adaptive_mutation_rate(0.9, 0.1, 10, max_generations, 
                                         min_rate = 0.1, max_rate = 0.8)
  expect_true(rate_bounded >= 0.1 && rate_bounded <= 0.8)
})

test_that("tournament_selection function works correctly", {
  set.seed(123)
  n_pop <- 20
  fitness <- rnorm(n_pop)  # Random fitness values
  
  # Create dummy population
  population <- list()
  for (i in 1:n_pop) {
    population[[i]] <- paste0("solution_", i)
  }
  
  # Test basic tournament selection
  selected <- tournament_selection(population, fitness, tournament_size = 3, n_select = 5)
  expect_type(selected, "integer")
  expect_length(selected, 5)
  expect_true(all(selected >= 1 && selected <= n_pop))
  
  # Test single selection
  selected_one <- tournament_selection(population, fitness, tournament_size = 5, n_select = 1)
  expect_length(selected_one, 1)
  
  # Test that better fitness values are more likely to be selected
  # Create population with clear best and worst
  fitness_clear <- c(rep(-10, 5), rep(10, 15))  # Lower is better (minimization)
  selected_biased <- tournament_selection(population, fitness_clear, 
                                          tournament_size = 10, n_select = 10)
  
  # Should select more from the better half (indices 1-5)
  better_selected <- sum(selected_biased <= 5)
  expect_true(better_selected > 2)  # At least some bias toward better solutions
})

test_that("elite_selection function works correctly", {
  set.seed(123)
  n_pop <- 15
  fitness <- rnorm(n_pop)
  
  population <- list()
  for (i in 1:n_pop) {
    population[[i]] <- paste0("solution_", i)
  }
  
  # Test elite selection
  selected <- elite_selection(population, fitness, n_select = 5)
  expect_type(selected, "integer")
  expect_length(selected, 5)
  expect_true(all(selected >= 1 && selected <= n_pop))
  
  # Check that selected indices correspond to best fitness
  selected_fitness <- fitness[selected]
  expect_true(all(selected_fitness <= sort(fitness)[5]))  # All should be in top 5
  
  # Test edge case - select all
  all_selected <- elite_selection(population, fitness, n_pop)
  expect_equal(sort(all_selected), 1:n_pop)
})

test_that("crowding_replacement function works correctly", {
  set.seed(123)
  candidates <- paste0("ind", 1:50)
  
  # Create population and offspring
  population <- list()
  population_fitness <- numeric(10)
  for (i in 1:10) {
    population[[i]] <- sample(candidates, 12)
    population_fitness[i] <- runif(1)
  }
  
  offspring <- list()
  offspring_fitness <- numeric(3)
  for (i in 1:3) {
    offspring[[i]] <- sample(candidates, 12)
    offspring_fitness[i] <- runif(1, min = 0, max = 0.5)  # Better fitness
  }
  
  # Test crowding replacement
  result <- crowding_replacement(offspring, offspring_fitness, population, 
                                population_fitness, crowding_factor = 3)
  
  expect_type(result, "list")
  expect_true("population" %in% names(result))
  expect_true("fitness" %in% names(result))
  expect_length(result$population, 10)
  expect_length(result$fitness, 10)
  
  # Check that at least some offspring were integrated
  # (difficult to test precisely due to similarity calculations)
  expect_true(all(sapply(result$population, length) == 12))
})

test_that("fitness_sharing function works correctly", {
  set.seed(123)
  candidates <- paste0("ind", 1:30)
  
  # Create population with some similar solutions
  population <- list()
  fitness <- runif(8)
  
  # Create some similar solutions
  base_solution <- sample(candidates, 10)
  population[[1]] <- base_solution
  population[[2]] <- c(base_solution[1:8], sample(setdiff(candidates, base_solution), 2))
  
  # Add diverse solutions
  for (i in 3:8) {
    population[[i]] <- sample(candidates, 10)
  }
  
  # Test fitness sharing
  shared_fitness <- fitness_sharing(fitness, population, sharing_radius = 0.3, alpha = 1)
  
  expect_type(shared_fitness, "double")
  expect_length(shared_fitness, 8)
  expect_true(all(shared_fitness > 0))
  
  # Similar solutions should have reduced fitness
  expect_true(shared_fitness[1] > fitness[1])  # Fitness should be penalized
  expect_true(shared_fitness[2] > fitness[2])  # Fitness should be penalized
  
  # Test with no sharing (large radius)
  no_sharing <- fitness_sharing(fitness, population, sharing_radius = 2.0)
  expect_equal(no_sharing, fitness)  # Should be unchanged
})

test_that("edge cases and error handling work", {
  candidates <- paste0("ind", 1:20)
  
  # Test crossover with insufficient candidates
  expect_error(crossover(candidates[1:10], candidates[11:15], 25, candidates[1:15]))
  
  # Test generate_offspring with no elites
  expect_error(generate_offspring(list(), candidates, 5, 0.5))
  
  # Test empty candidate pool
  expect_error(crossover(candidates[1:5], candidates[6:10], 10, character(0)))
  
  # Test mutate_solution with no available candidates
  solution <- candidates[1:10]
  remaining_candidates <- candidates[11:20]
  # This should work but may not change much
  mutated <- mutate_solution(solution, remaining_candidates, 0.5, 10)
  expect_length(mutated, 10)
})

test_that("deterministic behavior with set seed", {
  # Test that functions give reproducible results
  set.seed(42)
  candidates <- paste0("ind", 1:50)
  parent1 <- sample(candidates, 10)
  parent2 <- sample(candidates, 10)
  
  set.seed(42)
  result1 <- crossover(parent1, parent2, 12, candidates, "uniform")
  
  set.seed(42)
  result2 <- crossover(parent1, parent2, 12, candidates, "uniform")
  
  expect_equal(result1, result2)
  
  # Test with generate_offspring
  elites <- list(sample(candidates, 10), sample(candidates, 10))
  
  set.seed(123)
  offspring1 <- generate_offspring(elites, candidates, 5, 0.3, mc.cores = 1)
  
  set.seed(123)
  offspring2 <- generate_offspring(elites, candidates, 5, 0.3, mc.cores = 1)
  
  expect_equal(offspring1, offspring2)
})