# STPGA v8 Development Session — Starter Prompt

## Context

STPGA is an R package for **Selection of Training Populations by Genetic Algorithm**.
Author: Deniz Akdemir (deniz.akdemir.work@gmail.com)
GitHub: https://github.com/denizakdemir/STPGA
Current CRAN version: 5.2.1
Current dev version: 7.0.1 (not yet on CRAN)
Source: `~/Dropbox/MasterFolder/RESEARCH/ModernClinicalDesignbyOptimization/SourceCode/STPGA/`

The package selects an optimal training population subset (from a pool of candidates)
to best predict a held-out test set, using genetic algorithm optimization over a
user-specified criterion (PEV, CD, A/D/E/G/I-optimality, distance-based, mixed-model).

## Goal for this session

Develop **STPGA v8.0.0**: a major quality release. **No interface changes** — the
public API (`subset_ga`, `subset_ga_single`, `criterion`, `GenAlgForSubsetSelection*`)
must remain fully backward-compatible. All improvements are internal.

### Priority improvements

#### 1. Computational efficiency
- Profile and optimize the inner GA loop (`ga-core.R`), especially `evaluate_population()`
- The evaluation cache in `evaluation.R` uses string concatenation for keys — benchmark
  against `digest::digest()` hashing
- `compute_prediction_core()` in `criteria.R` calls `get_adaptive_ridge()` on every
  fitness evaluation — this is expensive; cache the adaptive lambda per (P_train hash)
- For `pev_mean` / `cd_mean` with large populations, the GA calls `crossprod(P_train)` 
  repeatedly for the SAME training sets across generations; exploit this
- Consider using `RcppArmadillo` or at least LAPACK wrappers for the matrix inverses
  in `safe_matrix_inverse()`
- Parallel evaluation via `parallel::mclapply` already exists but is only triggered for
  >10 solutions and >1 core; tune the threshold

#### 2. Statistical correctness review
- **CD criterion sign**: `cd_mean()` returns `mean(leverage)` (positive), but the GA
  *minimizes* fitness. For training population selection the goal is to *maximize*
  coverage of the test set (higher leverage = better). Verify whether `cd_mean` should
  return `-mean(leverage)` so minimization gives optimal coverage. Check against the
  Rincent et al. (2012) and Isidro et al. (2015) papers.
- **PEV formula**: `compute_prediction_core()` computes
  `PEV = I + X_test (X'X)^{-1} X_test'`. The `I` term accounts for observation
  variance, which is correct for prediction of new observations but may not match
  the genomic-selection literature convention of just `X_test (X'X)^{-1} X_test'`.
  Provide both options (parameter `include_obs_variance = TRUE/FALSE`).
- **Mixed-model PEV**: validate `pev_mean_mm()` against `EMMREML` or `rrBLUP` output
  on the WheatData example; the Henderson (1984) formula involves G_22, G_21, V_11
  — check signs and matrix ordering carefully.
- **Distance criteria**: `distance_train_to_test_max` / `mean` minimize the
  *maximum/mean distance* from training to test. This is correct for spatial coverage
  but confusingly labeled — add documentation clarifying the direction.

#### 3. GA improvements to explore
- **Simulated annealing** as an alternative to pure GA — add `method = c("ga", "sa")` 
  option; SA is often faster for this problem class
- **Local search** post-GA: Fedorov exchange step is already optionally available;
  make it the default for a final polish pass
- **Adaptive population sizing**: reduce `npop` early when diversity collapses
- **Multi-start**: run several short GA runs and keep the best (addresses premature
  convergence better than restarts)
- **Check convergence window logic** in `check_multi_criteria_convergence()`:
  the circular buffer uses `%% convergence_window + 1` indexing but the buffer is
  pre-filled with zeros — early generations will have spurious zeros affecting the
  convergence signal; initialize with `NA` and filter

#### 4. Testing (TDD approach)
Use `testthat` (already in Suggests). Tests must be **fast** (< 30s total) and cover:

```
tests/testthat/
  test-criteria.R          # correctness of each criterion function
  test-criteria-signs.R    # verify GA goes in the right direction for each criterion
  test-ga-core.R           # subset_ga + subset_ga_single basic runs
  test-ga-convergence.R    # convergence detection, restart
  test-distance.R          # distance criteria with D vs dmds
  test-validation.R        # validate_matrix_params edge cases
  test-legacy.R            # GenAlgForSubsetSelection backward compat
  test-parallel.R          # mc.cores > 1 gives same result as sequential
```

Key design rule: use small n (30–50 individuals), small p (5 PCs), few iterations
(niterations=20, npop=20) in all tests. Never load WheatData in tests (too slow);
generate synthetic data with `set.seed()`.

Criterion correctness test template:
```r
test_that("a_optimality decreases with larger training set", {
  set.seed(1); P <- matrix(rnorm(50*5), 50, 5); rownames(P) <- paste0("ind", 1:50)
  test <- paste0("ind", 1:10); cands <- setdiff(rownames(P), test)
  a_small <- a_optimality(sample(cands, 10), test, P)
  a_large <- a_optimality(sample(cands, 25), test, P)
  expect_lt(a_large, a_small)
})
```

GA direction test template:
```r
test_that("subset_ga with cd_mean selects training set with high coverage", {
  set.seed(42); P <- matrix(rnorm(80*5), 80, 5); rownames(P) <- paste0("ind", 1:80)
  test <- paste0("ind", 1:20); cands <- setdiff(rownames(P), test)
  res <- subset_ga(P, cands, test, ntoselect=20, criterion="cd_mean",
                   niterations=30, npop=20, minitbefstop=15, verbose=FALSE)
  # Compare against a random subset
  random_cd <- cd_mean(sample(cands, 20), test, P)
  expect_lt(res$best_fitness, random_cd)  # GA should beat random (check sign!)
})
```

#### 5. Documentation
- Convert all Roxygen2 examples to use a shared tiny synthetic dataset (not WheatData)
- Write a short vignette `vignettes/quickstart.Rmd`: 10 minutes, one dataset,
  one run, one plot
- Write `vignettes/criteria-guide.Rmd`: table of all criteria, formula, use case,
  direction (minimize/maximize), recommended P input type (coordinates vs distance matrix)
- Remove or fix `vignettes/advanced_ga_features.Rmd` (currently broken)
- Add `pkgdown` config for GitHub Pages docs

#### 6. CRAN compliance fixes (from current R CMD check)
- Fix `WARNING: Undocumented code objects` for legacy aliases:
  `Amat.pieces`, `GenAlgForSubsetSelectionMO`, `GenAlgForSubsetSelectionMONoTest`,
  `calculate_population_distances`, `dist_to_test`, `dist_to_test2`, `disttoideal`,
  `evaluate_population_legacy`, `evaluate_population_optimized`, `neg_dist_in_train`,
  `neg_dist_in_train2`, `unified_distance_criterion`
  → Either add minimal `@rdname` docs or move them to non-exported internals
- Fix `WARNING: Documented arguments not in \usage: safe_matrix_inverse tol`
- Remove `vignettes/STPGA-package.log` (leftover)
- Ensure `inst/doc/` contains built vignettes before `R CMD build`

## Key files

| File | Role |
|------|------|
| `R/ga-core.R` | `subset_ga`, `subset_ga_single`, convergence, restart |
| `R/criteria.R` | All criterion functions + unified `criterion()` dispatcher |
| `R/distances.R` | Distance-based criteria |
| `R/evaluation.R` | `evaluate_population`, caching |
| `R/ga-operators.R` | Selection, crossover, mutation, diversity |
| `R/validation.R` | Input validation utilities |
| `R/performance.R` | Timing, benchmarking utilities |
| `R/matrix-operations.R` | Stable linear algebra (may overlap with criteria.R) |
| `tests/testthat/` | Test suite (expand significantly) |
| `NEWS.md` | Changelog (created in v7.0.1) |

## Before starting

1. `cd ~/Dropbox/MasterFolder/RESEARCH/ModernClinicalDesignbyOptimization/SourceCode/STPGA`
2. `git status` — review pending uncommitted changes (many in working tree from v7 dev)
3. `git log --oneline -5`
4. Run existing tests: `Rscript -e "devtools::test()"`
5. Run quick baseline check: `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --no-vignettes ../STPGA_7.0.1.tar.gz`

## Definition of done (v8.0.0 CRAN-ready)

- [ ] `R CMD check` passes with 0 ERRORs, 0 WARNINGs, ≤ 2 NOTEs
- [ ] All `testthat` tests pass in < 30s
- [ ] CD and PEV criteria confirmed statistically correct (with citations)
- [ ] GA direction verified for each criterion (test coverage)
- [ ] `NEWS.md` updated
- [ ] `DESCRIPTION` version = 8.0.0, date = today
- [ ] README updated with new API examples
- [ ] Vignette `quickstart.Rmd` builds cleanly
