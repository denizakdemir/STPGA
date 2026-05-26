# STPGA News

## STPGA 7.0.1 (2026-05-26)

### Bug fixes

* **`subset_ga()` / `subset_ga_single()` — "argument Test is missing" error fixed.**
  `Test` now defaults to `NULL` in `subset_ga()`. The root cause was an R-language
  subtlety: assigning `NULL` to a named list element removes the key, so
  `subset_ga_single()` was silently dropping `Test` before forwarding the call via
  `do.call()`. Giving `Test = NULL` a default in `subset_ga()` makes both the
  no-test wrapper and any direct call without a test set work correctly.

* **Distance criteria (`max_to_test`, `mean_to_test`, `neg_min_internal`,
  `neg_mean_internal`) caused a "missing value where TRUE/FALSE needed" crash.**
  These criteria were implemented in `distances.R` but were not registered in the
  unified `criterion()` dispatch function in `criteria.R`. Every call with a
  distance criterion was silently caught by `tryCatch()` and returned `Inf` for all
  fitness values, making `improvement_rate = |Inf − Inf| / Inf = NaN`, and
  `NaN < tolconv` evaluates to `NA` in R, triggering the error. The four distance
  criteria are now wired into the switch in `criterion()`.

### Notes for users upgrading from 5.x

* The main entry-points are now `subset_ga()` (with test set) and
  `subset_ga_single()` (whole-population, no test set).  The legacy functions
  `GenAlgForSubsetSelection()` and `GenAlgForSubsetSelectionNoTest()` are still
  available for backward compatibility.

* Criterion names have changed from ALL-CAPS (`CDMEAN2`, `PEVMEAN2`, …) to
  `snake_case` (`cd_mean`, `pev_mean`, …).  All old names are mapped via the
  unified `criterion()` function, so existing code continues to work.

* **Distance-based criteria (`max_to_test`, `mean_to_test`) require a square
  distance matrix as `P`** (e.g., `D` in the demo data), not an MDS coordinate
  matrix (e.g., `dmds`).  MDS / PCA coordinate matrices should be used with
  prediction-based criteria (`cd_mean`, `pev_mean`, etc.).

---

## STPGA 7.0.0 (2025-07-27)

### Major rewrite

* New clean API: `subset_ga()`, `subset_ga_single()`, multi-objective variants.
* Unified `criterion()` dispatcher covering classical (A, D, E, G, I optimality),
  PEV-based, CD-based, mixed-model, and distance-based criteria with full legacy
  name support.
* Multi-criteria convergence detection with configurable window and restart
  mechanism.
* Rank-based and hybrid selection, adaptive mutation, diversity preservation
  (crowding + fitness sharing).
* Numerically stable matrix operations throughout (`safe_matrix_inverse()`,
  adaptive ridge).
* Literature-corrected PEV and CD formulas based on Henderson (1984) and
  Fedorov (1972).

---

## STPGA 5.2.1

* Previous CRAN release.
