# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

STPGA (Selection of Training Populations by Genetic Algorithm) is an R package for optimizing training set selection in high-dimensional prediction problems, particularly useful in breeding value prediction and experimental design.

## Development Commands

### Build and Check
```bash
# Build the package
R CMD build .

# Check the package (runs tests, examples, and validates structure)
R CMD check STPGA_*.tar.gz

# Install the package locally
R CMD INSTALL .

# Load in R for interactive development
R -e "devtools::load_all()"
```

### Documentation
```bash
# Generate documentation from roxygen comments (if using roxygen2)
R -e "devtools::document()"

# Build PDF manual
R CMD Rd2pdf .
```

### Testing
Currently no formal test suite exists. To test functions interactively:
```r
# Load package
library(STPGA)

# Load example data
data(WheatData)

# Run example from documentation
example(GenAlgForSubsetSelection)
```

## Architecture

### Core Algorithm Structure
The package implements genetic algorithms for subset selection with two main variants:
1. **Single-objective optimization** (`GenAlgForSubsetSelection.R`, `GenAlgForSubsetSelectionNoTest.R`)
2. **Multi-objective optimization** (`GenAlgForSubsetSelectionMO.R`, `GenAlgForSubsetSelectionMONoTest.R`)

### Optimization Criteria
Multiple criteria functions in `/R/` implement different optimality measures:
- **Classic design criteria**: A-optimality (`AOPT.R`), D-optimality (`DOPT.R`), E-optimality (`EOPT.R`)
- **CD-based criteria**: Maximum (`CDMAX*.R`) and mean (`CDMEAN*.R`) coefficient of determination
- **PEV-based criteria**: Maximum (`PEVMAX*.R`) and mean (`PEVMEAN*.R`) prediction error variance

The numeric suffixes (0, 2) indicate different computational approaches or approximations.

### Key Design Patterns
1. **Criteria functions** follow a consistent interface: they accept design matrices and return scalar fitness values
2. **Genetic operations** are modularized in separate files (`GenerateCrossesfromElites.R`, `makeonecross.R`)
3. **Distance calculations** are centralized in utility functions (`DistanceBasedCriteria.R`, `disttoideal.R`)

### Integration Points
- Depends on `AlgDesign` for Fedorov exchange algorithm initialization
- Uses `emoa` for multi-objective optimization utilities
- Optional integration with `EMMREML` for mixed model calculations
- Can leverage `glmnet` for regularized regression in some criteria

## Important Considerations
- The package assumes marker/feature data is available for all individuals
- Designed for scenarios where p < n (features less than samples) for optimal performance
- Uses ridge regression internally when needed for numerical stability
- Multi-objective functions use Pareto optimization concepts