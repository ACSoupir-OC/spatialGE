# spatialGE Design Decisions

**Document Purpose:** Record architectural and design decisions for posterity

**Last Updated:** 2026-03-18  
**Version:** 2.0.0

---

## Table of Contents

1. [Package Architecture](#package-architecture)
2. [Data Structures](#data-structures)
3. [Function Design Patterns](#function-design-patterns)
4. [Performance Optimization](#performance-optimization)
5. [Testing Strategy](#testing-strategy)
6. [Parallelization Approach](#parallelization-approach)
7. [Error Handling](#error-handling)
8. [Documentation Standards](#documentation-standards)

---

## Package Architecture

### Decision: Modular Ingestor Pattern (v2.0.0)

**Date:** 2026-02-26  
**Status:** Adopted

**Context:**
The original `STList_legacy()` function (50KB) handled all data input formats in a single monolithic function. This made the code difficult to maintain, test, and extend.

**Decision:**
Implement a modular ingestor pattern with platform-specific functions:
- `import_visium.R` - 10X Visium data
- `import_xenium.R` - 10X Xenium data
- `import_smi.R` - SMI data
- `ingest_generic.R` - Generic CSV/TSV matrices

**Consequences:**
- ✅ Easier to add new platforms
- ✅ Better testability (test each ingestor independently)
- ✅ Clearer error messages (platform-specific)
- ⚠️ Slight increase in total lines of code
- ⚠️ Need to maintain dispatch logic in `STlist()`

**Implementation:**
```r
# Main dispatcher
STlist <- function(data, platform = "auto", ...) {
  if (platform == "auto") {
    platform <- detect_platform(data)
  }
  
  switch(platform,
    visium = ingest_visium(data, ...),
    xenium = ingest_xenium(data, ...),
    generic = ingest_generic(data, ...),
    stop("Unknown platform: ", platform)
  )
}
```

**References:** `R/STlist.R`, `R/ingest_*.R`

---

### Decision: Preserve Legacy Functions

**Date:** 2026-02-26  
**Status:** Adopted

**Context:**
Users have existing workflows and publications based on `STList_legacy()`. Breaking changes would disrupt reproducibility.

**Decision:**
Keep `STList_legacy()` as a deprecated but functional alternative. Mark with `@deprecated` tag and issue warnings on use.

**Consequences:**
- ✅ Backward compatibility maintained
- ✅ Users can migrate at their own pace
- ⚠️ Two code paths to maintain
- ⚠️ Confusion for new users (which function to use?)

**Mitigation:**
- Clear documentation in README
- Deprecation warnings point to new function
- Vignettes use new `STlist()` function

---

## Data Structures

### Decision: S4 Class for STlist Objects

**Date:** 2022 (original package)  
**Status:** Retained

**Context:**
Spatial transcriptomics data has complex structure: multiple samples, count matrices, spatial coordinates, metadata, and analysis results.

**Decision:**
Use S4 class definition with slots:
- `data` - List of count matrices (one per sample)
- `gene_meta` - Gene-level metadata
- `cell_meta` - Cell/spot-level metadata
- `spatial` - Spatial coordinates
- `results` - Analysis results (SThet, STclust, etc.)

**Consequences:**
- ✅ Type safety (slots have defined classes)
- ✅ Clear structure for users
- ✅ Can define methods (show, plot, summary)
- ⚠️ S4 has steeper learning curve than S3
- ⚠️ More verbose than simple lists

**Rationale:**
The complexity of spatial transcriptomics data justifies S4's rigor. The class structure is documented in `R/classDefinitions.R`.

---

### Decision: Store Results in STlist Object

**Date:** 2022  
**Status:** Retained with modifications

**Context:**
Users need to track analysis results alongside original data. Separate objects lead to synchronization issues.

**Decision:**
Store analysis results in `results` slot as nested list:
```r
stlist@results$SThet$Moran
stlist@results$STclust$domains
stlist@results$STdiff$table
```

**Consequences:**
- ✅ Single object contains all analysis state
- ✅ Results stay with data
- ⚠️ Object can become very large
- ⚠️ Need clear naming conventions

**Mitigation:**
- Document result structure in function help
- Provide `summarize_STlist()` to inspect results
- Consider adding `save_results()` to export separately

---

## Function Design Patterns

### Decision: Consistent Return Values

**Date:** 2026-02-26  
**Status:** Adopted

**Context:**
Different functions returned different structures, making it hard to chain operations or write generic code.

**Decision:**
All analysis functions return the modified STlist object with results stored in `@results` slot. Additionally, they return results invisibly for immediate inspection.

**Pattern:**
```r
SThet <- function(stlist, ...) {
  # ... computation ...
  
  stlist@results$SThet <- results
  
  invisible(list(
    stlist = stlist,
    results = results
  ))
}
```

**Consequences:**
- ✅ Consistent API across functions
- ✅ Easy to chain: `stlist %>% SThet() %>% STclust()`
- ✅ Results accessible both ways

---

### Decision: Parameter Naming Conventions

**Date:** 2026-02-26  
**Status:** Adopted

**Context:**
Inconsistent parameter names across functions caused confusion.

**Decision:**
Standardize parameter names:
- `stlist` - STlist object (always first parameter)
- `assay` - Which assay to use ("counts", "scaled", etc.)
- `spots` / `cells` - Cell/spot identifiers
- `genes` - Gene names or indices
- `nsim` / `nperm` - Number of simulations/permutations
- `B` - Number of permutations (alternative to nsim)
- `verbose` - Logical, print progress
- `BPPARAM` - Parallelization parameters

**Consequences:**
- ✅ Easier to learn API
- ✅ Better autocomplete in IDEs
- ⚠️ Breaking change for some parameters

---

## Performance Optimization

### Decision: fastcluster Over DynamicTreeCut

**Date:** 2026-03-05  
**Status:** In Progress

**Context:**
Profiling showed DynamicTreeCut consumes 45% of STclust runtime. It's a pure R implementation without optimization.

**Decision:**
Replace with `fastcluster` package (C++ backend) for 20-35% speedup.

**Consequences:**
- ✅ Significant performance improvement
- ✅ Drop-in replacement (same API)
- ✅ No changes to output
- ⚠️ New dependency

**Implementation:**
```r
# Before
hc <- hclust(dist_matrix, method = "ward.D2")
clusters <- DynamicTreeCut::cutreeDynamic(hc, ...)

# After
library(fastcluster)
hc <- fastcluster::hclust(dist_matrix, method = "ward.D2")
clusters <- DynamicTreeCut::cutreeDynamic(hc, ...)  # Still need DTC for cutting
```

**Note:** Only `hclust()` is replaced; tree cutting still uses DynamicTreeCut.

**References:** `OPTIMIZATION_SUMMARY.md`, `PERFORMANCE_ANALYSIS.md`

---

### Decision: Lazy Distance Matrix Computation

**Date:** 2026-03-05  
**Status:** Planned

**Context:**
Distance matrix calculation consumes 35% of STclust runtime. For large datasets (>10K cells), this becomes prohibitive.

**Decision:**
Implement lazy computation:
1. Pre-compute only when needed
2. Cache results for reuse across parameter sweeps
3. Use HDF5Array for out-of-memory storage

**Consequences:**
- ✅ Faster for parameter sweeps
- ✅ Enables larger datasets
- ⚠️ More complex implementation
- ⚠️ Memory management overhead

---

### Decision: Avoid gene-by-gene Loops

**Date:** 2026-02-26  
**Status:** Adopted where possible

**Context:**
R loops are slow. Gene-by-gene operations (e.g., Moran's I for each gene) can be vectorized.

**Decision:**
Use matrix operations instead of loops:
```r
# Slow
for (g in genes) {
  moran[g] <- calculate_moran(gene_expr[, g])
}

# Fast (where mathematically possible)
moran <- calculate_moran_vectorized(gene_expr, weights)
```

**Consequences:**
- ✅ 10-100x speedup for vectorizable operations
- ⚠️ Not all operations can be vectorized
- ⚠️ Higher memory usage for intermediate matrices

---

## Testing Strategy

### Decision: testthat Framework

**Date:** 2026-02-27  
**Status:** Adopted

**Context:**
Package needs comprehensive test coverage for CRAN submission and regression prevention.

**Decision:**
Use `testthat` (>= 3.0.0) with:
- Unit tests for individual functions
- Integration tests for workflows
- Regression tests with saved reference data
- Helper functions for common test data

**Structure:**
```
tests/testthat/
├── setup.R              # Global setup
├── helper-data.R        # Test data generators
├── helper-regression.R  # Regression test helpers
├── test-stlist.R        # STlist tests
├── test-sthet.R         # SThet tests
├── test-stclust.R       # STclust tests
├── test-stdiff.R        # STdiff tests
└── ...
```

**Coverage Goal:** >75% by end of Q2 2026

---

### Decision: Regression Testing with Saved Data

**Date:** 2026-02-27  
**Status:** Adopted

**Context:**
Need to ensure refactoring doesn't change results.

**Decision:**
Save reference outputs for key functions and compare:
```r
test_that("SThet produces consistent results", {
  expected <- readRDS("reference/sthet_result.rds")
  actual <- SThet(test_stlist)
  expect_equal(actual@results$SThet, expected, tolerance = 1e-10)
})
```

**Consequences:**
- ✅ Catches unintended changes
- ✅ Documents expected behavior
- ⚠️ Reference files need updating when behavior intentionally changes

---

## Parallelization Approach

### Decision: mclapply as Default, foreach as Fallback

**Date:** 2026-02-26  
**Status:** Adopted with known limitation

**Context:**
Spatial analyses are computationally intensive. Parallelization is essential for large datasets.

**Decision:**
Use `parallel::mclapply()` as default (fastest on Unix), with `foreach` + `doParallel` as Windows-compatible fallback.

**Pattern:**
```r
if (.Platform$OS.type == "unix") {
  results <- mclapply(genes, process_gene, mc.cores = BPPARAM$workers)
} else {
  registerDoParallel(BPPARAM$workers)
  results <- foreach(g = genes) %dopar% process_gene(g)
}
```

**Consequences:**
- ✅ Best performance on Unix/Linux/Mac
- ✅ Works on Windows (slower)
- ⚠️ Code duplication
- ⚠️ Need to test both paths

**Future:** Consider `BiocParallel` for unified interface.

---

### Decision: User-Controlled Parallelization

**Date:** 2026-02-26  
**Status:** Adopted

**Context:**
Users may want to control parallelization (e.g., on shared systems).

**Decision:**
All parallel functions accept `BPPARAM` parameter:
```r
STdiff(stlist, BPPARAM = BiocParallel::SnowParam(workers = 4))
```

Default: `BiocParallel::bpparam()` (respects system settings)

**Consequences:**
- ✅ Flexible for different environments
- ✅ Users can tune performance
- ⚠️ More complex function signatures

---

## Error Handling

### Decision: Consistent Error Messages

**Date:** 2026-02-26  
**Status:** Adopted

**Context:**
Inconsistent error messages made debugging difficult.

**Decision:**
Use `cli` package for styled error messages:
```r
if (!inherits(stlist, "STlist")) {
  cli::cli_abort(c(
    "x" = "Input must be an {.cls STlist} object",
    "i" = "Use {.fn STlist} to create one"
  ))
}
```

**Consequences:**
- ✅ Clear, actionable errors
- ✅ Consistent formatting
- ✅ Helpful suggestions
- ⚠️ New dependency (cli)

---

### Decision: Input Validation at Function Entry

**Date:** 2026-02-26  
**Status:** Adopted

**Context:**
Better to fail fast with clear errors than crash deep in computation.

**Decision:**
Validate all inputs at function start:
1. Check class/type of required arguments
2. Check for missing required parameters
3. Check for incompatible parameter combinations
4. Check data dimensions and structure

**Pattern:**
```r
STclust <- function(stlist, assay = "counts", ws = 0.025, ...) {
  # Validation
  check_stlist(stlist)
  check_assay(stlist, assay)
  check_scalar_numeric(ws, min = 0, max = 1)
  
  # ... computation ...
}
```

---

## Documentation Standards

### Decision: Roxygen2 for All Functions

**Date:** 2022  
**Status:** Retained

**Context:**
CRAN requires comprehensive documentation.

**Decision:**
All exported functions must have:
- `@title` - One-line description
- `@description` - Detailed explanation
- `@param` - All parameters documented
- `@return` - Return value description
- `@examples` - Working examples
- `@export` - Export tag

**Consequences:**
- ✅ CRAN-compliant
- ✅ Auto-generates manual pages
- ✅ IDE autocomplete support
- ⚠️ Documentation must be maintained

---

### Decision: Vignettes for Workflows

**Date:** 2022  
**Status:** Retained

**Context:**
Function documentation is insufficient for complex workflows.

**Decision:**
Provide vignettes for:
1. Getting started (basic workflow)
2. Data import (all supported formats)
3. Spatial autocorrelation (SThet)
4. Tissue domain detection (STclust)
5. Differential expression (STdiff)
6. Gene set enrichment (STenrich)
7. Expression gradients (STgradient)

**Consequences:**
- ✅ Users can learn by example
- ✅ Demonstrates best practices
- ⚠️ Vignettes must be maintained
- ⚠️ Increases build time

---

## References

### Related Documents
- `RefactoringPlan.md` - Detailed refactoring strategy
- `OPTIMIZATION_SUMMARY.md` - Performance analysis
- `PERFORMANCE_ANALYSIS.md` - Profiling data
- `DEVELOPMENT.md` - Development guidelines
- `AGENTS/TASKS.md` - Current task list

### External References
- Ospina et al. (2022). spatialGE. Bioinformatics. https://doi.org/10.1093/bioinformatics/btac145
- Soupir et al. (2025). Permutation-Free Nearest Neighbor G Function. bioRxiv. https://doi.org/10.1101/2025.06.11.659088

---

**Document Maintained By:** Alex Soupir  
**Contact:** Alex.Soupir@moffitt.org

---

*This document should be updated when significant architectural decisions are made.*
