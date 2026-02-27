# spatialGE Refactoring Plan

## Overview

This document provides a comprehensive review of the `spatialGE` R package codebase. It documents each function, its purpose, current implementation, and potential refactoring opportunities. The goal is to improve code maintainability while preserving all existing functionality.

## Package Structure

The package consists of:
- **Core functions**: Main user-facing functions
- **Helper functions**: Internal functions supporting core functionality  
- **Class definitions**: S4 class for STlist objects
- **Utility functions**: General-purpose helpers
- **Plotting functions**: Visualization helpers

---

## Core Functions

### 1. STlist / STList_legacy (STlist.R)

**Purpose**: Creates STlist objects from spatial transcriptomics data sources.

**Key Responsibilities**:
- Ingest data from multiple sources (Visium, Xenium, CosMx, Seurat, generic matrices)
- Match count and coordinate data
- Create S4 STlist objects with proper structure

**Current Issues**:
- STList_legacy is large (50KB) and mixes multiple data source handling
- New modular STlist implementation exists but legacy still documented
- Inconsistent error handling across different input types

**Refactoring Suggestions**:
- Consolidate data source detection into separate, smaller files
- Create dedicated modules for each platform (visium, xenium, cosmx, seurat)
- Standardize error messaging and codes
- Remove or deprecate STList_legacy in favor of STlist

---

### 2. SThet (SThet.R)

**Purpose**: Computes global spatial autocorrelation (Moran's I, Geary's C).

**Key Responsibilities**:
- Calculate spatial weights using k-nearest neighbors or Euclidean distances
- Compute Moran's I and Geary's C statistics
- Store results in STlist object

**Current Issues**:
- Two helper functions (gene_moran_i_notest, gene_geary_c_notest) are verbose (~400 lines combined)
- Parallelization uses mclapply (Unix-only) with commented-out alternatives
- Hard-coded column names in gene_meta access
- Complex nested loops for storing results

**Refactoring Suggestions**:
- Extract weight calculation to separate function/module
- Simplify result storage with vectorized operations
- Add Windows-compatible parallelization (foreach + doParallel)
- Use helper functions for gene_meta access

---

### 3. STclust (STclust.R)

**Purpose**: Spatially-informed hierarchical clustering for tissue domain detection.

**Key Responsibilities**:
- Calculate variable genes using Seurat's FindVariableFeatures
- Compute scaled expression and spatial distance matrices
- Apply weighted combination of distances
- Perform hierarchical clustering with DynamicTreeCut or fixed k

**Current Issues**:
- Weight calculation (calculate_weighted_dist) is verbose
- Multiple clustering methods (dtc vs ks) handled in same function
- Hard-coded default weights (0.025)
- Results stored with complex column naming scheme

**Refactoring Suggestions**:
- Separate distance matrix calculation from clustering logic
- Create dedicated functions for DTC and fixed-k clustering
- Parameterize weight defaults with better documentation
- Simplify result column naming conventions

---

### 4. STdiff (STdiff.R)

**Purpose**: Differential gene expression analysis with spatial mixed models.

**Key Responsibilities**:
- Run non-spatial tests (mixed models, t-tests, Wilcoxon)
- Select top DE genes based on adjusted p-values
- Fit spatial models using spaMM for selected genes
- Compile and format results

**Current Issues**:
- Extremely large function (42KB, ~900 lines)
- Complex annotation handling with coded annotations
- Multiple test types mixed in one function
- Heavy parallelization with mclapply
- Sparse matrix handling scattered throughout

**Refactoring Suggestions**:
- **High priority**: Split into smaller functions:
  - `STdiff_run_nonspatial` - Non-spatial tests
  - `STdiff_run_spatial` - Spatial mixed models
  - `STdiff_compile_results` - Result formatting
- Separate test type logic (mm, t_test, wilcoxon)
- Simplify annotation recoding
- Create helper for sparse matrix expansion

---

### 5. STenrich (STenrich.R)

**Purpose**: Randomization test for gene set spatial enrichment.

**Key Responsibilities**:
- Calculate gene set expression (average or GSVA scores)
- Identify high-expression spots above threshold
- Compute sum of pairwise distances
- Perform permutation testing

**Current Issues**:
- GSVA calculation path vs average expression path mixed
- Distance matrix computation commented out (HDF5Array)
- Complex permutation loop
- Memory management unclear

**Refactoring Suggestions**:
- Separate GSVA calculation to dedicated module
- Implement out-of-memory distance computation if needed
- Simplify permutation testing with clearer structure
- Add better memory profiling

---

### 6. STgradient (STgradient.R)

**Purpose**: Detect gene expression spatial gradients relative to reference domain.

**Key Responsibilities**:
- Calculate distances from non-reference spots to reference domain
- Compute distance summaries (min or average)
- Fit linear models and Spearman correlations
- Handle outlier removal and robust regression

**Current Issues**:
- Large function (~21KB, ~500 lines)
- Multiple data processing paths (robust, outlier removal, log transform)
- Complex distance calculation logic
- Hard-coded neighborhood tolerance values

**Refactoring Suggestions**:
- Split into:
  - `STgradient_distance_calc` - Distance computation
  - `STgradient_correlation` - Model fitting
  - `STgradient_handle_outliers` - Outlier management
- Parameterize neighborhood tolerance
- Simplify log-transform option handling

---

### 7. STplot (STplot.R)

**Purpose**: Plot gene expression, clusters, or metadata in spatial context.

**Key Responsibilities**:
- Dispatch to appropriate plotting function based on input
- Handle Visium y-axis flipping
- Manage plot aesthetics

**Current Issues**:
- Small wrapper function (reasonable size)
- Dispatches to other plotting functions

**Refactoring Suggestions**:
- Good modular design, minimal changes needed
- Consider adding more plot types

---

### 8. STplot_interpolation (STplot_interpolation.R)

**Purpose**: Interpolate gene expression onto continuous spatial surface.

**Key Responsibilities**:
- Perform kriging interpolation on gene expression
- Generate interpolated plots

**Current Issues**:
- Depends on gene_interpolation.R

**Refactoring Suggestions**:
- Review kriging implementation for efficiency

---

## Class Definitions

### STlist Class (classDefinitions.R)

**Purpose**: S4 class defining STlist object structure.

**Slots**:
- `counts`: List of sparse count matrices per sample
- `spatial_meta`: List of spot/cell coordinates per sample
- `gene_meta`: List of gene-level statistics per sample
- `sample_meta`: Data frame with sample-level metadata
- `tr_counts`: List of transformed counts
- `gene_krige`: List of kriging results
- `misc`: Parameters and images

**Current Issues**:
- Basic class definition with show/summary methods
- Methods duplicated between show and summary

**Refactoring Suggestions**:
- Remove duplicate summary method (same as show)
- Add methods for accessing/modifying slots
- Document slot purposes more clearly

---

## Helper Functions

### compare_SThet (compare_SThet.R)

**Purpose**: Compare spatial autocorrelation statistics across samples.

**Key Responsibilities**:
- Extract Moran's I/Geary's C from STlist
- Create comparison plots colored by metadata

**Current Issues**:
- Relatively small function (~7KB)
- Uses ggplot2 with complex facet wrapping

**Refactoring Suggestions**:
- Extract color palette logic to separate function
- Simplify metadata handling

---

## Utility Functions

### detect_input (detect_input.R)

**Purpose**: Detect data source type from user inputs.

**Current Issues**:
- Large function handling multiple input types
- Complex nested conditionals

**Refactoring Suggestions**:
- Split detection by input type
- Create dispatch table for different sources

---

### filter_data (filter_data.R)

**Purpose**: Filter STlist data by various criteria.

**Current Issues**:
- Gene filtering logic
- Sample filtering logic

**Refactoring Suggestions**:
- Separate gene and sample filtering
- Add validation functions

---

### transform_data (transform_data.R)

**Purpose**: Normalize and transform count data.

**Current Issues**:
- VST calculation
- Log transformation
- Scaling

**Refactoring Suggestions**:
- Separate transformation steps
- Add option selection interface

---

### utils.R

**Purpose**: General utility functions.

**Key Functions**:
- `count_cores` - Detect available cores
- `color_parse` - Parse color palettes
- `raise_err` - Standardized error handling
- `expandSparse` / `makeSparse` - Sparse matrix conversion

**Refactoring Suggestions**:
- Good centralized utility location
- Add more comprehensive documentation

---

## Input Processing Functions

### ingest_generic, ingest_visium, ingest_xenium, ingest_smi

**Purpose**: Platform-specific data ingestion.

**Current Issues**:
- Each platform has separate file handling
- Some functions use HDF5 reading, others use MEX
- Complex file path discovery

**Refactoring Suggestions**:
- Create platform-specific modules
- Standardize file discovery logic
- Add input validation

---

## Plotting Helper Functions

### plot_counts, plot_image, plot_spatial_expression, plot_spatial_geneset, plot_spatial_meta

**Purpose**: Specialized plotting functions.

**Current Issues**:
- Multiple similar functions with overlapping code
- Color handling inconsistent across functions
- Some use ggplot2 directly, some use wrappers

**Refactoring Suggestions**:
- Consolidate common plotting logic
- Create base plot function with parameters
- Standardize color handling

---

## Code Quality Issues

### 1. Code Duplication

**Locations**:
- `gene_moran_i_notest` and `gene_geary_c_notest` share ~80% code
- `STdiff` has similar model fitting patterns
- Multiple plotting functions repeat ggplot2 setup

**Action**: Extract shared logic to helper functions.

### 2. Hard-coded Values

**Locations**:
- Default spatial weight: 0.025
- Default k values and deepSplit
- Neighborhood tolerance: c(0.75, 1.25) for Visium

**Action**: Move to configuration constants at top of files.

### 3. Comments and Documentation

**Locations**:
- Many functions have commented-out code blocks
- Parallelization alternatives commented out
- Some functions lack roxygen2 documentation

**Action**: Remove commented code, document all functions properly.

### 4. Error Handling

**Locations**:
- `raise_err` function exists but not consistently used
- Error codes (error0001, error0002, etc.) not documented
- Some functions use stop() directly

**Action**: Create error documentation, use raise_err consistently.

### 5. Memory Management

**Locations**:
- `rm()` calls scattered throughout
- No explicit memory profiling
- Large matrices not always handled efficiently

**Action**: Add memory profiling, document memory-intensive operations.

### 6. Parallelization

**Locations**:
- `parallel::mclapply` used extensively (Unix only)
- Windows alternatives commented out
- No fallback for parallelization failures

**Action**: Add Windows-compatible parallelization, error handling for parallel failures.

---

## Recommended Refactoring Order

### Phase 1: Foundation (High Priority)
1. **consolidate error handling** - Document and standardize error codes
2. **extract utility functions** - Move common helpers to utils.R
3. **fix parallelization** - Add Windows support, error handling

### Phase 2: Core Functions (Medium Priority)
4. **split STdiff** - Break into smaller focused functions
5. **simplify SThet** - Reduce complexity of Moran/Geary calculations
6. **clean up STclust** - Separate distance calculation from clustering

### Phase 3: Input/Output (Medium Priority)
7. **modularize input processing** - Platform-specific modules
8. **consolidate plotting** - Base plot function with parameters

### Phase 4: Documentation (Ongoing)
9. **add roxygen2 documentation** - All public functions
10. **add unit tests** - Test each refactored component

---

## Testing Strategy

After refactoring:
1. Verify all existing tests pass
2. Add tests for new helper functions
3. Test on sample datasets (melanoma example)
4. Validate output matches pre-refactoring results

---

## Notes

- This refactoring preserves all existing functionality
- Changes are additive and backward compatible where possible
- Legacy functions (STList_legacy) should be maintained for reproducibility
- Consider versioning changes that affect API

---

*Last updated: 2026-02-26*
*Reviewed by: Claw (spatialGE refactoring assistant)*