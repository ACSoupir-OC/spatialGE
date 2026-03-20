# STclust Refactoring Plan

## Overview
Refactor `STclust` function into modular components while maintaining backward compatibility through a `STclust_legacy()` function.

**Status:** ✅ SThet refactored (commit 3fdd76e)
**Target:** STclust refactoring (next step)

---

## Current State Analysis

### File: `R/STclust.R` (350 lines)

**Main function:** `STclust()` (lines 1-208)
- Calculates spatially-informed clusters
- Uses weighted distance matrices (expression + spatial)
- Supports two clustering methods:
  - DTC (DynamicTreeCut) for adaptive k
  - Fixed k values for range testing

**Helper functions (lines 210-350):**
1. `calculate_dist_matrices()` - Calculate and scale distance matrices (35 lines)
2. `calculate_weighted_dist()` - Create weighted distance matrices (24 lines)
3. `get_hier_clusters_dtc()` - DTC-based clustering (45 lines)
4. `get_hier_clusters_ks()` - Fixed-k clustering (52 lines)

**Total:** ~156 lines of code (194 lines including documentation and whitespace)

---

## Refactoring Goals

### 1. Modular Structure
Split into separate files with clear separation of concerns:

```
R/STclust.R              # Thin wrapper + public API
R/STclust_core.R         # Main modular functions
R/STclust_helpers.R      # Helper/utility functions
R/STclust_legacy.R       # Original implementation (renamed)
```

### 2. Modular Functions

#### R/STclust_core.R (New)
- **`STclust_select_genes()`** - Variable gene selection and count filtering
  - Input: STlist, samples, topgenes
  - Output: Filtered expression matrices per sample
  
- **`STclust_calculate_distances()`** - Distance matrix calculation
  - Input: expression matrix, coordinates, dist_metric
  - Output: list of scaled distance matrices (expression + spatial)
  
- **`STclust_weight_distances()`** - Weighted distance combination
  - Input: scaled distance matrices, ws vector
  - Output: list of weighted distance matrices
  
- **`STclust_hierarchical()`** - Core clustering logic
  - Input: weighted distances, ws, ks, linkage
  - Output: list of cluster assignments per sample/weight
  
- **`STclust()`** - Public wrapper maintaining original API

#### R/STclust_helpers.R (New)
- **`calculate_dist_matrices()`** - Current helper (35 lines)
- **`calculate_weighted_dist()`** - Current helper (24 lines)
- **`get_hier_clusters_dtc()`** - DTC clustering (45 lines)
- **`get_hier_clusters_ks()`** - Fixed-k clustering (52 lines)

### 3. Legacy Preservation

Rename current `STclust()` to `STclust_legacy()` with:
- Same signature and behavior
- Exported for reproducibility
- Users can explicitly call `STclust_legacy()` for old behavior

---

## Implementation Steps

### Step 1: Baseline Testing
- Run existing test suite to establish baseline
- Document current test results (PASS/FAIL/WARN counts)

### Step 2: Create Helper Functions File
- Extract all helper functions to `R/STclust_helpers.R`
- Keep original function names for backward compatibility
- Add documentation to each function

### Step 3: Create Core Modular Functions
Create `R/STclust_core.R` with:

#### Function 1: `STclust_select_genes()`
```r
STclust_select_genes = function(x, samples, topgenes, cores){
  # 1. Calculate VST for variable genes
  # 2. Subset to topgenes per sample
  # 3. Return filtered expression matrices
}
```

#### Function 2: `STclust_calculate_distances()`
```r
STclust_calculate_distances = function(expr_dat, coord_dat, dist_metric){
  # Use existing calculate_dist_matrices logic
  # Return scaled distance matrices
}
```

#### Function 3: `STclust_weight_distances()`
```r
STclust_weight_distances = function(scaled_dists, ws){
  # Use existing calculate_weighted_dist logic
  # Return weighted distance matrices
}
```

#### Function 4: `STclust_hierarchical()`
```r
STclust_hierarchical = function(weighted_dists, ws, ks, linkage, deepSplit, verbose){
  # Branch based on ks parameter:
  # - If 'dtc': call get_hier_clusters_dtc
  # - If numeric: call get_hier_clusters_ks
  # Return cluster assignments
}
```

#### Function 5: `STclust()` (Public Wrapper)
```r
STclust = function(x=NULL, samples=NULL, ws=0.025, dist_metric='euclidean',
                   linkage='ward.D2', ks='dtc', topgenes=2000, deepSplit=FALSE,
                   cores=NULL, verbose=TRUE){
  
  # Call modular functions in sequence:
  # 1. x <- STclust_select_genes(x, samples, topgenes, cores)
  # 2. res_ls <- parallel::mclapply(samples, function(i){
  #      scaled_dists <- STclust_calculate_distances(...)
  #      weighted_dists <- STclust_weight_distances(scaled_dists, ws)
  #      hierclusters <- STclust_hierarchical(weighted_dists, ws, ks, linkage, deepSplit, verbose)
  #    })
  # 3. Update STlist with results
  # 4. Return x
}
```

### Step 4: Create Legacy Version
- Rename current `STclust` to `STclust_legacy` in new file
- Keep original implementation unchanged
- Add deprecation note in documentation

### Step 5: Update Main STclust.R File
```r
##
# STclust: Spatial clustering (modular implementation)
#
# @importFrom magrittr %>%
# @importFrom methods as is new

source('R/STclust_helpers.R')
source('R/STclust_core.R')

# Legacy function is defined in STclust_legacy.R
if(file.exists('R/STclust_legacy.R')){
  source('R/STclust_legacy.R')
}
```

### Step 6: Testing Strategy

#### Existing Tests to Verify:
1. `tests/testthat/test-analysis.R` - STclust test (already exists)
2. Tests in test-differential.R that use STclust
3. Tests in test-enrichment.R that use STclust

#### New Tests to Add:
1. **`test-stclust.R`** - Dedicated STclust tests
   - Basic clustering with ks=2:3, ws=0.025
   - DTC clustering (adaptive k)
   - Fixed k clustering
   - Multiple samples
   - Edge cases (single sample, single gene, etc.)

2. **Comparison test: STclust() vs STclust_legacy()**
   - Verify identical results
   - Test with various parameters

3. **Module-level tests:**
   - `STclust_select_genes()` - gene selection correctness
   - `STclust_calculate_distances()` - distance matrix accuracy
   - `STclust_weight_distances()` - weight calculation
   - `STclust_hierarchical()` - clustering logic

### Step 7: Documentation Updates

#### README.md
Add section about STclust refactoring:
```markdown
## STclust Refactoring

The `STclust` function has been refactored into modular components for better maintainability:

- **`STclust()`** - New modular implementation (default)
- **`STclust_legacy()`** - Original implementation (for reproducibility)

### Usage

The new `STclust()` function maintains the same API as the original:

```r
melanoma <- STclust(melanoma, ks=2:3, ws=0.025)
```

If you need the original behavior for reproducibility, use `STclust_legacy()`:

```r
melanoma <- STclust_legacy(melanoma, ks=2:3, ws=0.025)
```

Both functions produce identical results.
```

#### roxygen2 Documentation
- Update `STclust()` documentation with refactoring notes
- Add deprecation note for `STclust_legacy()` (optional, or keep as "legacy for reproducibility")
- Document each new modular function

### Step 8: Commit and Push Strategy

**Commit 1: Baseline & Structure**
- Run baseline tests
- Create helper functions file
- Verify baseline still passes

**Commit 2: Core Modular Functions**
- Add core modular functions
- Test each module independently

**Commit 3: Wrapper & Legacy**
- Create STclust wrapper
- Create STclust_legacy
- Add comparison test

**Commit 4: Documentation**
- Update README.md
- Update roxygen2 docs
- Final test suite run

**Total Commits:** 4 commits (or combine as appropriate)

---

## Testing Checklist

### Before Refactoring:
- [ ] Run full test suite, record baseline
- [ ] Document test results (PASS/FAIL/WARN/SKIP counts)

### During Refactoring:
- [ ] Each helper function testable in isolation
- [ ] Each core function testable in isolation
- [ ] Wrapper function produces same output as legacy

### After Refactoring:
- [ ] Full test suite passes
- [ ] Comparison test: STclust() == STclust_legacy()
- [ ] No regressions in dependent tests (SThet, STdiff, STenrich)

### Test Coverage:
- [x] Basic clustering (ks=2:3, ws=0.025)
- [ ] DTC clustering (adaptive k)
- [ ] Fixed k clustering
- [ ] Multiple samples
- [ ] Edge cases (single sample, boundary conditions)

---

## Complexity Assessment

**Current STclust:** 350 lines, 4 helper functions

**Estimated new structure:**
- `STclust.R`: ~50 lines (wrapper + sourcing)
- `STclust_core.R`: ~200 lines (modular functions)
- `STclust_helpers.R`: ~150 lines (helper functions)
- `STclust_legacy.R`: ~350 lines (original)

**Total lines:** ~750 lines (increase due to modularization and documentation)

**Benefit:** Better maintainability, testability, and understandability

---

## Risk Assessment

**Low Risk:**
- Helper functions are self-contained and well-tested by existing tests
- Modular functions can be tested independently
- Legacy function preserves original behavior

**Medium Risk:**
- Wrapper function needs careful parameter passing
- Parallelization logic needs verification

**Mitigation:**
- Extensive testing at each step
- Comparison tests with legacy implementation
- Incremental commits for easy rollback

---

## Timeline Estimate

| Step | Estimated Time |
|------|----------------|
| Baseline testing | 15 min |
| Helper functions file | 20 min |
| Core modular functions | 45 min |
| Wrapper & legacy | 15 min |
| Testing & validation | 30 min |
| Documentation | 15 min |
| Commits & push | 10 min |
| **Total** | **~2.5 hours** |

---

## Next Steps

1. ✅ Review and approve this plan
2. Begin implementation
3. Follow incremental approach with check-ins
4. Test thoroughly at each step
5. Document and commit when ready

**Ready to proceed?** Let me know if you'd like to start implementation or modify the plan.