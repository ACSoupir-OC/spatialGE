# spatialGE - Refactoring Tasks

**Version**: 1.0 (In Progress - Modular Refactoring)  
**Created**: 2026-03-19 13:35 UTC  
**Last Updated**: 2026-03-19 13:35 UTC  
**Status**: 🏃 STclust & STdiff Complete, STenrich/STgradient/SThet Pending  
**Package Status**: Ready for further development

---

## 🎯 Current State

### ✅ Completed Functions

**1. STclust** - Spatial clustering
- ✅ `R/STclust.R` - Public interface
- ✅ `R/STclust_core.R` - Core implementation (modular)
- ✅ `R/STclust_helpers.R` - Helper functions
- ✅ `R/STclust_legacy.R` - Legacy (identical to original)
- ✅ Tests: 30/30 passing

**2. STdiff** - Spatial differential analysis
- ✅ `R/STdiff.R` - Public interface
- ✅ `R/STdiff_core.R` - Core implementation (modular)
- ✅ `R/STdiff_helpers.R` - Helper functions
- ✅ `R/STdiff_nonspatial.R` - Non-spatial tests
- ✅ `R/STdiff_spatial.R` - Spatial tests
- ✅ `R/STdiff_results.R` - Results handling
- ✅ `R/STdiff_volcano.R` - Volcano plots
- ✅ `R/STdiff_legacy.R` - Legacy (backward compatible)

### 🏃 In Progress / Pending

**3. STenrich** - Spatial enrichment
- 📄 `R/STenrich.R` - Public interface (complete)
- 📄 `R/STenrich_helpers.R` - Helper functions (complete)
- 📄 `R/STenrich_legacy.R` - Legacy implementation (complete)
- ❌ **Missing**: `R/STenrich_core.R` - Core implementation
- ❌ **Missing**: Test suite

**4. STgradient** - Spatial gradient analysis
- 📄 `R/STgradient.R` - Public interface
- 📄 `R/STgradient_helpers.R` - Helper functions
- 📄 `R/STgradient_legacy.R` - Legacy implementation
- ❌ **Missing**: `R/STgradient_core.R` - Core implementation
- ❌ **Missing**: Test suite

**5. SThet** - Spatial heterogeneity
- 📄 `R/SThet.R` - Public interface
- 📄 `R/SThet_invdist_test.R` - Distance-based tests
- ❌ **Missing**: Helper functions file
- ❌ **Missing**: Legacy implementation
- ❌ **Missing**: Core implementation (modular)
- ❌ **Missing**: Test suite

**6. STplot** - Spatial plotting
- 📄 `R/STplot.R` - Basic plotting
- 📄 `R/STplot_interpolation.R` - Interpolation-based
- ❌ **Missing**: Refactoring into modular structure
- ❌ **Missing**: Comprehensive test suite

---

## 📋 Priority Tasks

### Priority 1: Complete STenrich Refactoring 🔴 CRITICAL

**Status**: Implementation exists but not modularized

**Gap**: `STenrich_core.R` missing - the core logic is currently in `STenrich.R`

**Deliverables**:
- [ ] Create `R/STenrich_core.R`
  - Extract core logic from `STenrich.R`
  - Create internal functions:
    - `STenrich_core()` - Main core function
    - `STenrich_validate_input()` - Move from helpers
    - `STenrich_prepare_data()` - Move from helpers
    - `STenrich_calculate_gs_mean_exp()` - Move from helpers
    - `STenrich_calculate_gs_gsva_score()` - Move from helpers
    - `STenrich_permutation_test()` - Move from helpers
    - `STenrich_format_results()` - Move from helpers
- [ ] Update `R/STenrich.R` to call `STenrich_core()`
- [ ] Verify backward compatibility with legacy
- [ ] Create `tests/testthat/test-STenrich.R`
  - Core functionality tests
  - Input validation tests
  - Edge case tests
  - Comparison with legacy output

**Timeline**: ~3 days

---

### Priority 2: Complete STgradient Refactoring 🔴 CRITICAL

**Status**: Implementation exists but not modularized

**Gap**: `STgradient_core.R` missing

**Deliverables**:
- [ ] Create `R/STgradient_core.R`
  - Extract core logic from `STgradient.R`
  - Create internal functions:
    - `STgradient_core()` - Main core function
    - `STgradient_validate_input()` - Move from helpers
    - `STgradient_prepare_distances()` - Move from helpers
    - `STgradient_filter_neighbors()` - Move from helpers
    - `STgradient_summarize_distances()` - Move from helpers
    - `STgradient_identify_variable_genes()` - Move from helpers
    - `STgradient_extract_expression()` - Move from helpers
    - `STgradient_detect_outliers()` - Move from helpers
    - `STgradient_calculate_correlations()` - Move from helpers
    - `STgradient_format_results()` - Move from helpers
    - `STgradient_cleanup()` - Move from helpers
- [ ] Update `R/STgradient.R` to call `STgradient_core()`
- [ ] Verify backward compatibility with legacy
- [ ] Create `tests/testthat/test-STgradient.R`
  - Core functionality tests
  - Input validation tests
  - Edge case tests
  - Comparison with legacy output

**Timeline**: ~2 days

---

### Priority 3: Refactor SThet 🔴 CRITICAL

**Status**: Monolithic implementation, needs full modularization

**Deliverables**:
- [ ] Create `R/SThet.R` (modular version)
  - Public interface function
  - Call `SThet_core()` from core file
- [ ] Create `R/SThet_core.R`
  - Core implementation
  - `SThet_core()` - Main core function
  - `SThet_validate_input()` - Input validation
  - `SThet_calculate_moran()` - Moran's I calculation
  - `SThet_calculate_geary()` - Geary's C calculation
  - `SThet_format_results()` - Results formatting
- [ ] Create `R/SThet_legacy.R`
  - Replicate original implementation for comparison
- [ ] Update `R/SThet_invdist_test.R` to use new core
- [ ] Create `tests/testthat/test-SThet.R`
  - Core functionality tests
  - Input validation tests
  - Comparison with legacy output

**Timeline**: ~3 days

---

### Priority 4: Complete Test Suite 🟡 HIGH

**Status**: Only STclust has comprehensive tests (30/30 passing)

**Deliverables**:
- [ ] `tests/testthat/test-STenrich.R` - STenrich tests (~20-30 tests)
- [ ] `tests/testthat/test-STgradient.R` - STgradient tests (~15-20 tests)
- [ ] `tests/testthat/test-SThet.R` - SThet tests (~15-20 tests)
- [ ] `tests/testthat/test-STplot.R` - STplot tests (~10-15 tests)
- [ ] `tests/testthat/test-integration.R` - End-to-end workflow tests

**Test Strategy**:
1. **Unit tests** - Test individual functions
2. **Integration tests** - Test full workflows
3. **Regression tests** - Compare with legacy outputs
4. **Edge case tests** - Invalid inputs, missing data, etc.
5. **Performance tests** - Benchmark refactored vs legacy

**Timeline**: ~4 days

---

### Priority 5: Update NAMESPACE 🟡 HIGH

**Status**: Partially complete

**Deliverables**:
- [ ] Review and update `NAMESPACE` exports
- [ ] Ensure all public functions are exported
- [ ] Ensure all internal functions are NOT exported
- [ ] Add proper roxygen2 directives
- [ ] Verify `useDynLib(spatialGE, .registration = TRUE)` is correct

**Current Exports** (from NAMESPACE review):
- `STclust`, `STclust_legacy` ✅
- `STdiff`, `STdiff_*` (all subfunctions) ✅
- `STenrich`, `STenrich_legacy` ✅
- `STgradient`, `STgradient_legacy` ✅
- `SThet` ✅ (missing legacy export?)
- `STList`, `STList_legacy` ✅
- `STplot`, `STplot_interpolation` ✅
- Helper functions: `calculate_dist_matrices`, `calculate_weighted_dist`, etc.

**Timeline**: ~1 day

---

### Priority 6: Documentation 🟡 HIGH

**Status**: Basic roxygen2 docs exist

**Deliverables**:
- [ ] Update `README.md` with modular architecture
- [ ] Create vignette: "Spatial Clustering with STclust"
- [ ] Create vignette: "Spatial Differential Analysis with STdiff"
- [ ] Create vignette: "Spatial Enrichment with STenrich"
- [ ] Create vignette: "Gradient Analysis with STgradient"
- [ ] Create vignette: "Heterogeneity Analysis with SThet"
- [ ] Update function documentation (roxygen2)

**Timeline**: ~4 days

---

### Priority 7: Performance Optimization 🟢 MEDIUM

**Status**: Baseline established

**Deliverables**:
- [ ] Benchmark STenrich (modular vs legacy)
- [ ] Benchmark STgradient (modular vs legacy)
- [ ] Benchmark SThet (modular vs legacy)
- [ ] Identify bottlenecks
- [ ] Implement optimizations (parallelization, C++)

**Timeline**: ~2 days

---

### Priority 8: Seurat Integration 🟢 MEDIUM

**Status**: Basic integration exists (`seurat_helpers.cpp`)

**Deliverables**:
- [ ] Review and enhance Seurat integration
- [ ] Create wrapper functions for Seurat objects
- [ ] Add conversion utilities (STlist ↔ Seurat)
- [ ] Create examples using Seurat data

**Timeline**: ~2 days

---

## 📊 Success Metrics (End of Phase)

| Metric | Current | Target | Status |
|--------|---------|--------|--------|
| Refactored functions | 2/7 (STclust, STdiff) | 7/7 (100%) | 🏃 29% |
| Test coverage | ~40% (STclust only) | >80% | 🏃 In Progress |
| Vignettes | 0 | 5+ | ⏳ Pending |
| Documentation | Basic | Comprehensive | ⏳ Pending |
| Performance | Baseline | Optimized | ⏳ Pending |
| Backward compatibility | Verified (STclust, STdiff) | All functions | 🏃 In Progress |

---

## 🛠️ Refactoring Guidelines

### Core Function Pattern

```r
# Public interface (R/STfunction.R)
STfunction = function(x, samples, ...) {
  # Validation and prep
  validate_input <- STfunction_validate_input(...)
  prepared_data <- STfunction_prepare_data(...)
  
  # Call core
  core_result <- STfunction_core(prepared_data, ...)
  
  # Format and return
  format_results(core_result, ...)
}

# Core implementation (R/STfunction_core.R)
STfunction_core = function(data, ...) {
  # Main algorithm logic
  # No validation, no formatting
  # Just computation
}

# Helper functions (R/STfunction_helpers.R)
STfunction_validate_input = function(...) { ... }
STfunction_prepare_data = function(...) { ... }
STfunction_calculate_step1 = function(...) { ... }
STfunction_calculate_step2 = function(...) { ... }

# Legacy (R/STfunction_legacy.R)
STfunction_legacy = function(x, samples, ...) {
  # Exact copy of original FridleyLab/spatialGE implementation
}
```

### Testing Pattern

```r
# tests/testthat/test-STfunction.R
test_that("STfunction basic functionality", {
  # Test with example data
  result <- STfunction(st_obj, ...)
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
})

test_that("STfunction matches legacy output", {
  result_modern <- STfunction(st_obj, ...)
  result_legacy <- STfunction_legacy(st_obj, ...)
  expect_equal(result_modern, result_legacy, tolerance = 1e-10)
})

test_that("STfunction input validation", {
  expect_error(STfunction(st_obj, invalid_param = TRUE))
})

test_that("STfunction edge cases", {
  # Missing data, empty inputs, etc.
})
```

---

## 📝 Notes

1. **Legacy functions are ground truth** - All legacy implementations must match FridleyLab/spatialGE exactly for reproducibility
2. **Test before refactor** - Always establish test baseline before modifying code
3. **Validate outputs** - Compare refactored outputs with legacy outputs (should be identical within numerical tolerance)
4. **Commit frequently** - Push changes after each subagent completion
5. **Document decisions** - Record any deviations from original implementation

---

*This TASKS.md serves as the master development roadmap for spatialGE refactoring. Refer to it throughout development.*

*Last updated: 2026-03-19 13:35 UTC*
