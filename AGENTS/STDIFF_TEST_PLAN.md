# Comprehensive STdiff Tests

**Test File:** `tests/testthat/test-stdiff-comprehensive.R`  
**Created:** 2026-03-18  
**Coverage Target:** >80% for STdiff module

---

## Test Structure

### 1. Core Function Tests
- `STdiff_select_genes()` - Gene filtering and non-spatial testing
- `STdiff_run_nonspatial()` - Non-spatial model fitting
- `STdiff_run_spatial()` - Spatial mixed models
- `STdiff_compile_results()` - Result compilation

### 2. Helper Function Tests
- `prepare_stdiff_combo()` - Combination preparation
- `count_cores()` - Core detection
- `expandSparse()` - Sparse matrix expansion
- `non_spatial_de()` - Non-spatial DE testing
- `spatial_de()` - Spatial model fitting
- `raise_err()` - Error handling

### 3. Integration Tests
- Full STdiff workflow
- Different test types (mm, t_test, wilcoxon)
- Pairwise vs reference-based testing
- Parallelization

### 4. Edge Cases
- Missing annotations
- Invalid parameters
- Small sample sizes
- Convergence failures

### 5. Regression Tests
- Compare with STdiff_legacy()
- Saved reference results

---

## Implementation Plan

### Subagent Tasks

**Subagent 1: Core Function Tests** (5-7 min)
- Test `STdiff_select_genes()` with various parameters
- Test gene filtering logic
- Test annotation validation
- Test combination preparation

**Subagent 2: Helper Function Tests** (5-7 min)
- Test `non_spatial_de()` with different test types
- Test `spatial_de()` model fitting
- Test `expandSparse()` matrix conversion
- Test error handling functions

**Subagent 3: Integration Tests** (10-15 min)
- Full STdiff workflow tests
- Compare test types (mm, t_test, wilcoxon)
- Test pairwise comparisons
- Test parallelization

**Subagent 4: Edge Cases & Error Handling** (5-7 min)
- Invalid parameter handling
- Missing data handling
- Convergence failure handling
- Boundary conditions

**Subagent 5: Regression Tests** (10-15 min)
- Compare STdiff() vs STdiff_legacy()
- Load and compare saved reference results
- Ensure reproducibility across runs

---

## Test Data Requirements

**Minimal Test Data:**
- 2-3 samples
- 100-200 spots per sample
- 500-1000 genes
- Cluster annotations (STclust results)

**Can Use:**
- Melanoma Thrane dataset (subset)
- TNBC Bassiouni dataset (subset)
- Simulated data for edge cases

---

## Expected Coverage

| Function | Lines | Target Coverage | Priority |
|----------|-------|-----------------|----------|
| STdiff_core.R | ~400 | 85% | High |
| STdiff_helpers.R | ~300 | 85% | High |
| STdiff_volcano.R | ~200 | 75% | Medium |
| STdiff_legacy.R | ~900 | 60% | Low (legacy) |
| **Total** | **~1800** | **~80%** | **High** |

---

## Success Criteria

- [ ] All tests pass on CI/CD
- [ ] Coverage >80% for core functions
- [ ] No false positives/negatives
- [ ] Tests complete in <5 minutes
- [ ] Edge cases properly handled
- [ ] Regression tests confirm identical results to legacy

---

## Notes

- Use `devtools::load_all()` for testing development version
- Skip slow tests on CI with `skip_on_cis()`
- Use temporary directories for test outputs
- Clean up after tests (remove temporary files)
