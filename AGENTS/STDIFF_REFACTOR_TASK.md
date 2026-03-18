# STdiff Refactoring Task

**Priority:** High  
**Estimated Effort:** 2-3 days  
**Risk:** Medium (complex function, many edge cases)

---

## Objective

Refactor `STdiff.R` (42KB monolith) into modular, maintainable functions while preserving all existing functionality and ensuring backward compatibility.

---

## Current State

**File:** `R/STdiff.R` (~42KB, ~900 lines)

**Issues:**
- Monolithic function mixing multiple test types
- Complex annotation handling scattered throughout
- Hard to test individual components
- Difficult to maintain and extend

---

## Target Architecture

### New File Structure

```
R/
├── STdiff.R                  # Main dispatcher (thin wrapper)
├── STdiff_core.R             # Already exists - keep
├── STdiff_helpers.R          # Already exists - keep
├── STdiff_volcano.R          # Already exists - keep
├── STdiff_nonspatial.R       # NEW: Non-spatial tests
├── STdiff_spatial.R          # NEW: Spatial mixed models
├── STdiff_results.R          # NEW: Result compilation
└── STdiff_legacy.R           # Keep for backward compatibility
```

### New Functions

**STdiff_nonspatial.R:**
- `STdiff_run_nonspatial()` - Main entry point for non-spatial tests
- `STdiff_fit_mm()` - Mixed models
- `STdiff_fit_ttest()` - T-tests
- `STdiff_fit_wilcoxon()` - Wilcoxon tests
- `STdiff_adjust_pvalues()` - P-value adjustment

**STdiff_spatial.R:**
- `STdiff_run_spatial()` - Main entry point for spatial models
- `STdiff_fit_spatial_models()` - Fit spaMM models with Matern covariance
- `STdiff_check_convergence()` - Check model convergence
- `STdiff_extract_spatial_pvalues()` - Extract spatial p-values

**STdiff_results.R:**
- `STdiff_compile_results()` - Compile all results into final format
- `STdiff_format_output()` - Format for user consumption
- `STdiff_add_metadata()` - Add cluster/sample metadata

---

## Subagent Tasks

### Subagent 1: Create STdiff_nonspatial.R (10-15 min)
**Task:** Extract non-spatial testing logic into dedicated module

**Deliverables:**
- `R/STdiff_nonspatial.R` with 4-5 functions
- All non-spatial test logic (mm, t_test, wilcoxon)
- P-value adjustment logic
- Roxygen2 documentation

**Verification:**
- [ ] File exists at `R/STdiff_nonspatial.R`
- [ ] Contains `STdiff_run_nonspatial()` function
- [ ] Contains test-specific functions (mm, ttest, wilcoxon)
- [ ] All functions have Roxygen2 docs
- [ ] Code imports correctly: `source("R/STdiff_nonspatial.R")`

---

### Subagent 2: Create STdiff_spatial.R (10-15 min)
**Task:** Extract spatial model fitting logic into dedicated module

**Deliverables:**
- `R/STdiff_spatial.R` with 3-4 functions
- Spatial model fitting with spaMM
- Convergence checking
- P-value extraction

**Verification:**
- [ ] File exists at `R/STdiff_spatial.R`
- [ ] Contains `STdiff_run_spatial()` function
- [ ] Contains `STdiff_fit_spatial_models()`
- [ ] Handles convergence failures gracefully
- [ ] All functions have Roxygen2 docs

---

### Subagent 3: Create STdiff_results.R (5-10 min)
**Task:** Extract result compilation logic into dedicated module

**Deliverables:**
- `R/STdiff_results.R` with 2-3 functions
- Result compilation from non-spatial and spatial tests
- Output formatting
- Metadata addition

**Verification:**
- [ ] File exists at `R/STdiff_results.R`
- [ ] Contains `STdiff_compile_results()` function
- [ ] Produces identical output to current implementation
- [ ] All functions have Roxygen2 docs

---

### Subagent 4: Update STdiff.R Dispatcher (5-10 min)
**Task:** Refactor main STdiff.R to use new modular functions

**Deliverables:**
- Updated `R/STdiff.R` (thin wrapper, ~100 lines)
- Calls new modular functions
- Maintains same API for users
- Backward compatible

**Verification:**
- [ ] `STdiff.R` is significantly smaller (<200 lines)
- [ ] Same function signature as before
- [ ] Calls `STdiff_run_nonspatial()`, `STdiff_run_spatial()`, `STdiff_compile_results()`
- [ ] Backward compatible with existing code

---

### Subagent 5: Update NAMESPACE and Exports (5 min)
**Task:** Update NAMESPACE with new function exports

**Deliverables:**
- Updated `NAMESPACE` file
- Export all new public functions
- Keep internal functions unexported

**Verification:**
- [ ] NAMESPACE updated with new exports
- [ ] `STdiff_run_nonspatial` exported (if public)
- [ ] `STdiff_run_spatial` exported (if public)
- [ ] `STdiff_compile_results` exported (if public)
- [ ] Internal helper functions not exported

---

### Subagent 6: Verification Testing (15-20 min)
**Task:** Run comprehensive tests to ensure identical behavior

**Deliverables:**
- Test results showing identical output
- All existing tests pass
- New modular tests pass

**Verification:**
- [ ] Run `tests/testthat/test-stdiff-comprehensive.R`
- [ ] Compare output with saved reference results
- [ ] Verify regression tests pass (STdiff vs STdiff_legacy)
- [ ] Document any differences

---

## Success Criteria

1. **Functionality Preserved:** All tests pass, identical output to original
2. **Code Reduction:** `STdiff.R` reduced from ~900 lines to <200 lines
3. **Modularity:** Clear separation of concerns (non-spatial, spatial, results)
4. **Documentation:** All new functions have Roxygen2 docs
5. **Backward Compatibility:** Existing code using `STdiff()` continues to work
6. **Test Coverage:** >80% coverage for new modules

---

## Risk Mitigation

**Risk:** Breaking changes to API  
**Mitigation:** Keep main `STdiff()` signature unchanged, only refactor internals

**Risk:** Loss of functionality  
**Mitigation:** Comprehensive regression tests comparing old vs new

**Risk:** Convergence issues in spatial models  
**Mitigation:** Preserve existing error handling, add logging

**Risk:** Performance regression  
**Mitigation:** Benchmark before/after, ensure no significant slowdown

---

## References

- `R/STdiff.R` - Current implementation
- `R/STdiff_core.R` - Existing core functions
- `R/STdiff_helpers.R` - Existing helpers
- `tests/testthat/test-stdiff-comprehensive.R` - Test suite
- `AGENTS/TASKS.md` - Overall project tasks

---

**Start Date:** 2026-03-18  
**Target Completion:** 2026-03-20
