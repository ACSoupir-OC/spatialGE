# spatialGE Legacy Compatibility Verification

**Date:** 2026-04-09  
**Goal:** Ensure all refactored functions match legacy implementations, legacy functions remain exposed, documentation references legacy functions, and all tests pass

---

## Overview

spatialGE v2.0.0 refactored core functions (STclust, STdiff) into modular architecture while maintaining backward compatibility. This task verifies:
1. Refactored functions produce identical results to legacy
2. Legacy functions remain exported and functional
3. Documentation mentions legacy alternatives
4. All tests pass

---

## Subagent Tasks

### Subagent 1: Audit Function Exports

**Task:** Verify all legacy functions are exported in NAMESPACE and documented

**Steps:**
1. Read `NAMESPACE` file - list all exported functions
2. Read `man/` directory - list all .Rd files
3. Identify legacy functions: `STclust_legacy`, `STdiff_legacy`, `STenrich_legacy`, `STgradient_legacy`, `SThet_legacy`
4. Verify each legacy function is in NAMESPACE (exported)
5. Verify each legacy function has .Rd documentation
6. Check DESCRIPTION for version 2.0.0

**Verification Checklist:**
- [ ] `STclust_legacy` exported in NAMESPACE
- [ ] `STdiff_legacy` exported in NAMESPACE
- [ ] `STenrich_legacy` exported in NAMESPACE (if exists)
- [ ] `STgradient_legacy` exported in NAMESPACE (if exists)
- [ ] `SThet_legacy` exported in NAMESPACE (if exists)
- [ ] All legacy functions have .Rd files in `man/`
- [ ] All legacy functions have @export tag in source

**Output:** `verification_exports.md`

---

### Subagent 2: Compare Refactored vs Legacy Implementations

**Task:** Verify refactored functions produce identical results to legacy

**Steps:**
1. Read refactored source files:
   - `R/STclust_core.R`, `R/STclust_helpers.R`
   - `R/STdiff_core.R`, `R/STdiff_helpers.R`, `R/STdiff_spatial.R`, `R/STdiff_nonspatial.R`
2. Read legacy source files:
   - `R/STclust_legacy.R`
   - `R/STdiff_legacy.R`
3. Compare algorithmic steps (not line-by-line, but logical flow)
4. Identify any differences in:
   - Default parameters
   - Statistical methods
   - Edge case handling
   - Return value structure

**Verification Checklist:**
- [ ] STclust: Core algorithm matches legacy
- [ ] STclust: Default parameters identical
- [ ] STclust: Return value structure identical
- [ ] STdiff: Non-spatial tests match legacy
- [ ] STdiff: Spatial tests match legacy
- [ ] STdiff: Default parameters identical
- [ ] STdiff: Return value structure identical

**Output:** `verification_algorithm_comparison.md`

---

### Subagent 3: Run Test Suite

**Task:** Execute all tests and verify they pass

**Steps:**
1. Load conda R environment: `/home/node/.openclaw/workspace/miniconda3/envs/spatialge-r`
2. Run `devtools::load_all()` to load package
3. Run `testthat::test_local()` or `devtools::test()`
4. Capture test results (pass/fail/skip)
5. Identify any failing tests
6. For failing tests, determine if:
   - Test is outdated
   - Implementation changed intentionally
   - Bug introduced in refactoring

**Verification Checklist:**
- [ ] All tests load without errors
- [ ] Test count matches expected (~25-30 test files)
- [ ] All tests pass (or document expected failures)
- [ ] Legacy function tests pass
- [ ] Refactored function tests pass

**Output:** `verification_test_results.md` with full test output

---

### Subagent 4: Documentation Audit

**Task:** Verify documentation mentions legacy functions and provides migration guidance

**Steps:**
1. Read .Rd files for refactored functions:
   - `man/STclust.Rd`
   - `man/STdiff.Rd`
   - `man/STenrich.Rd`
   - `man/STgradient.Rd`
   - `man/SThet.Rd`
2. Check if each mentions:
   - Legacy function name (e.g., "For the legacy implementation, see `STclust_legacy`")
   - When to use legacy vs refactored
   - Any behavioral differences
3. Read source .R files for @seealso, @inheritParams tags
4. Verify NEWS.md mentions refactoring and backward compatibility

**Verification Checklist:**
- [ ] STclust documentation mentions STclust_legacy
- [ ] STdiff documentation mentions STdiff_legacy
- [ ] STenrich documentation mentions STenrich_legacy (if exists)
- [ ] STgradient documentation mentions STgradient_legacy (if exists)
- [ ] NEWS.md documents refactoring changes
- [ ] Migration guide or notes exist for users

**Output:** `verification_documentation.md`

---

### Subagent 5: Numerical Equivalence Test

**Task:** Run direct numerical comparison between legacy and refactored functions

**Steps:**
1. Create test script that:
   - Loads TNBC or melanoma test data
   - Runs both legacy and refactored versions
   - Compares outputs numerically
2. Test cases:
   - STclust: Same clusters produced?
   - STdiff: Same p-values, logFC, statistics?
   - Allow for floating-point tolerance (1e-10)
3. Document any differences found

**R Script Structure:**
```r
library(spatialGE)
data(tnbc)  # or load test data

# STclust comparison
result_legacy <- STclust_legacy(tnbc, ...)
result_refactored <- STclust(tnbc, ...)
identical(result_legacy$clusters, result_refactored$clusters)

# STdiff comparison
result_legacy <- STdiff_legacy(tnbc, ...)
result_refactored <- STdiff(tnbc, ...)
all.equal(result_legacy$table, result_refactored$table, tolerance=1e-10)
```

**Verification Checklist:**
- [ ] STclust produces identical clusters
- [ ] STdiff produces identical statistics (within tolerance)
- [ ] Any differences documented and justified

**Output:** `verification_numerical_equivalence.md` with test script and results

---

## Pre-Subagent Setup

**Environment:**
```bash
# Use spatialge-r conda environment
export R_LIBS=/home/node/.openclaw/workspace/miniconda3/envs/spatialge-r/lib/R/library
```

**Working Directory:** `/home/node/.openclaw/workspace/spatialGE`

**Test Data:** TNBC dataset (should be in `tests/testthat/data/` or via Git LFS)

---

## Post-Subagent Verification

After all subagents complete:

1. **Review all verification reports**
2. **Fix any issues found:**
   - Missing exports → Add to NAMESPACE
   - Missing documentation → Add @seealso tags
   - Failing tests → Fix or document
   - Numerical differences → Investigate and fix
3. **Commit fixes** with message: "verify: Legacy compatibility verification + fixes"
4. **Update TASK.md** with completion status

---

## Success Criteria

- [ ] All 5 subagents complete successfully
- [ ] All legacy functions exported and documented
- [ ] All tests pass
- [ ] Documentation references legacy functions
- [ ] Numerical equivalence verified (or differences documented)
- [ ] All fixes committed to git

---

## Timeline

**Estimated:** 2-3 hours total
- Subagent 1 (Exports): 10 min
- Subagent 2 (Algorithm): 20 min
- Subagent 3 (Tests): 30-45 min
- Subagent 4 (Documentation): 15 min
- Subagent 5 (Numerical): 30-45 min
- Review + fixes: 30-60 min

---

*Created: 2026-04-09 | spatialGE v2.0.0 Legacy Verification*
