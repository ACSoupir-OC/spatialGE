# STclust Test Expansion Plan

**Goal:** Expand STclust test suite from 6 tests (216 lines) to 15+ tests (>500 lines) matching STdiff comprehensiveness

**Timeline:** ~15 minutes (3 subagents × 5 min each)

---

## Current State

- **File:** `tests/testthat/test-stclust.R`
- **Tests:** 6 tests, 216 lines
- **Coverage:** ~75%
- **Missing:** Edge cases, parameter variations, helper function isolation

---

## Subagent Tasks

### Subagent 1: Core Function Tests (150 lines)

**File:** `tests/testthat/test-stclust-comprehensive.R` (Section 1)

**Tests to add:**
1. **STclust_select_genes()** - Gene selection correctness
   - Test VST calculation
   - Test topgenes filtering
   - Test multi-sample gene selection
   - Test count filtering

2. **STclust_calculate_distances()** - Distance matrix accuracy
   - Test euclidean distance
   - Test manhattan distance
   - Test correlation distance
   - Test spatial distance calculation
   - Test scaling correctness

3. **STclust_weight_distances()** - Weight combination
   - Test ws=0 (expression only)
   - Test ws=1 (spatial only)
   - Test ws=0.025 (balanced)
   - Test multiple ws values

**Verification:**
- [ ] File created with 3+ tests
- [ ] All helper functions tested in isolation
- [ ] Tests run without errors
- [ ] Commit made: "test: Add STclust core function tests"

---

### Subagent 2: Clustering & Parameter Tests (200 lines)

**File:** `tests/testthat/test-stclust-comprehensive.R` (Section 2)

**Tests to add:**
4. **STclust_hierarchical()** - Clustering logic
   - Test DTC method (dynamicTreeCut)
   - Test fixed k method
   - Test different linkage methods (ward.D2, average, complete)
   - Test deepSplit parameter variations

5. **Parameter combinations** - Integration tests
   - Test ks=2 (single fixed k)
   - Test ks=2:5 (range of k values)
   - Test ws vector (multiple weights)
   - Test ks='dtc' with deepSplit=2,3,4

6. **Multiple samples** - Batch processing
   - Test single sample
   - Test 2-3 samples
   - Test sample-specific clustering

**Verification:**
- [ ] File appended with 3+ tests
- [ ] All parameter combinations tested
- [ ] Tests run without errors
- [ ] Commit made: "test: Add STclust parameter and clustering tests"

---

### Subagent 3: Edge Cases & Error Handling (200 lines)

**File:** `tests/testthat/test-stclust-comprehensive.R` (Section 3)

**Tests to add:**
7. **Edge cases** - Boundary conditions
   - Single gene (minimum)
   - Two spots (minimum for clustering)
   - All zeros in expression
   - Identical spots (degenerate case)
   - Very large ws (spatial dominated)
   - Very small ws (expression dominated)

8. **Error handling** - Invalid inputs
   - Invalid ks values (negative, zero)
   - Invalid ws values (<0, >1)
   - Missing samples
   - Invalid dist_metric
   - Invalid linkage method
   - NULL STlist

9. **Regression tests** - Consistency
   - Deterministic results (same seed = same output)
   - STclust() vs STclust_legacy() with various parameters
   - Reproducibility across runs

10. **Integration tests** - Full workflow
    - STclust → STdiff pipeline
    - STclust → STenrich pipeline
    - STclust → plotting functions

**Verification:**
- [ ] File appended with 4+ tests
- [ ] All edge cases covered
- [ ] Error handling validated
- [ ] Commit made: "test: Add STclust edge cases and error handling"

---

## Expected Final State

**File:** `tests/testthat/test-stclust-comprehensive.R`
- **Tests:** 15+ (up from 6)
- **Lines:** 500-600 (up from 216)
- **Coverage:** >80% (up from ~75%)

**Sections:**
1. Core function tests (3 tests)
2. Clustering & parameter tests (3 tests)
3. Edge cases & error handling (4+ tests)
4. Existing tests preserved (6 tests)

---

## Success Criteria

- [ ] All 15+ tests pass
- [ ] No regressions in existing tests
- [ ] Code coverage increases to >80%
- [ ] All commits pushed to GitHub
- [ ] TASKS.md updated

---

## Execution Protocol

**After user approval:**
1. Spawn Subagent 1 → Wait → Verify → Commit
2. Spawn Subagent 2 → Wait → Verify → Commit
3. Spawn Subagent 3 → Wait → Verify → Commit
4. Update TASKS.md
5. Create completion summary

**Do NOT ask for approval between subagents** - proceed autonomously after initial approval.
