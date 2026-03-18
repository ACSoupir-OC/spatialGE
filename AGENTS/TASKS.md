# spatialGE Development Tasks

**Package Version:** 2.0.0  
**Repository:** https://github.com/ACSoupir-OC/spatialGE  
**Last Review:** 2026-03-18  
**Status:** Active Development - Refactoring Phase

---

## 🎯 Current Priorities

### Priority 1: Complete Refactoring (Q2 2026)

| Task | Status | Effort | Notes |
|------|--------|--------|-------|
| Modularize STdiff (42KB monolith) | 📋 Planned | High | Split into 3-4 functions |
| Replace DynamicTreeCut with fastcluster | ⚠️ In Progress | Low | 20-35% speedup expected |
| Vectorize distance calculations | 📋 Planned | Medium | Matrix ops instead of loops |
| Standardize error handling | 📋 Planned | Medium | Across all modules |
| Windows-compatible parallelization | 📋 Planned | Medium | Replace mclapply |

### Priority 2: Test Coverage (Q2 2026)

| Task | Status | Target | Current |
|------|--------|--------|---------|
| Core functions (STlist, SThet, STclust) | ✅ Done | 80% | ~75% |
| STdiff (all test types) | 📋 Planned | 80% | ~40% |
| STenrich (GSVA + avg paths) | 📋 Planned | 80% | ~30% |
| STgradient (all modes) | ✅ Done | 80% | ~70% |
| Plotting functions | 📋 Planned | 60% | ~20% |
| **Overall Coverage** | 🏃 In Progress | **>75%** | **~50%** |

### Priority 3: Performance Optimization (Q3 2026)

| Optimization | Expected Speedup | Status |
|--------------|------------------|--------|
| fastcluster integration | 20-35% | ⚠️ Partial |
| Pre-compute dendrograms | 5-10% | 📋 Planned |
| Vectorized distance calc | 15-20% | 📋 Planned |
| Rcpp distance calculation | 20-30% | 📋 Future |
| Sparse matrix optimization | 3-5% | 📋 Planned |

---

## 📋 Detailed Task Breakdown

### Module: STdiff (High Priority)

**Current State:** 42KB monolithic function (~900 lines)

**Tasks:**
- [ ] **STdiff-01:** Split `STdiff.R` into modular functions
  - `STdiff_run_nonspatial()` - Non-spatial tests (mm, t_test, wilcoxon)
  - `STdiff_run_spatial()` - Spatial mixed models (spaMM)
  - `STdiff_compile_results()` - Result formatting
  - `STdiff_helpers.R` - Annotation recoding, sparse matrix helpers
- [ ] **STdiff-02:** Add comprehensive test coverage
  - Test each test type separately
  - Test annotation handling
  - Test sparse matrix expansion
  - Test parallelization
- [ ] **STdiff-03:** Simplify annotation recoding logic
  - Extract to dedicated helper function
  - Add unit tests for edge cases
- [ ] **STdiff-04:** Document function parameters consistently
  - Roxygen2 documentation for all new functions
  - Examples for each test type

**Estimated Effort:** 2-3 days  
**Risk:** Medium (complex function, many edge cases)

---

### Module: STenrich (Medium Priority)

**Current State:** Mixed GSVA/average paths, memory concerns

**Tasks:**
- [ ] **STenrich-01:** Separate GSVA calculation to dedicated module
  - `STenrich_gsva.R` - GSVA score calculation
  - `STenrich_average.R` - Average expression path
- [ ] **STenrich-02:** Implement out-of-memory distance computation
  - Use HDF5Array for large datasets
  - Test with >10K spots
- [ ] **STenrich-03:** Simplify permutation testing
  - Clearer loop structure
  - Progress bar integration
- [ ] **STenrich-04:** Add memory profiling
  - Benchmark with different dataset sizes
  - Document memory requirements

**Estimated Effort:** 1-2 days  
**Risk:** Low

---

### Module: STclust (Medium Priority)

**Current State:** Functional, performance bottlenecks identified

**Tasks:**
- [ ] **STclust-01:** Integrate fastcluster replacement
  - Replace `hclust()` with `fastcluster::hclust()`
  - Validate identical output
  - Benchmark performance gain
- [ ] **STclust-02:** Pre-compute dendrogram for multiple ws/k
  - Reuse across parameter sweeps
  - Add parameter to control behavior
- [ ] **STclust-03:** Separate distance matrix calculation
  - `STclust_distances.R` - Distance computation
  - `STclust_cluster.R` - Clustering logic
- [ ] **STclust-04:** Parameterize weight defaults
  - Better documentation for `ws` parameter
  - Add vignette with tuning guidance

**Estimated Effort:** 1 day  
**Risk:** Low

---

### Module: SThet (Low Priority)

**Current State:** Functional, minor cleanup needed

**Tasks:**
- [ ] **SThet-01:** Extract weight calculation to separate function
  - `SThet_weights.R` - k-NN and distance-based weights
- [ ] **SThet-02:** Simplify result storage
  - Vectorized operations where possible
- [ ] **SThet-03:** Add Windows-compatible parallelization
  - `foreach` + `doParallel` fallback
- [ ] **SThet-04:** Standardize gene_meta access
  - Helper functions for column access
  - Reduce hard-coded column names

**Estimated Effort:** 0.5 days  
**Risk:** Low

---

### Module: STgradient (Low Priority)

**Current State:** Recently refactored, well-documented

**Tasks:**
- [ ] **STgradient-01:** Complete test coverage
  - Test all distance modes (min, average)
  - Test outlier removal
  - Test robust regression
- [ ] **STgradient-02:** Document edge cases
  - Small sample sizes
  - Missing reference domain
  - Extreme outliers

**Estimated Effort:** 0.5 days  
**Risk:** Low

---

### Infrastructure (Ongoing)

**Tasks:**
- [ ] **INFRA-01:** Set up codecov integration
  - Add coverage badge to README
  - Configure GitHub Actions
- [ ] **INFRA-02:** Add benchmark regression testing
  - Track performance over time
  - Alert on >10% regression
- [ ] **INFRA-03:** Update pkgdown documentation
  - Add function reference pages
  - Include performance notes
- [ ] **INFRA-04:** Create developer vignette
  - Package architecture overview
  - Contributing guidelines
  - Testing best practices

**Estimated Effort:** 1-2 days  
**Risk:** Low

---

## 📊 Test Coverage Goals

### Current Coverage by Module

| Module | Files | Lines | Covered | % |
|--------|-------|-------|---------|---|
| STlist | 3 | ~1200 | ~900 | 75% |
| SThet | 1 | ~300 | ~225 | 75% |
| STclust | 4 | ~800 | ~600 | 75% |
| STdiff | 5 | ~2000 | ~800 | 40% |
| STenrich | 3 | ~900 | ~270 | 30% |
| STgradient | 3 | ~700 | ~490 | 70% |
| Plotting | 8 | ~1200 | ~240 | 20% |
| Utils | 4 | ~600 | ~480 | 80% |
| **Total** | **31** | **~7700** | **~4005** | **~52%** |

### Target Coverage (End of Q2 2026)

| Module | Target % | Priority |
|--------|----------|----------|
| STlist | 85% | High |
| SThet | 85% | Medium |
| STclust | 85% | High |
| STdiff | 80% | High |
| STenrich | 80% | Medium |
| STgradient | 85% | Medium |
| Plotting | 60% | Low |
| Utils | 90% | Low |
| **Overall** | **>75%** | **High** |

---

## 🚀 Performance Goals

### Current Bottlenecks (100 cells × 500 genes)

| Function | Time (ms) | % of Total | Optimization Target |
|----------|-----------|------------|---------------------|
| DynamicTreeCut | 6.30 | 45.5% | Replace with fastcluster |
| calculate_dist_matrices | 4.80 | 34.6% | Vectorize or Rcpp |
| hclust | 0.34 | 2.4% | fastcluster (2-5x) |
| Other | 2.46 | 17.5% | Minor optimizations |

### Performance Targets (End of Q3 2026)

| Metric | Current | Target | Improvement |
|--------|---------|--------|-------------|
| STclust (100 cells) | ~14 ms | ~6 ms | 55% faster |
| STclust (1000 cells) | ~1.4 s | ~0.6 s | 55% faster |
| STdiff (100 genes) | ~30 s | ~20 s | 33% faster |
| STenrich (100 perms) | ~15 s | ~10 s | 33% faster |

---

## 📅 Timeline

### Q2 2026 (April - June)
- **Week 1-2:** STdiff refactoring
- **Week 3-4:** Test coverage expansion (STdiff, STenrich)
- **Week 5-6:** fastcluster integration
- **Week 7-8:** Vectorization of distance calculations
- **Week 9-10:** Infrastructure (codecov, benchmarks)
- **Week 11-12:** Documentation, bug fixes, release prep

### Q3 2026 (July - September)
- **Week 1-4:** Advanced optimizations (Rcpp if needed)
- **Week 5-8:** New features (based on user feedback)
- **Week 9-12:** CRAN submission prep

---

## 📝 Notes

### Known Issues
1. **mclapply incompatibility:** Windows users cannot use parallel features
2. **Memory usage:** STenrich can consume >8GB for large datasets
3. **Hard-coded paths:** Some test files have absolute paths (being fixed)

### Recent Improvements (2026-03-18)
- ✅ Cloned from GitHub
- ✅ Reviewed package structure
- ✅ Created AGENTS/ and archive/ folders
- ✅ Added archive/ to .gitignore
- ✅ Fixed hardcoded paths in test scripts
- ✅ Created comprehensive task list

### Dependencies
- **R (>= 4.3.0)** - Required for latest features
- **spaMM** - Spatial mixed models (STdiff)
- **DynamicTreeCut** - Being replaced with fastcluster
- **GSVA** - Gene set enrichment (STenrich)
- **hdf5r** - HDF5 support for large datasets

---

**Next Review:** 2026-04-01  
**Contact:** Alex Soupir <Alex.Soupir@moffitt.org>

---

*This file is updated as tasks are completed. Check git history for changes.*
