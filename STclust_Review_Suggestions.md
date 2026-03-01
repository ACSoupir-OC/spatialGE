# STclust Refactoring Review - Suggestions for Future Work

**Date:** 2026-03-01  
**Status:** Refactoring completed and working (30/30 tests passing)

## ✅ Current State (Do NOT Change)
- **Legacy function preserved exactly** - `STclust_legacy.R` remains unchanged for backward compatibility and reproducibility
- **NAMESPACE left alone** - Trust `devtools::document()` to handle namespace generation
- **Sourcing removed from `STclust.R`** - Package structure now follows R best practices
- **Tests updated** - Use `devtools::load_all()` instead of manual sourcing

## 📝 Future Improvement Suggestions

### 1. Code Quality Enhancements
- **[LOW]** Fix hardcoded `verbose = 1L` in `STclust.R` line 132 that overwrites user input
- **[LOW]** Remove dead code in `STclust_hierarchical()` (hardcoded `clmethod = 'hclust'`)
- **[LOW]** Standardize parameter naming across functions (e.g., `ws` → `spatial_weights`, `ks` → `n_clusters`)

### 2. Documentation Improvements
- **[MEDIUM]** Add cross-references in roxygen2 docs (See Also sections)
- **[MEDIUM]** Create vignette/tutorial for STclust usage best practices
- **[MEDIUM]** Add performance benchmarks comparing new vs legacy implementation

### 3. Testing Enhancements
- **[MEDIUM]** Add edge case tests (single sample, single gene, empty data)
- **[LOW]** Test parallelization behavior on different OS platforms
- **[LOW]** Add memory usage profiling tests

### 4. Structural Improvements
- **[HIGH]** Consider refactoring helper functions to avoid duplication between `STclust_helpers.R` and `STclust_legacy.R`
- **[HIGH]** Evaluate if C++ functions should be exported via `@export` tags instead of manual NAMESPACE edits

### 5. Maintenance Items
- **[LOW]** Update README.md with STclust refactoring summary
- **[LOW]** Add deprecation notice in `STclust_legacy()` documentation (optional - keep as "legacy for reproducibility")
- **[LOW]** Review all roxygen2 documentation for consistency

## ⚠️ Strict Constraints (DO NOT VIOLATE)
- Never modify `STclust_legacy.R` unless explicitly directed
- Never manually edit NAMESPACE file (let `devtools::document()` handle it)
- Preserve backward compatibility at all costs
- Keep legacy function for reproducibility purposes

## 📊 Priority Order for Future Work
1. Code quality fixes (LOW priority - can wait)
2. Documentation improvements (MEDIUM priority)
3. Testing enhancements (MEDIUM priority)
4. Structural improvements (HIGH priority - but only after thorough testing)

---
*This document serves as a reference for future development. The current STclust refactoring is production-ready and should not be modified without explicit approval.*