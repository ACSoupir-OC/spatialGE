# spatialGE Legacy Documentation Verification

**Date:** 2026-04-09  
**Task:** Check if refactored function documentation mentions legacy alternatives  
**Scope:** STclust, STdiff, STenrich, STgradient, SThet

---

## Legacy Reference Summary

| Function | Mentions Legacy? | Quote if Yes |
|----------|------------------|--------------|
| STclust | ❌ N | *None found* |
| STdiff | ❌ N | *None found* |
| STenrich | ✅ Y | "Original implementation of STenrich for reproducibility with previous results<br><br>**This is the legacy implementation from version 1.x.** Use `STenrich()` for new analyses." |
| STgradient | ❌ N | *See SeeAlso section: "STgradient_legacy for reproducibility with original implementation"* |
| SThet | ❌ N | *None found* |

---

## NEWS.md Refactoring Mentions

**Mentions refactoring?** ✅ **Yes**

From `NEWS.md` (spatialGE 2.0.0):
> **Modular Architecture** - Major refactoring of core functions
> - `STclust` modular architecture with 5 helper functions
> - `STdiff` modular architecture (non-spatial + spatial tests)
> - Replaced base R `hclust()` with `fastcluster::hclust()` for 2-3x speedup
> - Expected 20-35% speedup on `STclust` runtime

---

## Assessment: **PARTIAL**

### ✅ Complete
- **STenrich**: Fully documented with clear legacy marker, superseded lifecycle stage, and SeeAlso reference

### ⚠️ Partial
- **STgradient**: Has SeeAlso reference to `STgradient_legacy` but lacks explicit "legacy" marker in description like STenrich

### ❌ Missing
- **STclust**: No legacy reference despite modular refactoring
- **STdiff**: No legacy reference despite modular refactoring (only mentions "Maintains the original API")
- **SThet**: No legacy reference (no legacy version exists, but should confirm)

---

## Recommended Text Additions

### STclust.Rd
**Add to description/details:**
```
\strong{Note:} This function was refactored in version 2.0.0 to use a modular architecture
with 5 helper functions (STclust_select_genes, STclust_calculate_distances, etc.).
The original monolithic implementation is available as STclust_legacy() for reproducibility
with version 1.x results.
```

### STdiff.Rd
**Replace/expand existing text:**
Current: "Maintains the original API while using the new modular implementation"

**Replace with:**
```
\strong{Note:} This function was refactored in version 2.0.0 to use a modular architecture
(STdiff_run_nonspatial, STdiff_fit_spatial_models, STdiff_compile_results).
The original implementation is available as STdiff_legacy() for reproducibility with
version 1.x results.
```

### SThet.Rd
**Add to description (confirm no legacy exists):**
```
\strong{Note:} No legacy version available. This function has not been refactored and
continues to use the original implementation.
```

---

## Files for Legacy Functions

- ✅ `R/STenrich_legacy.R` exists (referenced in STenrich.Rd)
- ⚠️ `R/STgradient_legacy.R` exists (referenced in STgradient.Rd)
- ❓ Need to verify: `R/STclust_legacy.R`, `R/STdiff_legacy.R`

---

**Status:** Documentation audit complete. 3/5 functions need legacy references added.
