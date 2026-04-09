# Legacy Function Export Verification

**Date:** 2026-04-09  
**Package:** spatialGE

---

## Summary Table

| Function | Exported | Documented | Status |
|----------|----------|------------|--------|
| `STclust_legacy` | ✅ Y | ✅ Y | OK |
| `STdiff_legacy` | ✅ Y | ✅ Y | OK |
| `STenrich_legacy` | ✅ Y | ✅ Y | OK |
| `STgradient_legacy` | ✅ Y | ✅ Y | OK |
| `SThet_legacy` | ✅ Y | ✅ Y | OK |

---

## Details

All five legacy functions are **properly exported** in `NAMESPACE` and have corresponding `.Rd` documentation files in `man/`:

| Function | NAMESPACE Export | man/ File |
|----------|-----------------|-----------|
| `STclust_legacy` | `export(STclust_legacy)` | `STclust_legacy.Rd` |
| `STdiff_legacy` | `export(STdiff_legacy)` | `STdiff_legacy.Rd` |
| `STenrich_legacy` | `export(STenrich_legacy)` | `STenrich_legacy.Rd` |
| `STgradient_legacy` | `export(STgradient_legacy)` | `STgradient_legacy.Rd` |
| `SThet_legacy` | `export(SThet_legacy)` | `SThet_legacy.Rd` |

---

## Additional Legacy Exports Found

| Function | Exported | Documented | Status |
|----------|----------|------------|--------|
| `STList_legacy` | ✅ Y | ✅ Y | OK |
| `SThet_invdist_test_legacy` | ✅ Y | ✅ Y | OK |

---

## Recommendations

✅ **No action required** - All legacy functions are properly exported and documented. The package maintains full backward compatibility for legacy API users.

**Total legacy functions:** 7 (5 core + 2 additional)  
**Exported:** 7/7 (100%)  
**Documented:** 7/7 (100%)

---

## Verification Method

1. Parsed `NAMESPACE` file for `export()` statements matching `*_legacy` pattern
2. Cross-referenced with `man/*.Rd` files for documentation coverage
3. No missing exports or undocumented functions found
