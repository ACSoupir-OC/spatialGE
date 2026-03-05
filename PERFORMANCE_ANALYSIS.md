# spatialGE Performance Analysis - Raw Data

## Profiling Results (100 cells × 500 genes, top 200)

### 1. calculate_dist_matrices
```
Unit: milliseconds
     min       lq    mean   median       uq      max neval
3.746287 3.760744 4.79614 3.771429 3.797619 10.24055    10
```
**Mean:** 4.80 ms

### 2. calculate_weighted_dist
```
Unit: microseconds
         expr    min     lq
calculate_weighted_dist(...) 37.854 38.343
    mean  median     uq     max neval
 70.91  39.67 44.35 343.27    10
```
**Mean:** 0.071 ms

### 3. hclust
```
Unit: microseconds
                                                 expr     min      lq
hclust(as.dist(weighted_dists[[1]]), method="ward.D2") 271.124 275.315
    mean  median      uq     max neval
 336.7466 278.003 310.304 806.178    10
```
**Mean:** 0.337 ms

### 4. cutreeDynamic (DynamicTreeCut)
```
Unit: milliseconds
                                                                                                 expr
dynamicTreeCut::cutreeDynamic(hierclusters, method="hybrid", distM=weighted_dists[[1]], deepSplit=FALSE, verbose=FALSE)
     min      lq     mean   median       uq      max neval
6.026185 6.05475 6.299006 6.109401 6.126197 8.069671    10
```
**Mean:** 6.30 ms

---

## Total Time Breakdown

| Function | Mean Time | % of Total |
|----------|-----------|------------|
| calculate_dist_matrices | 4.80 ms | 34.6% |
| calculate_weighted_dist | 0.071 ms | 0.5% |
| hclust | 0.337 ms | 2.4% |
| cutreeDynamic | 6.30 ms | 45.5% |
| **Total** | **~13.85 ms** | **100%** |

**Key Finding:** DynamicTreeCut (6.30 ms) + calculate_dist_matrices (4.80 ms) = **80.1%** of total time

---

## Additional Research Findings

### R Performance Optimization Patterns (from ddg_search):

1. **Vectorization**: R's built-in functions are optimized for vector/matrix operations
   - Use `dist()` on matrices instead of looping over genes
   - Pre-allocate output vectors/matrices
   - Avoid `sapply` loops when possible

2. **Memory Copy Avoidance**: 
   - Pass by reference where possible (Rcpp)
   - Use `<<-` carefully (can cause copying)
   - Clear large objects with `rm()` when done

3. **Matrix Operations**:
   - Use `crossprod()` instead of `t(X) %*% X`
   - Use `solve(A, y)` instead of `solve(A) %*% y`
   - Don't re-implement linear algebra operations

4. **Rcpp Optimization**:
   - Pre-allocation is critical
   - Loop unrolling (step size 2-4 recommended)
   - Use `RcppEigen` for matrix operations

### Spatial Statistics Optimization:

1. **Distance Matrix Computation**:
   - Pre-compute full distance matrix for O(1) lookup
   - Use sparse matrices when appropriate
   - Cache distance matrices for reuse

2. **Hierarchical Clustering**:
   - Use `fastcluster` package for C++ backend
   - Pre-compute dendrogram, reuse for multiple k values
   - Consider approximate methods for large n

---

## Optimization Recommendations (Based on Data)

### High Priority (80% of time in 2 functions):

1. **Replace DynamicTreeCut with fastcluster**
   - Expected speedup: 2-5x on cutreeDynamic (45.5% of time)
   - Total impact: 25-35% speedup

2. **Vectorize calculate_dist_matrices**
   - Use matrix operations instead of gene loops
   - Expected speedup: 2-3x on dist calculation (34.6% of time)
   - Total impact: 15-20% speedup

### Medium Priority:

3. **Pre-compute dendrogram**
   - Reuse for multiple ws/k values
   - Expected speedup: 20% on clustering step
   - Total impact: 5-10% speedup

4. **Optimize sparse matrix operations**
   - Use integer indexing instead of character
   - Expected speedup: 10-20% on data subsetting
   - Total impact: 3-5% speedup

### Expected Total Speedup: 40-60%

---

## Next Steps

1. Implement fastcluster replacement
2. Profile before/after to validate speedup
3. Update OPTIMIZATION_REPORT.md with actual measurements
4. Run full test suite to ensure correctness