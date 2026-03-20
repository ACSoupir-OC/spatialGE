# Changelog

## spatialGE 2.0.0

- Major refactoring of `STclust` to modular architecture.
- Replaced base R [`hclust()`](https://rdrr.io/r/stats/hclust.html) with
  [`fastcluster::hclust()`](https://rdrr.io/pkg/fastcluster/man/hclust.html)
  for 2-5x faster hierarchical clustering.
- Expected 20-35% speedup on `STclust` runtime (dynamicTreeCut accounts
  for ~45.5% of runtime).
- New modular functions:
  [`STclust_select_genes()`](https://acsoupir-oc.github.io/spatialGE/reference/STclust_select_genes.md),
  [`STclust_calculate_distances()`](https://acsoupir-oc.github.io/spatialGE/reference/STclust_calculate_distances.md),
  [`STclust_weight_distances()`](https://acsoupir-oc.github.io/spatialGE/reference/STclust_weight_distances.md),
  [`STclust_hierarchical()`](https://acsoupir-oc.github.io/spatialGE/reference/STclust_hierarchical.md).
- Core functions exported for flexibility:
  [`calculate_dist_matrices()`](https://acsoupir-oc.github.io/spatialGE/reference/calculate_dist_matrices.md),
  [`calculate_weighted_dist()`](https://acsoupir-oc.github.io/spatialGE/reference/calculate_weighted_dist.md),
  [`get_hier_clusters_dtc()`](https://acsoupir-oc.github.io/spatialGE/reference/get_hier_clusters_dtc.md),
  [`get_hier_clusters_ks()`](https://acsoupir-oc.github.io/spatialGE/reference/get_hier_clusters_ks.md).

## spatialGE 1.2.2

CRAN release: 2025-06-04

## spatialGE 1.2.1

CRAN release: 2025-05-26

- Added C++ code to speed up `STenrich` computations (Thank you,
  Dr. Soupir!).
- The `STenrich` function can now calculate gene set enrichment scores
  (via `GSVA`) in addition to gene set average expression.
- The `STenrich` function can now test for gene set “hot-spots” within
  an specific tissue domain (`annot` and `domain` arguments).
- Spatial distances can now be log-transformed for `STgradient` (`log`
  argument).
- Example data sets are now in the spatialGE_Data GitHub repository.
- New functions `spatial_metadata` and `tissue_names` to quickly access
  the names of spot/cell annotations and sample names.
- Re-assessed package dependencies.
- Use of `DelayedArray` for some calculations.
- Several bug fixes.

## spatialGE 1.2.0

CRAN release: 2025-05-13

- Functions `STenrich`, `STgradient`, and `STdiff` available.
- Added C++ code for Seurat’s implementation of `FindVariableFeatures`

## spatialGE 1.1.0

- Multiple, significant changes to STlist (incompatible with STlist
  objects from version 1.0).
