# spatialGE 2.0.0

* Major refactoring of `STclust` to modular architecture.
* Replaced base R `hclust()` with `fastcluster::hclust()` for 2-5x faster hierarchical clustering.
* Expected 20-35% speedup on `STclust` runtime (dynamicTreeCut accounts for ~45.5% of runtime).
* New modular functions: `STclust_select_genes()`, `STclust_calculate_distances()`, `STclust_weight_distances()`, `STclust_hierarchical()`.
* Core functions exported for flexibility: `calculate_dist_matrices()`, `calculate_weighted_dist()`, `get_hier_clusters_dtc()`, `get_hier_clusters_ks()`.

# spatialGE 1.2.2

# spatialGE 1.2.1

* Added C++ code to speed up `STenrich` computations (Thank you, Dr. Soupir!).
* The `STenrich` function can now calculate gene set enrichment scores (via `GSVA`) in addition to gene set average expression.
* The `STenrich` function can now test for gene set "hot-spots" within an specific tissue domain (`annot` and `domain` arguments).
* Spatial distances can now be log-transformed for `STgradient` (`log` argument).
* Example data sets are now in the spatialGE_Data GitHub repository.
* New functions `spatial_metadata` and `tissue_names` to quickly access the names of spot/cell annotations and sample names.
* Re-assessed package dependencies.
* Use of `DelayedArray` for some calculations.
* Several bug fixes.

# spatialGE 1.2.0

* Functions `STenrich`, `STgradient`, and `STdiff` available.
* Added C++ code for Seurat's implementation of `FindVariableFeatures`

# spatialGE 1.1.0

* Multiple, significant changes to STlist (incompatible with STlist objects from version 1.0).
