# spatialGE 2.0.0 (2026-03-22)

## Major Changes

* **Comprehensive Documentation** - All 7 core functions now include @examples with TNBC test data
  - STlist: Data loading from Visium, CosMx, CSV/TSV
  - STplot: Gene expression and cluster visualization
  - SThet: Moran's I spatial heterogeneity analysis
  - STclust: Tissue domain detection with DynamicTreeCut
  - STdiff: Differential expression (Wilcoxon, t-test, spatial mixed models)
  - STenrich: Gene set enrichment with permutation testing
  - STgradient: Spatial gradient detection with Spearman correlation

* **Modular Architecture** - Major refactoring of core functions
  - `STclust` modular architecture with 5 helper functions
  - `STdiff` modular architecture (non-spatial + spatial tests)
  - Replaced base R `hclust()` with `fastcluster::hclust()` for 2-5x speedup
  - Expected 20-35% speedup on `STclust` runtime

* **New Modular Functions**:
  - `STclust_select_genes()`, `STclust_calculate_distances()`, `STclust_weight_distances()`
  - `STclust_hierarchical()`, `STclust_get_hier_clusters_*()`
  - `STdiff_run_nonspatial()`, `STdiff_fit_spatial_models()`, `STdiff_compile_results()`
  - Core functions exported for flexibility

* **C++ Optimizations**:
  - C++ code for Seurat's `FindVariableFeatures` implementation
  - C++ code for `STenrich` computations (GSVA integration)

## New Features

* `STenrich` can now calculate gene set enrichment scores via GSVA
* `STenrich` can test for gene set "hot-spots" within specific tissue domains
* `STgradient` supports log-transformed spatial distances
* New functions: `spatial_metadata()`, `tissue_names()` for quick metadata access
* Example datasets now in spatialGE_Data repository
* Use of `DelayedArray` for memory-efficient calculations

## Bug Fixes & Improvements

* Re-assessed package dependencies
* Multiple bug fixes in STlist object handling
* Improved error messages and validation

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
