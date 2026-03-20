# Seurat Integration: Bidirectional Workflows

## Overview

This vignette demonstrates how to integrate **spatialGE** with
**Seurat** for seamless spatial transcriptomics workflows. The
`spatialGE` package now provides bidirectional conversion utilities that
allow you to:

1.  **Convert Seurat objects to STlist** - Run spatialGE analyses on
    Seurat-processed data
2.  **Convert STlist to Seurat** - Leverage Seurat’s visualization and
    downstream analysis tools
3.  **Run spatialGE directly from Seurat** - One-step analysis workflows
4.  **Transfer results back to Seurat** - Add cluster assignments,
    spatial statistics to Seurat metadata

## Installation

Make sure both packages are installed:

``` r

# Install spatialGE
# install.packages("spatialGE")  # CRAN version
# devtools::install_github("ACSoupir-OC/spatialGE")  # Development version

# Install Seurat (if not already installed)
install.packages("Seurat")
```

Load the required libraries:

``` r

library(spatialGE)

# Seurat is required for the examples in this vignette
# Install with: install.packages("Seurat")
if (requireNamespace("Seurat", quietly = TRUE)) {
  library(Seurat)
  has_seurat <- TRUE
} else {
  has_seurat <- FALSE
  message("Seurat not installed. Examples will be skipped.")
}
```

## Workflow 1: Seurat → STlist → spatialGE Analysis

### Starting with a Seurat Object

Assume you have a Seurat object that has been processed through the
standard Seurat workflow:

``` r

# Example: Create or load a Seurat object
# seurat_obj <- readRDS("my_seurat_object.rds")

# Or create from 10X Visium data
# seurat_obj <- Load10X_Spatial(data.dir = "path/to/visium/data")

# Standard Seurat preprocessing
# seurat_obj <- SCTransform(seurat_obj)
# seurat_obj <- RunPCA(seurat_obj)
# seurat_obj <- FindNeighbors(seurat_obj)
# seurat_obj <- FindClusters(seurat_obj)
```

For this vignette, we’ll show code examples. To run them, you’ll
need: 1. spatialGE installed: `install.packages("spatialGE")` 2. Seurat
installed: `install.packages("Seurat")` 3. Example data (see code below)

``` r

# Load example data that comes with spatialGE
data("melanoma_thrane", package = "spatialGE")

# melanoma_thrane is an STlist object
# Let's see what's in it
print(melanoma_thrane)
summary(melanoma_thrane)

# Or load your own STlist object
# st_obj <- readRDS("my_spatial_data.rds")
```

### Convert Seurat to STlist

Use
[`as.STlist.Seurat()`](https://acsoupir-oc.github.io/spatialGE/reference/as.STlist.Seurat.md)
to convert:

``` r

# Convert Seurat object to STlist
st_obj <- as.STlist.Seurat(
  seurat_obj,
  assay = "RNA",      # Which assay to use
  slot = "counts",    # Which slot (counts, data, scale.data)
  use.spatial = TRUE, # Extract spatial coordinates
  verbose = TRUE
)

# Inspect the STlist object
print(st_obj)
summary(st_obj)
```

**Key Parameters:**

| Parameter     | Description                   | Default    |
|---------------|-------------------------------|------------|
| `assay`       | Which Seurat assay to extract | `"RNA"`    |
| `slot`        | Which slot within the assay   | `"counts"` |
| `use.spatial` | Extract spatial coordinates   | `TRUE`     |
| `verbose`     | Print progress messages       | `TRUE`     |

### Run spatialGE Analysis

Now run any spatialGE analysis on the converted object:

``` r

# Run spatial heterogeneity analysis
genes_of_interest <- c("MLANA", "CD37", "TP53")
st_obj <- SThet(st_obj, genes = genes_of_interest)

# View results
print(st_obj@gene_meta[[1]])

# Run spatial clustering
st_obj <- STclust(st_obj)

# View cluster assignments
head(st_obj@spatial_meta[[1]])
```

## Workflow 2: STlist → Seurat → Visualization

### Convert STlist to Seurat

After running spatialGE analysis, convert back to Seurat for
visualization:

``` r

# Convert STlist to Seurat object
seurat_obj <- as.Seurat.STlist(
  st_obj,
  samples = NULL,           # All samples, or specify c("sample1", "sample2")
  assay.name = "RNA",       # Name for the assay
  add.spatial.info = TRUE,  # Include spatial coordinates
  verbose = TRUE
)

# Inspect the Seurat object
print(seurat_obj)
```

### Visualize in Seurat

Leverage Seurat’s powerful visualization tools:

``` r

# FeaturePlot for gene expression
FeaturePlot(seurat_obj, features = c("MLANA", "CD37"))

# DimPlot for clustering (if clusters were added)
DimPlot(seurat_obj, group.by = "spatialGE_cluster")

# Spatial plotting (for Visium data)
SpatialDimPlot(seurat_obj, group.by = "spatialGE_cluster")
SpatialFeaturePlot(seurat_obj, features = "MLANA")
```

## Workflow 3: Direct Analysis with `spatialGE_from_seurat()`

For quick analyses, use the wrapper function:

``` r

# One-step heterogeneity analysis
result <- spatialGE_from_seurat(
  seurat_obj,
  analysis = "SThet",
  genes = c("MLANA", "CD37", "TP53"),
  assay = "RNA",
  verbose = TRUE
)

# One-step clustering
result <- spatialGE_from_seurat(
  seurat_obj,
  analysis = "STclust",
  assay = "RNA"
)

# One-step differential expression
result <- spatialGE_from_seurat(
  seurat_obj,
  analysis = "STdiff",
  annot = "cell_type"  # Metadata column for grouping
)
```

**Available Analyses:**

- `"SThet"` - Spatial heterogeneity (Moran’s I, Geary’s C)
- `"STclust"` - Spatial domain detection
- `"STdiff"` - Spatial differential expression
- `"STenrich"` - Spatial enrichment analysis
- `"STgradient"` - Expression gradient analysis

## Workflow 4: Adding Results Back to Seurat

Use
[`add_spatialGE_to_seurat()`](https://acsoupir-oc.github.io/spatialGE/reference/add_spatialGE_to_seurat.md)
to transfer analysis results:

``` r

# After running SThet
st_obj <- SThet(st_obj, genes = c("MLANA", "CD37"))

# Add Moran's I values to Seurat feature metadata
seurat_obj <- add_spatialGE_to_seurat(
  seurat_obj,
  st_obj,
  result_type = "moran",  # or "geary"
  verbose = TRUE
)

# View in Seurat
head(seurat_obj@assays$RNA@meta.features)

# After running STclust
st_obj <- STclust(st_obj)

# Add cluster assignments to Seurat cell metadata
seurat_obj <- add_spatialGE_to_seurat(
  seurat_obj,
  st_obj,
  result_type = "clusters",
  verbose = TRUE
)

# View cluster assignments
head(seurat_obj$spatialGE_cluster)
```

**Result Types:**

| Type | Description | Location in Seurat |
|----|----|----|
| `"clusters"` | Spatial domain assignments | `seurat_obj$spatialGE_cluster` |
| `"moran"` | Moran’s I statistics | `seurat_obj@assays$RNA@meta.features` |
| `"geary"` | Geary’s C statistics | `seurat_obj@assays$RNA@meta.features` |
| `"custom"` | Custom metadata | User-specified |

## Workflow 5: Gene Set Analysis Integration

Convert spatialGE enrichment results for Seurat GSEA:

``` r

# Run enrichment analysis
pathway_list <- list(
  "EMT" = c("VIM", "CDH2", "FN1", "SNAI1"),
  "CellCycle" = c("MKI67", "TOP2A", "PCNA")
)

enrich_result <- STenrich(st_obj, gene_sets = pathway_list)

# Convert to Seurat gene sets
gene_sets <- spatialGE_to_seurat_genesets(
  enrich_result,
  pval.thr = 0.05,
  return.type = "list"
)

# Add module scores to Seurat
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = gene_sets,
  name = "spatialGE_"
)

# View module scores
head(seurat_obj$spatialGE_1)
```

## Complete Example: End-to-End Workflow

Here’s a complete workflow from start to finish:

``` r

# ============================================================
# Step 1: Start with Seurat object (already preprocessed)
# ============================================================
# seurat_obj <- Load10X_Spatial(data.dir = "visium_data/")
# seurat_obj <- SCTransform(seurat_obj)

# ============================================================
# Step 2: Convert to STlist
# ============================================================
st_obj <- as.STlist.Seurat(seurat_obj, verbose = TRUE)

# ============================================================
# Step 3: Run spatial heterogeneity analysis
# ============================================================
genes <- c("MLANA", "CD37", "TP53", "GAPDH", "ACTB")
st_obj <- SThet(st_obj, genes = genes)

# ============================================================
# Step 4: Run spatial clustering
# ============================================================
st_obj <- STclust(st_obj, k = 5)

# ============================================================
# Step 5: Convert back to Seurat
# ============================================================
seurat_obj <- as.Seurat.STlist(st_obj, verbose = TRUE)

# ============================================================
# Step 6: Add results to Seurat metadata
# ============================================================
seurat_obj <- add_spatialGE_to_seurat(
  seurat_obj, st_obj,
  result_type = "clusters"
)
seurat_obj <- add_spatialGE_to_seurat(
  seurat_obj, st_obj,
  result_type = "moran"
)

# ============================================================
# Step 7: Visualize in Seurat
# ============================================================
# SpatialDimPlot(seurat_obj, group.by = "spatialGE_cluster")
# FeaturePlot(seurat_obj, features = c("MLANA", "CD37"))

# ============================================================
# Step 8: Save results
# ============================================================
# saveRDS(seurat_obj, "seurat_with_spatialge_results.rds")
```

## Handling Multi-Sample Objects

### Merged Seurat Objects

When working with merged Seurat objects (multiple samples):

``` r

# Seurat object with multiple samples
# seurat_merged <- merge(seurat1, seurat2, seurat3)

# Convert to STlist (preserves sample structure)
st_obj <- as.STlist.Seurat(seurat_merged)

# Check samples
names(st_obj@counts)  # Sample names preserved

# Run analysis on all samples
st_obj <- SThet(st_obj, genes = c("MLANA", "CD37"))

# Convert back (will merge automatically)
seurat_merged <- as.Seurat.STlist(st_obj)

# Or specify specific samples
seurat_sample1 <- as.Seurat.STlist(st_obj, samples = "sample_1")
```

### Batch Processing

``` r

# Process samples individually
results_list <- list()

for (sample_name in names(st_obj@counts)) {
  # Single sample analysis
  st_single <- st_obj[, sample_name]  # Subset if supported
  st_single <- SThet(st_single, genes = genes_of_interest)
  results_list[[sample_name]] <- st_single@gene_meta
}

# Combine results
combined_results <- do.call(rbind, results_list)
```

## Troubleshooting

### Issue: “Seurat package is required”

``` r

# Install Seurat
install.packages("Seurat")

# Or development version
# remotes::install_github("satijalab/seurat")
```

### Issue: No spatial coordinates found

``` r

# Check if Seurat object has images
names(seurat_obj@images)

# If empty, coordinates won't be extracted
# Solution: Load spatial data properly with Load10X_Spatial()

# Or add coordinates manually after conversion
st_obj@spatial_meta[[1]] <- your_coordinates_df
```

### Issue: Cell barcode mismatch

When converting between formats, cell barcodes may not match exactly:

``` r

# Check barcodes
head(colnames(seurat_obj))
head(st_obj@spatial_meta[[1]]$barcode)

# The conversion functions handle this automatically
# by matching common barcodes

# If issues persist, check for sample prefixes
# Seurat merged objects: SampleID_CellBarcode
# STlist: CellBarcode only
```

### Issue: Memory errors with large objects

``` r

# Process samples individually instead of merged
for (sample in names(seurat_obj@images)) {
  seurat_single <- subset(seurat_obj, images = sample)
  st_single <- as.STlist.Seurat(seurat_single)
  # ... analysis ...
}

# Or use sparse matrices (default in spatialGE)
# Check memory usage
pryr::mem_used()
```

## Best Practices

1.  **Always validate after conversion**

    ``` r

    st_obj <- as.STlist.Seurat(seurat_obj)
    print(st_obj)
    summary(st_obj)
    ```

2.  **Preserve original objects**

    ``` r

    seurat_original <- seurat_obj  # Keep backup
    seurat_obj <- as.Seurat.STlist(st_obj)  # Overwrite carefully
    ```

3.  **Check assay compatibility**

    ``` r

    # Use raw counts for spatialGE
    Seurat::Assays(seurat_obj)  # Check available assays
    st_obj <- as.STlist.Seurat(seurat_obj, assay = "RNA", slot = "counts")
    ```

4.  **Document your workflow**

    ``` r

    # Save conversion metadata
    st_obj@misc$conversion_info <- list(
      source = "seurat",
      assay = "RNA",
      timestamp = Sys.time(),
      seurat_version = packageVersion("Seurat")
    )
    ```

## See Also

- `vignette("basic_functions_vignette", package = "spatialGE")` - Core
  spatialGE functions
- `vignette("spatial_differential_expression", package = "spatialGE")` -
  STdiff tutorial
- `vignette("spatial_enrichment_gradients_smi", package = "spatialGE")` -
  STenrich and STgradient
- Seurat website: <https://satijalab.org/seurat/>

## References

- Hao Y, et al. (2024). “Seurat v5: A comprehensive toolkit for
  single-cell and spatial data analysis.” *Nature Biotechnology*.
- Stuart T, et al. (2019). “Comprehensive Integration of Single-Cell
  Data.” *Cell*.
- Soupir AC, et al. (2024). “spatialGE: Spatial transcriptomics analysis
  tools.” *bioRxiv*.
