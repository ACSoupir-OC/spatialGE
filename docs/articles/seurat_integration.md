<div id="main" class="col-md-9" role="main">

# Seurat Integration: Bidirectional Workflows

<div class="section level2">

## Overview

This vignette demonstrates how to integrate **spatialGE** with
**Seurat** for seamless spatial transcriptomics workflows. The
`spatialGE` package provides bidirectional conversion utilities that
allow you to:

1.  **Convert Seurat objects to STlist** - Run spatialGE analyses on
    Seurat-processed data
2.  **Convert STlist to Seurat** - Leverage Seurat’s visualization and
    downstream analysis tools
3.  **Run spatialGE directly from Seurat** - One-step analysis workflows
4.  **Transfer results back to Seurat** - Add cluster assignments,
    spatial statistics to Seurat metadata

</div>

<div class="section level2">

## Installation

Make sure both packages are installed:

<div id="cb1" class="sourceCode">

``` r
# Install spatialGE
devtools::install_github("ACSoupir-OC/spatialGE")

# Install Seurat (if not already installed)
install.packages("Seurat")
```

</div>

Load the required libraries:

<div id="cb2" class="sourceCode">

``` r
library(spatialGE)
library(Seurat)
library(Matrix)
library(tibble)
```

</div>

</div>

<div class="section level2">

## Creating Example Data

For this vignette, we’ll create a synthetic STlist object to demonstrate
the workflow:

<div id="cb3" class="sourceCode">

``` r
set.seed(42)

# Create synthetic count matrix (20 genes x 100 cells)
n_genes <- 20
n_cells <- 100

counts_mat <- Matrix::Matrix(
  matrix(rnbinom(n_genes * n_cells, mu = 5, size = 2), 
         nrow = n_genes, ncol = n_cells), 
  sparse = TRUE
)
rownames(counts_mat) <- paste0("Gene", 1:n_genes)
colnames(counts_mat) <- paste0("Cell", 1:n_cells)

# Add some marker gene names
rownames(counts_mat)[1:5] <- c("MLANA", "CD37", "TP53", "GAPDH", "ACTB")

# Create spatial coordinates
coords_df <- data.frame(
  barcode = paste0("Cell", 1:n_cells),
  xpos = runif(n_cells, 0, 1000),
  ypos = runif(n_cells, 0, 1000)
)

# Create STlist object
st_obj <- methods::new("STlist",
  counts = list(test_sample = counts_mat),
  spatial_meta = list(test_sample = coords_df),
  sample_meta = tibble::tibble(sample_name = "test_sample"),
  gene_meta = list()
)

cat(sprintf("Created STlist: %d genes x %d cells\n", n_genes, n_cells))
#> Created STlist: 20 genes x 100 cells
```

</div>

</div>

<div class="section level2">

## Workflow 1: STlist to Seurat Conversion

Convert your STlist object to Seurat for visualization:

<div id="cb4" class="sourceCode">

``` r
# Convert STlist to Seurat object
seurat_obj <- as.Seurat.STlist(st_obj, verbose = FALSE)

# Inspect the Seurat object
cat(sprintf("Seurat object:\n"))
#> Seurat object:
cat(sprintf("  Cells: %d\n", ncol(seurat_obj)))
#>   Cells: 100
cat(sprintf("  Features: %d\n", nrow(seurat_obj)))
#>   Features: 20

# View metadata
head(seurat_obj@meta.data)
#>        orig.ident nCount_RNA nFeature_RNA         x        y barcode      xpos
#> Cell1 test_sample        148           19 948.00454 445.0047   Cell1 948.00454
#> Cell2 test_sample         53           16 197.91424 534.5598   Cell2 197.91424
#> Cell3 test_sample         76           17  89.54941 923.2419   Cell3  89.54941
#> Cell4 test_sample         73           16  67.94101 429.1973   Cell4  67.94101
#> Cell5 test_sample        116           18 673.35430 847.3947   Cell5 673.35430
#> Cell6 test_sample        103           18 483.77301 846.9969   Cell6 483.77301
#>           ypos
#> Cell1 445.0047
#> Cell2 534.5598
#> Cell3 923.2419
#> Cell4 429.1973
#> Cell5 847.3947
#> Cell6 846.9969
```

</div>

</div>

<div class="section level2">

## Workflow 2: Seurat to STlist Conversion

Convert back from Seurat to STlist (round-trip validation):

<div id="cb5" class="sourceCode">

``` r
# Convert Seurat back to STlist
st_obj2 <- as.STlist.Seurat(seurat_obj, verbose = FALSE)

cat(sprintf("Round-trip validation:\n"))
#> Round-trip validation:
cat(sprintf("  Original: %d genes x %d cells\n", 
            nrow(st_obj@counts[[1]]), ncol(st_obj@counts[[1]])))
#>   Original: 20 genes x 100 cells
cat(sprintf("  After round-trip: %d genes x %d cells\n", 
            nrow(st_obj2@counts[[1]]), ncol(st_obj2@counts[[1]])))
#>   After round-trip: 20 genes x 100 cells
cat(sprintf("  Match: %s\n", 
            all(dim(st_obj@counts[[1]]) == dim(st_obj2@counts[[1]]))))
#>   Match: TRUE
```

</div>

</div>

<div class="section level2">

## Workflow 3: Add Results to Seurat Metadata

You can add custom analysis results to Seurat metadata:

<div id="cb6" class="sourceCode">

``` r
# Add cluster assignments
set.seed(42)
seurat_obj$spatialGE_cluster <- sample(c("A", "B", "C"), ncol(seurat_obj), replace = TRUE)

# Add spatial statistics
seurat_obj$moran_I <- runif(ncol(seurat_obj), 0, 0.5)

cat("Metadata columns added:\n")
#> Metadata columns added:
head(seurat_obj@meta.data[, c("orig.ident", "spatialGE_cluster", "moran_I")])
#>        orig.ident spatialGE_cluster   moran_I
#> Cell1 test_sample                 A 0.2899104
#> Cell2 test_sample                 A 0.4107020
#> Cell3 test_sample                 A 0.0568593
#> Cell4 test_sample                 A 0.3822539
#> Cell5 test_sample                 B 0.3118067
#> Cell6 test_sample                 B 0.0742233
```

</div>

</div>

<div class="section level2">

## Workflow 4: Gene Set Integration

Convert spatialGE enrichment results for Seurat GSEA or AddModuleScore:

<div id="cb7" class="sourceCode">

``` r
# Define pathway gene sets
pathway_list <- list(
  "Pathway1" = c("MLANA", "CD37", "TP53"),
  "Pathway2" = c("GAPDH", "ACTB")
)

# Convert to Seurat format
gene_sets <- spatialGE_to_seurat_genesets(
  list(test_sample = data.frame(
    geneset = names(pathway_list),
    genes = I(pathway_list),
    pval = c(0.01, 0.05),
    stringsAsFactors = FALSE
  )),
  pval.thr = 0.05,
  return.type = "list"
)

cat(sprintf("Generated %d gene sets\n", length(gene_sets)))
#> Generated 1 gene sets
for (i in seq_along(gene_sets)) {
  cat(sprintf("  - %s: %d genes\n", names(gene_sets)[i], length(gene_sets[[i]])))
}
#>   - test_sample_Pathway1: 1 genes
```

</div>

</div>

<div class="section level2">

## Complete Example: End-to-End Workflow

Here’s a complete workflow from start to finish:

<div id="cb8" class="sourceCode">

``` r
# ============================================================
# Complete workflow from STlist to Seurat and back
# ============================================================

# Step 1: Create or load STlist object
# st_obj <- STlist(rnacounts = count_files, spotcoords = coord_files)

# Step 2: Run spatialGE analysis
st_obj <- SThet(st_obj, genes = c("MLANA", "CD37", "TP53"))
st_obj <- STclust(st_obj, k = 3)

# Step 3: Convert to Seurat for visualization
seurat_obj <- as.Seurat.STlist(st_obj, verbose = TRUE)

# Step 4: Add results to Seurat metadata
seurat_obj <- add_spatialGE_to_seurat(seurat_obj, st_obj, result_type = "clusters")

# Step 5: Visualize in Seurat
# SpatialDimPlot(seurat_obj, group.by = "spatialGE_cluster")
# FeaturePlot(seurat_obj, features = c("MLANA", "CD37"))

# Step 6: Save results
# saveRDS(seurat_obj, "seurat_with_spatialge_results.rds")
```

</div>

</div>

<div class="section level2">

## Troubleshooting

<div class="section level3">

### Issue: No spatial coordinates found

Check if your object has spatial information:

<div id="cb9" class="sourceCode">

``` r
# For Seurat objects, check images slot
names(seurat_obj@images)

# For STlist objects, check spatial_meta
names(st_obj@spatial_meta)
```

</div>

</div>

<div class="section level3">

### Issue: Memory errors with large objects

Process samples individually:

<div id="cb10" class="sourceCode">

``` r
# Process one sample at a time
for (sample in names(st_obj@counts)) {
  st_single <- st_obj[, sample]
  st_single <- SThet(st_single, genes = genes_of_interest)
}
```

</div>

</div>

</div>

<div class="section level2">

## Best Practices

1.  **Always validate after conversion**

    <div id="cb11" class="sourceCode">

    ``` r
    st_obj <- as.STlist.Seurat(seurat_obj)
    print(st_obj)
    summary(st_obj)
    ```

    </div>

2.  **Preserve original objects**

    <div id="cb12" class="sourceCode">

    ``` r
    seurat_original <- seurat_obj  # Keep backup
    seurat_obj <- as.Seurat.STlist(st_obj)  # Overwrite carefully
    ```

    </div>

3.  **Use raw counts for spatialGE**

    <div id="cb13" class="sourceCode">

    ``` r
    # Check available assays
    Seurat::Assays(seurat_obj)

    # Use raw counts
    st_obj <- as.STlist.Seurat(seurat_obj, assay = "RNA", slot = "counts")
    ```

    </div>

</div>

<div class="section level2">

## See Also

-   `vignette("basic_functions_vignette", package = "spatialGE")` - Core
    spatialGE functions
-   Seurat website: <https://satijalab.org/seurat/>

</div>

<div class="section level2">

## References

-   Hao Y, et al. (2024). “Seurat v5.” *Nature Biotechnology*.
-   Stuart T, et al. (2019). “Single-Cell Data Integration.” *Cell*.

</div>

</div>
