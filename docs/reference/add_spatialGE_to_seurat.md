# Add spatialGE results to Seurat object

**\[stable\]**

Adds spatial analysis results (e.g., cluster assignments, Moran's I) to
a Seurat object's metadata.

## Usage

``` r
add_spatialGE_to_seurat(
  seurat_obj,
  st_obj,
  result_type = "clusters",
  custom.metadata = NULL,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  original Seurat object

- st_obj:

  STlist object with analysis results

- result_type:

  type of results to add: 'clusters', 'moran', 'geary', 'custom'

- custom.metadata:

  named list of metadata columns to add

- verbose:

  logical, whether to print progress

## Value

Seurat object with added metadata

## Details

After running spatialGE analysis on an STlist, this function transfers
the results back to the original Seurat object for integrated
visualization and downstream analysis.

## See also

[`as.STlist.Seurat`](https://acsoupir-oc.github.io/spatialGE/reference/as.STlist.Seurat.md),
[`as.Seurat.STlist`](https://acsoupir-oc.github.io/spatialGE/reference/as.Seurat.STlist.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert to STlist and run clustering
st_obj <- as.STlist.Seurat(seurat_obj)
st_obj <- STclust(st_obj)

# Add cluster assignments back to Seurat
seurat_obj <- add_spatialGE_to_seurat(seurat_obj, st_obj, result_type = 'clusters')

# Visualize in Seurat
Seurat::DimPlot(seurat_obj, group.by = 'spatialGE_cluster')
} # }
```
