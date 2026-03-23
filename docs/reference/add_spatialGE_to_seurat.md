<div id="main" class="col-md-9" role="main">

# Add spatialGE results to Seurat object

<div class="ref-description section level2">

**\[stable\]**

Adds spatial analysis results (e.g., cluster assignments, Moran's I) to
a Seurat object's metadata.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
add_spatialGE_to_seurat(
  seurat_obj,
  st_obj,
  result_type = "clusters",
  custom.metadata = NULL,
  verbose = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   seurat_obj:

    original Seurat object

-   st_obj:

    STlist object with analysis results

-   result_type:

    type of results to add: 'clusters', 'moran', 'geary', 'custom'

-   custom.metadata:

    named list of metadata columns to add

-   verbose:

    logical, whether to print progress

</div>

<div class="section level2">

## Value

Seurat object with added metadata

</div>

<div class="section level2">

## Details

After running spatialGE analysis on an STlist, this function transfers
the results back to the original Seurat object for integrated
visualization and downstream analysis.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`as.STlist.Seurat`, `as.Seurat.STlist`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

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

</div>

</div>

</div>
