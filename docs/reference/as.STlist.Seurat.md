<div id="main" class="col-md-9" role="main">

# Convert Seurat Object to STlist

<div class="ref-description section level2">

**\[stable\]**

Converts a Seurat object (particularly spatial Seurat objects) to STlist
format.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
as.STlist.Seurat(
  x,
  assay = "RNA",
  slot = "counts",
  use.spatial = TRUE,
  verbose = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    a Seurat object to convert

-   assay:

    which assay to use for counts (default: 'RNA')

-   slot:

    which slot to use for counts (default: 'counts')

-   use.spatial:

    logical, whether to extract spatial coordinates

-   verbose:

    logical, whether to print progress messages

</div>

<div class="section level2">

## Value

an STlist object

</div>

<div class="section level2">

## Details

This function extracts data from a Seurat object and creates an STlist:

-   Counts from the default assay (usually 'RNA' or 'SCT')

-   Spatial coordinates from images slot (for Visium/Visium HD)

-   Sample metadata from Seurat meta.data

Works with both single-sample and multi-sample Seurat objects.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`as.Seurat.STlist`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
# Convert Seurat to STlist
st_obj <- as.STlist.Seurat(seurat_obj)

# Use in spatialGE workflow
st_obj <- SThet(st_obj, genes = c("GENE1", "GENE2"))
compare_SThet(st_obj)
} # }
```

</div>

</div>

</div>
