# Convert Seurat Object to STlist

**\[stable\]**

Converts a Seurat object (particularly spatial Seurat objects) to STlist
format.

## Usage

``` r
as.STlist.Seurat(
  x,
  assay = "RNA",
  slot = "counts",
  use.spatial = TRUE,
  verbose = TRUE
)
```

## Arguments

- x:

  a Seurat object to convert

- assay:

  which assay to use for counts (default: 'RNA')

- slot:

  which slot to use for counts (default: 'counts')

- use.spatial:

  logical, whether to extract spatial coordinates

- verbose:

  logical, whether to print progress messages

## Value

an STlist object

## Details

This function extracts data from a Seurat object and creates an STlist:

- Counts from the default assay (usually 'RNA' or 'SCT')

- Spatial coordinates from images slot (for Visium/Visium HD)

- Sample metadata from Seurat meta.data

Works with both single-sample and multi-sample Seurat objects.

## See also

[`as.Seurat.STlist`](https://acsoupir-oc.github.io/spatialGE/reference/as.Seurat.STlist.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert Seurat to STlist
st_obj <- as.STlist.Seurat(seurat_obj)

# Use in spatialGE workflow
st_obj <- SThet(st_obj, genes = c("GENE1", "GENE2"))
compare_SThet(st_obj)
} # }
```
