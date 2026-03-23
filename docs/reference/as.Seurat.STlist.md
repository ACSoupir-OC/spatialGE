<div id="main" class="col-md-9" role="main">

# Convert STlist to Seurat Object

<div class="ref-description section level2">

**\[stable\]**

Converts an STlist object to a Seurat object for integration with Seurat
workflows.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
as.Seurat.STlist(
  x,
  samples = NULL,
  assay.name = "RNA",
  add.spatial.info = TRUE,
  verbose = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    an STlist object to convert

-   samples:

    vector of sample names to include. If NULL, all samples are
    converted

-   assay.name:

    name for the assay in the Seurat object (default: 'RNA')

-   add.spatial.info:

    logical, whether to add spatial coordinates to Seurat object

-   verbose:

    logical, whether to print progress messages

</div>

<div class="section level2">

## Value

a Seurat object containing the STlist data

</div>

<div class="section level2">

## Details

This function creates a Seurat object from an STlist, preserving:

-   Count data in the 'RNA' assay

-   Spatial coordinates in the 'spatial' reduction

-   Sample metadata in Seurat meta.data

-   Gene metadata in feature-level metadata

For multi-sample STlist objects, returns a merged Seurat object with
sample barcodes preserved.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`as.STlist.Seurat`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
# Convert STlist to Seurat
seurat_obj <- as.Seurat.STlist(st_obj)

# Convert specific samples
seurat_obj <- as.Seurat.STlist(st_obj, samples = c("sample1", "sample2"))

# Use in Seurat workflow
seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
} # }
```

</div>

</div>

</div>
