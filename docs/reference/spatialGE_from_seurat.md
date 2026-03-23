<div id="main" class="col-md-9" role="main">

# Run spatialGE analysis on Seurat object

<div class="ref-description section level2">

**\[stable\]**

Convenience wrapper that converts Seurat to STlist and runs spatial
analysis.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
spatialGE_from_seurat(
  x,
  analysis = "SThet",
  genes = NULL,
  assay = "RNA",
  verbose = TRUE,
  ...
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    a Seurat object (spatial or non-spatial)

-   analysis:

    type of analysis: 'SThet', 'STclust', 'STdiff', 'STenrich',
    'STgradient'

-   genes:

    vector of gene names to analyze (required for SThet, STgradient)

-   assay:

    which assay to use for counts (default: 'RNA')

-   verbose:

    logical, whether to print progress

-   ...:

    additional arguments passed to the analysis function

</div>

<div class="section level2">

## Value

depends on analysis type:

-   SThet: STlist with gene_meta updated

-   STclust: STlist with spatial_meta updated (cluster assignments)

-   STdiff: list of differential expression results

-   STenrich: list of enrichment results

-   STgradient: list of gradient analysis results

</div>

<div class="section level2">

## Details

This function provides a seamless workflow from Seurat to spatialGE:

1.  Converts Seurat object to STlist

2.  Runs specified spatial analysis (SThet, STclust, STdiff, etc.)

3.  Optionally converts results back to Seurat format

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
# Run heterogeneity analysis
result <- spatialGE_from_seurat(seurat_obj, analysis = 'SThet', 
                                 genes = c('GENE1', 'GENE2'))

# Run clustering
result <- spatialGE_from_seurat(seurat_obj, analysis = 'STclust')

# Run differential expression
result <- spatialGE_from_seurat(seurat_obj, analysis = 'STdiff',
                                 annot = 'cell_type')
} # }
```

</div>

</div>

</div>
