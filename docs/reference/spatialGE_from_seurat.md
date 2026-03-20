# Run spatialGE analysis on Seurat object

**\[stable\]**

Convenience wrapper that converts Seurat to STlist and runs spatial
analysis.

## Usage

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

## Arguments

- x:

  a Seurat object (spatial or non-spatial)

- analysis:

  type of analysis: 'SThet', 'STclust', 'STdiff', 'STenrich',
  'STgradient'

- genes:

  vector of gene names to analyze (required for SThet, STgradient)

- assay:

  which assay to use for counts (default: 'RNA')

- verbose:

  logical, whether to print progress

- ...:

  additional arguments passed to the analysis function

## Value

depends on analysis type:

- SThet: STlist with gene_meta updated

- STclust: STlist with spatial_meta updated (cluster assignments)

- STdiff: list of differential expression results

- STenrich: list of enrichment results

- STgradient: list of gradient analysis results

## Details

This function provides a seamless workflow from Seurat to spatialGE:

1.  Converts Seurat object to STlist

2.  Runs specified spatial analysis (SThet, STclust, STdiff, etc.)

3.  Optionally converts results back to Seurat format

## See also

[`as.STlist.Seurat`](https://acsoupir-oc.github.io/spatialGE/reference/as.STlist.Seurat.md),
[`as.Seurat.STlist`](https://acsoupir-oc.github.io/spatialGE/reference/as.Seurat.STlist.md)

## Examples

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
