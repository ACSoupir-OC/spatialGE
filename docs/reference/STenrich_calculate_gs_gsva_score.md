# STenrich_calculate_gs_gsva_score: Calculate GSVA scores

Calculate GSVA enrichment scores for gene sets

## Usage

``` r
STenrich_calculate_gs_gsva_score(
  x,
  pw_genes,
  gene_sets,
  min_genes,
  cores,
  verbose = TRUE
)
```

## Arguments

- x:

  an STlist with transformed gene expression (as DelayedArray)

- pw_genes:

  list of available genes per gene set

- gene_sets:

  a named list of gene sets to test

- min_genes:

  the minimum number of genes of a gene set present in data (default: 5)

- cores:

  the number of cores for parallelization

- verbose:

  verbosity level

## Value

list of data frames with GSVA scores
