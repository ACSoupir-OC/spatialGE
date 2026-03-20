# STclust_select_genes: Select variable genes for clustering

Selects top variable genes based on Seurat's VST method and filters
expression data

## Usage

``` r
STclust_select_genes(x, samples, topgenes, cores)
```

## Arguments

- x:

  an STlist object with normalized expression data

- samples:

  vector of sample names to process

- topgenes:

  number of top variable genes to select per sample

- cores:

  number of cores for parallelization

## Value

an STlist with updated `@gene_meta` (VST variance) and `@tr_counts`
filtered to top genes
