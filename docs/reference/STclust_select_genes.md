<div id="main" class="col-md-9" role="main">

# STclust_select_genes: Select variable genes for clustering

<div class="ref-description section level2">

Selects top variable genes based on Seurat's VST method and filters
expression data

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STclust_select_genes(x, samples, topgenes, cores)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    an STlist object with normalized expression data

-   samples:

    vector of sample names to process

-   topgenes:

    number of top variable genes to select per sample

-   cores:

    number of cores for parallelization

</div>

<div class="section level2">

## Value

an STlist with updated `@gene_meta` (VST variance) and `@tr_counts`
filtered to top genes

</div>

</div>
