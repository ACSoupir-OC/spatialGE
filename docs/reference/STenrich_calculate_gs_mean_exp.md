# STenrich_calculate_gs_mean_exp: Calculate average gene set expression

Calculate mean expression across genes in each gene set

## Usage

``` r
STenrich_calculate_gs_mean_exp(x, combo, pw_genes, min_genes, cores)
```

## Arguments

- x:

  an STlist with transformed gene expression (as DelayedArray)

- combo:

  data frame with combinations of samples and gene sets

- pw_genes:

  list of available genes per gene set

- min_genes:

  the minimum number of genes of a gene set present in data (default: 5)

- cores:

  the number of cores for parallelization

## Value

list of data frames with average gene set expression
