# STdiff_fit_spatial: Fit spatial mixed models for selected genes

Fits spatial mixed models with Matern covariance structure for genes
selected from non-spatial testing

## Usage

``` r
STdiff_fit_spatial(
  x = NULL,
  prep = NULL,
  sp_topgenes = 0.2,
  cores = NULL,
  verbose = 1L
)
```

## Arguments

- x:

  an STlist object

- prep:

  list from STdiff_select_genes containing combo_df, meta_dict,
  non_spatial_results

- sp_topgenes:

  proportion of top DE genes to fit spatial models on

- cores:

  number of cores for parallelization

- verbose:

  verbosity level

## Value

list containing spatial model results
