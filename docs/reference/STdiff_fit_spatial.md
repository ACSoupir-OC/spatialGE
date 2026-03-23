<div id="main" class="col-md-9" role="main">

# STdiff_fit_spatial: Fit spatial mixed models for selected genes

<div class="ref-description section level2">

Fits spatial mixed models with Matern covariance structure for genes
selected from non-spatial testing

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STdiff_fit_spatial(
  x = NULL,
  prep = NULL,
  sp_topgenes = 0.2,
  cores = NULL,
  verbose = 1L
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    an STlist object

-   prep:

    list from STdiff_select_genes containing combo_df, meta_dict,
    non_spatial_results

-   sp_topgenes:

    proportion of top DE genes to fit spatial models on

-   cores:

    number of cores for parallelization

-   verbose:

    verbosity level

</div>

<div class="section level2">

## Value

list containing spatial model results

</div>

</div>
