<div id="main" class="col-md-9" role="main">

# STclust_calculate_distances: Calculate scaled distance matrices

<div class="ref-description section level2">

Calculate and scale gene expression and spatial distance matrices for a
single sample

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STclust_calculate_distances(trcounts_df, coord_dat, dist_metric)
```

</div>

</div>

<div class="section level2">

## Arguments

-   trcounts_df:

    a filtered expression matrix for one sample

-   coord_dat:

    a data frame with spot coordinates for one sample

-   dist_metric:

    distance metric to use

</div>

<div class="section level2">

## Value

a list with scaled distance matrices (expression + spatial)

</div>

</div>
