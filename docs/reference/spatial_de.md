<div id="main" class="col-md-9" role="main">

# spatial_de: Run spatial linear models with Matern covariance

<div class="ref-description section level2">

Fits spatial mixed models with Matern covariance structure to test for
spatial autocorrelation in gene expression between clusters

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
spatial_de(expr_dat = NULL, non_sp_mods = NULL, annot_dict = NULL, verb = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   expr_dat:

    a data frame with expression data and coordinates

-   non_sp_mods:

    a list of non-spatial models to update with covariance structures

-   annot_dict:

    a data frame mapping coded annotations to original annotation names

-   verb:

    verbosity level (0, 1, or 2)

</div>

<div class="section level2">

## Value

a list of spatial models with convergence status

</div>

</div>
