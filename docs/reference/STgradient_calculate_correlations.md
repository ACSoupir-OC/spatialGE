<div id="main" class="col-md-9" role="main">

# STgradient_calculate_correlations: Calculate correlations

<div class="ref-description section level2">

For each gene, fit linear models and calculate Spearman correlations

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STgradient_calculate_correlations(
  vargenes_expr,
  outs_dist2ref,
  robust,
  log_dist,
  sample_name
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   vargenes_expr:

    data.frame with expression and distance data

-   outs_dist2ref:

    list of outlier barcodes per gene

-   robust:

    logical whether to use robust regression

-   log_dist:

    logical whether to log-transform distances

-   sample_name:

    name of sample

</div>

<div class="section level2">

## Value

data.frame with correlation results

</div>

</div>
