# STgradient_calculate_correlations: Calculate correlations

For each gene, fit linear models and calculate Spearman correlations

## Usage

``` r
STgradient_calculate_correlations(
  vargenes_expr,
  outs_dist2ref,
  robust,
  log_dist,
  sample_name
)
```

## Arguments

- vargenes_expr:

  data.frame with expression and distance data

- outs_dist2ref:

  list of outlier barcodes per gene

- robust:

  logical whether to use robust regression

- log_dist:

  logical whether to log-transform distances

- sample_name:

  name of sample

## Value

data.frame with correlation results
