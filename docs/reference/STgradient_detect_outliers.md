# STgradient_detect_outliers: Detect expression outliers

Identify outlier spots using IQR-based method for each gene

## Usage

``` r
STgradient_detect_outliers(vargenes_expr, out_rm)
```

## Arguments

- vargenes_expr:

  data.frame with expression and distance data

- out_rm:

  logical whether to detect outliers

## Value

list with outlier barcodes per gene
