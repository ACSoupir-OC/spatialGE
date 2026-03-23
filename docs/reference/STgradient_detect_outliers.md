<div id="main" class="col-md-9" role="main">

# STgradient_detect_outliers: Detect expression outliers

<div class="ref-description section level2">

Identify outlier spots using IQR-based method for each gene

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STgradient_detect_outliers(vargenes_expr, out_rm)
```

</div>

</div>

<div class="section level2">

## Arguments

-   vargenes_expr:

    data.frame with expression and distance data

-   out_rm:

    logical whether to detect outliers

</div>

<div class="section level2">

## Value

list with outlier barcodes per gene

</div>

</div>
