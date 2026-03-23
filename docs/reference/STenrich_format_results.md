<div id="main" class="col-md-9" role="main">

# STenrich_format_results: Format final results

<div class="ref-description section level2">

Compile results and adjust p-values

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STenrich_format_results(pval_res, samples, pval_adj_method)
```

</div>

</div>

<div class="section level2">

## Arguments

-   pval_res:

    data frame with p-values per sample and gene set

-   samples:

    a vector with sample names to run analysis

-   pval_adj_method:

    the method for multiple comparison adjustment (default: 'BH')

</div>

<div class="section level2">

## Value

list of data frames with p-values and adjusted p-values per sample

</div>

</div>
