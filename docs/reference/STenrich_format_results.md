# STenrich_format_results: Format final results

Compile results and adjust p-values

## Usage

``` r
STenrich_format_results(pval_res, samples, pval_adj_method)
```

## Arguments

- pval_res:

  data frame with p-values per sample and gene set

- samples:

  a vector with sample names to run analysis

- pval_adj_method:

  the method for multiple comparison adjustment (default: 'BH')

## Value

list of data frames with p-values and adjusted p-values per sample
