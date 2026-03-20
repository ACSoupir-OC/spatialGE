# STdiff_adjust_pvalues: Adjust p-values using multiple testing correction

Applies multiple testing correction to p-values using various methods

## Usage

``` r
STdiff_adjust_pvalues(pvals = NULL, method = "fdr")
```

## Arguments

- pvals:

  numeric vector of p-values to adjust

- method:

  adjustment method (default: 'fdr' for Benjamini-Hochberg)

## Value

numeric vector of adjusted p-values
