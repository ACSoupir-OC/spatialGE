# STdiff_fit_wilcoxon: Fit non-spatial Wilcoxon tests

Performs Wilcoxon rank-sum tests for testing differential expression
between groups of spots/cells (non-parametric alternative to t-test)

## Usage

``` r
STdiff_fit_wilcoxon(expr_data = NULL, combo = NULL, pairwise = NULL)
```

## Arguments

- expr_data:

  a data frame with expression data and metadata

- combo:

  a data frame with columns (samplename, meta1, meta2, gene)

- pairwise:

  whether to perform pairwise tests (TRUE) or reference-based (FALSE)

## Value

a list containing test results with p-values, log-fold changes, and
metadata
