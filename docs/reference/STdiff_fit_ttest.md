# STdiff_fit_ttest: Fit non-spatial t-tests

Performs Welch's t-tests for testing differential expression between
groups of spots/cells

## Usage

``` r
STdiff_fit_ttest(expr_data = NULL, combo = NULL, pairwise = NULL)
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
