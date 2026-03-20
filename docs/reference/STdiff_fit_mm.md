# STdiff_fit_mm: Fit non-spatial mixed models

Fits non-spatial linear mixed models using spaMM for testing
differential expression between groups of spots/cells

## Usage

``` r
STdiff_fit_mm(expr_data = NULL, combo = NULL, pairwise = NULL)
```

## Arguments

- expr_data:

  a data frame (rows=spots/cells) with normalized expression data

- combo:

  a data frame with columns (samplename, meta1, meta2, gene) containing
  all combinations of gene x cluster to test

- pairwise:

  whether to perform pairwise tests (TRUE) or reference-based (FALSE)

## Value

a list containing non-spatial model results with p-values, log-fold
changes, and metadata
