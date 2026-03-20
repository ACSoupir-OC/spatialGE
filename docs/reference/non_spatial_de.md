# non_spatial_de: Run non-spatial linear models for gene x cluster combinations

Runs non-spatial linear models (or t-tests/Wilcoxon) for testing
differentially expressed genes between groups of spots/cells

## Usage

``` r
non_spatial_de(expr_data = NULL, combo = NULL, pairwise = NULL)
```

## Arguments

- expr_data:

  a data frame (rows=spots/cells) with normalized expression data as
  well as coordinates (x,y), and a group variable for random effects.

- combo:

  a data frame with columns (samplename, meta1, meta2, gene) containing
  all combinations of gene x cluster to test.

- pairwise:

  whether to perform pairwise tests (TRUE) or reference-based (FALSE)

## Value

a list containing non-spatial model results with p-values, log-fold
changes, and metadata
