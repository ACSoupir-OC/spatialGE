# stdiff_mean_test: Perform t-test or Wilcoxon test for differential expression

Performs simple mean-based tests (t-test or Wilcoxon) between groups of
spots/cells for a given gene

## Usage

``` r
stdiff_mean_test(
  expr_data = NULL,
  combo = NULL,
  pairwise = NULL,
  test_type = NULL
)
```

## Arguments

- expr_data:

  a data frame with expression data and metadata

- combo:

  a data frame with columns (samplename, meta1, meta2, gene)

- pairwise:

  whether to perform pairwise tests (TRUE) or reference-based (FALSE)

- test_type:

  type of test: 't_test' or 'wilcoxon'

## Value

a list containing test results with p-values, log-fold changes, and
metadata
