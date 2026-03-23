<div id="main" class="col-md-9" role="main">

# stdiff_mean_test: Perform t-test or Wilcoxon test for differential expression

<div class="ref-description section level2">

Performs simple mean-based tests (t-test or Wilcoxon) between groups of
spots/cells for a given gene

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
stdiff_mean_test(
  expr_data = NULL,
  combo = NULL,
  pairwise = NULL,
  test_type = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   expr_data:

    a data frame with expression data and metadata

-   combo:

    a data frame with columns (samplename, meta1, meta2, gene)

-   pairwise:

    whether to perform pairwise tests (TRUE) or reference-based (FALSE)

-   test_type:

    type of test: 't_test' or 'wilcoxon'

</div>

<div class="section level2">

## Value

a list containing test results with p-values, log-fold changes, and
metadata

</div>

</div>
