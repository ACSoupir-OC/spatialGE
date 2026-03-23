<div id="main" class="col-md-9" role="main">

# STdiff_fit_ttest: Fit non-spatial t-tests

<div class="ref-description section level2">

Performs Welch's t-tests for testing differential expression between
groups of spots/cells

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STdiff_fit_ttest(expr_data = NULL, combo = NULL, pairwise = NULL)
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

</div>

<div class="section level2">

## Value

a list containing test results with p-values, log-fold changes, and
metadata

</div>

</div>
