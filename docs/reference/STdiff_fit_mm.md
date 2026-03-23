<div id="main" class="col-md-9" role="main">

# STdiff_fit_mm: Fit non-spatial mixed models

<div class="ref-description section level2">

Fits non-spatial linear mixed models using spaMM for testing
differential expression between groups of spots/cells

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STdiff_fit_mm(expr_data = NULL, combo = NULL, pairwise = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   expr_data:

    a data frame (rows=spots/cells) with normalized expression data

-   combo:

    a data frame with columns (samplename, meta1, meta2, gene)
    containing all combinations of gene x cluster to test

-   pairwise:

    whether to perform pairwise tests (TRUE) or reference-based (FALSE)

</div>

<div class="section level2">

## Value

a list containing non-spatial model results with p-values, log-fold
changes, and metadata

</div>

</div>
