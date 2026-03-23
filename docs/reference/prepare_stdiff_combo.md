<div id="main" class="col-md-9" role="main">

# prepare_stdiff_combo: Create combinations for differential testing

<div class="ref-description section level2">

Creates a data frame with all combinations of genes x cluster pairs to
test for differential expression

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
prepare_stdiff_combo(
  to_expand = NULL,
  user_clusters = NULL,
  pairwise = NULL,
  verbose = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   to_expand:

    a data frame with columns (samplename, orig_annot, gene, meta)
    containing all gene x annotation combinations per sample

-   user_clusters:

    optional vector of specific clusters to test

-   pairwise:

    whether to perform pairwise tests (TRUE) or reference-based (FALSE)

-   verbose:

    verbosity level for progress messages

</div>

<div class="section level2">

## Value

a data frame with all combinations to test (samplename, meta1, meta2,
gene)

</div>

</div>
