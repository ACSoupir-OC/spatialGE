<div id="main" class="col-md-9" role="main">

# STdiff_run_nonspatial: Run non-spatial differential expression tests

<div class="ref-description section level2">

Main entry point for non-spatial DE testing. Runs linear models,
t-tests, or Wilcoxon tests for selected genes between groups of
spots/cells

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STdiff_run_nonspatial(
  x = NULL,
  samples = NULL,
  annot = NULL,
  w = NULL,
  k = NULL,
  deepSplit = NULL,
  topgenes = 5000,
  pval_thr = 0.05,
  pval_adj = "fdr",
  test_type = "mm",
  clusters = NULL,
  pairwise = FALSE,
  verbose = 1L,
  cores = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    an STlist object

-   samples:

    vector of sample names to test

-   annot:

    column name in spatial_meta for cluster annotations

-   w:

    spatial weight parameter (used if annot is NULL)

-   k:

    number of clusters (used if annot is NULL)

-   deepSplit:

    deepSplit parameter for dynamicTreeCut clusters

-   topgenes:

    number of top variable genes to select

-   pval_thr:

    p-value threshold for selecting DE genes

-   pval_adj:

    p-value adjustment method

-   test_type:

    type of test: 'mm', 't_test', or 'wilcoxon'

-   clusters:

    optional vector of specific clusters to test

-   pairwise:

    whether to perform pairwise tests (TRUE) or reference-based (FALSE)

-   verbose:

    verbosity level (0, 1, or 2)

-   cores:

    number of cores for parallelization

</div>

<div class="section level2">

## Value

list containing: combo_df, meta_dict, non_spatial_results, pval_thr,
test_type, pairwise

</div>

</div>
