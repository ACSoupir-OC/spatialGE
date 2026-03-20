# STdiff_select_genes: Gene filtering and non-spatial differential testing

Selects variable genes and performs non-spatial differential expression
testing (linear models, t-tests, or Wilcoxon tests)

## Usage

``` r
STdiff_select_genes(
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

## Arguments

- x:

  an STlist object

- samples:

  vector of sample names to test

- annot:

  column name in spatial_meta for cluster annotations

- w:

  spatial weight parameter (used if annot is NULL)

- k:

  number of clusters (used if annot is NULL)

- deepSplit:

  deepSplit parameter for dynamicTreeCut clusters

- topgenes:

  number of top variable genes to select

- pval_thr:

  p-value threshold for selecting DE genes

- pval_adj:

  p-value adjustment method

- test_type:

  type of test: 'mm', 't_test', or 'wilcoxon'

- clusters:

  optional vector of specific clusters to test

- pairwise:

  whether to perform pairwise tests

- verbose:

  verbosity level

- cores:

  number of cores for parallelization

## Value

list containing: combo_df, meta_dict, non_spatial_results, pval_thr,
test_type
