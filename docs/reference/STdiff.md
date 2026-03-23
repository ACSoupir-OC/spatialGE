<div id="main" class="col-md-9" role="main">

# STdiff: Differential gene expression analysis for spatial transcriptomics

<div class="ref-description section level2">

**\[stable\]**

Tests for differentially expressed genes between groups of spots/cells
in spatial transcriptomics data. First runs non-spatial tests (linear
models, t-tests, or Wilcoxon tests) to detect DE genes, then optionally
fits spatial mixed models with Matern covariance structure on a subset
of DE genes.

Public wrapper function that calls the modular STdiff functions

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STdiff(
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
  sp_topgenes = 0.2,
  clusters = NULL,
  pairwise = FALSE,
  verbose = 1L,
  cores = NULL
)

STdiff(
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
  sp_topgenes = 0.2,
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

    an STlist object containing spatial transcriptomics data

-   samples:

    vector of sample names or indices to test. If NULL, uses all samples

-   annot:

    column name in spatial_meta for cluster annotations. Required if w/k
    not specified

-   w:

    spatial weight parameter for STclust (required if annot is NULL)

-   k:

    number of clusters for STclust, or 'dtc' for dynamicTreeCut
    (required if annot is NULL)

-   deepSplit:

    deepSplit parameter for dynamicTreeCut clusters (required if
    k='dtc')

-   topgenes:

    number of top variable genes to select (default=5000). If NULL, all
    genes used

-   pval_thr:

    p-value threshold for selecting DE genes from non-spatial tests
    (default=0.05)

-   pval_adj:

    p-value adjustment method (default='fdr', passed to stats::p.adjust)

-   test_type:

    type of test: 'mm' (linear models), 't_test', or 'wilcoxon'

-   sp_topgenes:

    proportion (0-1) of DE genes to fit spatial models on. If 0, no
    spatial tests

-   clusters:

    optional vector of specific cluster names to test (vs all clusters)

-   pairwise:

    whether to perform pairwise tests (TRUE) or reference-based (FALSE)

-   verbose:

    verbosity level (0=silent, 1=progress, 2=detailed)

-   cores:

    number of cores for parallelization. If NULL, auto-detected

</div>

<div class="section level2">

## Value

a list with one data frame per sample containing differential expression
results. Each data frame includes columns: sample, gene, cluster_1,
cluster_2 (if pairwise), avg_log2fc,
mm_p\_val/ttest_p\_val/wilcox_p\_val, adj_p\_val, exp_p\_val,
exp_adj_p\_val, comments (indicating spatial model status:
'no_spatial_test', 'time_out', 'no_convergence', or NA for successful
spatial models)

</div>

<div class="section level2">

## Details

Tests for genes with significantly higher or lower expression in one
group of spots/cells with respect to the rest. Users can specify cluster
assignments via metadata columns or via STclust parameters (w, k,
deepSplit). Uses spaMM::fitme and is computationally expensive. Set
sp_topgenes=0 to run non-spatial tests only.

Maintains the original API while using the new modular implementation

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
# Load TNBC test data
data_dir <- "tests/testthat/data/tnbc_bassiouni"
count_files <- list.files(data_dir, pattern='counts', full.names=TRUE)[1:2]
coord_files <- list.files(data_dir, pattern='mapping', full.names=TRUE)[1:2]
clin_file <- file.path(data_dir, "bassiouni_clinical.csv")

# Create STlist
tnbc <- STlist(rnacounts=count_files, spotcoords=coord_files, samples=clin_file)
tnbc <- transform_data(tnbc)

# Run non-spatial DE (Wilcoxon)
tnbc <- STdiff(tnbc, test='wilcoxon', group='tissue_type')

# Extract and view results
res <- STdiff_compile_results(tnbc)
head(res)

# Alternative: Run non-spatial tests only
result = STdiff(x=stlist_obj, annot='cluster', test_type='mm', sp_topgenes=0)
} # }
```

</div>

</div>

</div>
