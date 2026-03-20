# STdiff_run_spatial: Run spatial differential expression analysis

Main entry point for spatial mixed model fitting. This function
orchestrates the fitting of spatial models with Matern covariance
structure for genes selected from non-spatial differential expression
testing.

The spatial models test for spatial autocorrelation in gene expression
between clusters using the formula: exprval ~ meta + Matern(1\|xpos +
ypos) with fixed nu parameter (nu=0.5) for Matern covariance.

## Usage

``` r
STdiff_run_spatial(
  x = NULL,
  prep = NULL,
  sp_topgenes = 0.2,
  cores = NULL,
  verbose = 1L
)
```

## Arguments

- x:

  an STlist object containing the spatial transcriptomics data

- prep:

  list from STdiff_select_genes containing combo_df, meta_dict,
  non_spatial_results, pval_thr, test_type, and pairwise

- sp_topgenes:

  proportion (0-1) of top DE genes to fit spatial models on (default 0.2
  = 20%)

- cores:

  number of CPU cores for parallelization (NULL = auto-detect)

- verbose:

  verbosity level (0 = silent, 1 = progress, 2 = detailed)

## Value

list containing:

- sp_models: nested list of fitted spatial models (HLfit objects) or
  status strings

- meta_dict: annotation dictionary mapping coded to original cluster
  names

- combo_df: combinations dataframe for reference

## Examples

``` r
if (FALSE) { # \dontrun{
# Run spatial DE analysis on top 20% of DE genes
spatial_res = STdiff_run_spatial(x=stlist_obj, prep=prep_res,
                                 sp_topgenes=0.2, cores=4, verbose=1)
} # }
```
