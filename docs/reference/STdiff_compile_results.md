# STdiff_compile_results: Compile non-spatial and spatial results

Main entry point for compiling differential expression results. Merges
non-spatial model results (linear models, t-tests, or Wilcoxon tests)
with spatial mixed model results containing Matern covariance structure.
Handles convergence issues, time-outs, and applies appropriate p-value
adjustments. Returns a list of data frames with complete DE statistics
per sample.

## Usage

``` r
STdiff_compile_results(prep = NULL, spatial = NULL)
```

## Arguments

- prep:

  list from STdiff_select_genes containing:

  - non_spatial_results: list of data frames with non-spatial model
    results

  - meta_dict: data frame mapping coded annotations to original names

  - pairwise: logical indicating test type

  - pval_adj: p-value adjustment method

- spatial:

  list from STdiff_fit_spatial containing:

  - sp_models: list of spatial model results (or NULL if not run)

  - meta_dict: same meta_dict from prep

  - combo_df: combinations data frame

## Value

list of data frames with final differential expression results per
sample. Each data frame contains columns: sample, gene, cluster_1,
cluster_2 (if pairwise), avg_log2fc, mm_p_val/ttest_p_val/wilcox_p_val,
adj_p_val, exp_p_val, exp_adj_p_val, comments (indicating spatial model
status: 'no_spatial_test', 'time_out', 'no_convergence', or NA for
successful spatial models)

## Examples

``` r
# Example usage (internal function, called by STdiff())
# result_de = STdiff_compile_results(prep=prep, spatial=spatial_res)
```
