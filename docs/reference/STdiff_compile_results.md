<div id="main" class="col-md-9" role="main">

# STdiff_compile_results: Compile non-spatial and spatial results

<div class="ref-description section level2">

Main entry point for compiling differential expression results. Merges
non-spatial model results (linear models, t-tests, or Wilcoxon tests)
with spatial mixed model results containing Matern covariance structure.
Handles convergence issues, time-outs, and applies appropriate p-value
adjustments. Returns a list of data frames with complete DE statistics
per sample.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STdiff_compile_results(prep = NULL, spatial = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   prep:

    list from STdiff_select_genes containing:

    -   non_spatial_results: list of data frames with non-spatial model
        results

    -   meta_dict: data frame mapping coded annotations to original
        names

    -   pairwise: logical indicating test type

    -   pval_adj: p-value adjustment method

-   spatial:

    list from STdiff_fit_spatial containing:

    -   sp_models: list of spatial model results (or NULL if not run)

    -   meta_dict: same meta_dict from prep

    -   combo_df: combinations data frame

</div>

<div class="section level2">

## Value

list of data frames with final differential expression results per
sample. Each data frame contains columns: sample, gene, cluster_1,
cluster_2 (if pairwise), avg_log2fc,
mm_p\_val/ttest_p\_val/wilcox_p\_val, adj_p\_val, exp_p\_val,
exp_adj_p\_val, comments (indicating spatial model status:
'no_spatial_test', 'time_out', 'no_convergence', or NA for successful
spatial models)

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# Example usage (internal function, called by STdiff())
# result_de = STdiff_compile_results(prep=prep, spatial=spatial_res)
```

</div>

</div>

</div>
