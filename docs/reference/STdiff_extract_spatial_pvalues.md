<div id="main" class="col-md-9" role="main">

# STdiff_extract_spatial_pvalues: Extract spatial p-values from model results

<div class="ref-description section level2">

Extracts spatial p-values and associated statistics from fitted spatial
models, handling convergence failures and generating adjusted p-values.

This is a convenience function that wraps STdiff_check_convergence and
provides a simplified output focused on p-values for downstream
analysis.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STdiff_extract_spatial_pvalues(
  sp_models = NULL,
  meta_dict = NULL,
  combo_df = NULL,
  pairwise = NULL,
  pval_adj = "fdr"
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   sp_models:

    list of spatial models from STdiff_run_spatial()

-   meta_dict:

    annotation dictionary from STdiff_select_genes()

-   combo_df:

    combinations dataframe from STdiff_select_genes()

-   pairwise:

    whether tests are pairwise (TRUE) or reference-based (FALSE)

-   pval_adj:

    p-value adjustment method (default: 'fdr')

</div>

<div class="section level2">

## Value

data frame with columns:

-   sample: sample name

-   gene: gene name

-   cluster_1: first cluster (original annotation name)

-   cluster_2: second cluster (only if pairwise=TRUE)

-   exp_p\_val: raw spatial p-value

-   exp_adj_p\_val: adjusted spatial p-value

-   exp_estimate: coefficient estimate for spatial effect

-   convergence_status: model convergence status

-   comments_spatial: warning/error message if any

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
# Extract spatial p-values for DE gene ranking
spatial_pvals = STdiff_extract_spatial_pvalues(
  sp_models=spatial_res$sp_models,
  meta_dict=prep$meta_dict,
  combo_df=prep$combo_df,
  pairwise=FALSE,
  pval_adj='fdr'
)
# View top spatial DE genes
head(spatial_pvals[order(spatial_pvals$exp_adj_p_val), ], 10)
} # }
```

</div>

</div>

</div>
