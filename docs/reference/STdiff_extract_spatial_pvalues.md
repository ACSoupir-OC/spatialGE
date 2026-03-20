# STdiff_extract_spatial_pvalues: Extract spatial p-values from model results

Extracts spatial p-values and associated statistics from fitted spatial
models, handling convergence failures and generating adjusted p-values.

This is a convenience function that wraps STdiff_check_convergence and
provides a simplified output focused on p-values for downstream
analysis.

## Usage

``` r
STdiff_extract_spatial_pvalues(
  sp_models = NULL,
  meta_dict = NULL,
  combo_df = NULL,
  pairwise = NULL,
  pval_adj = "fdr"
)
```

## Arguments

- sp_models:

  list of spatial models from STdiff_run_spatial()

- meta_dict:

  annotation dictionary from STdiff_select_genes()

- combo_df:

  combinations dataframe from STdiff_select_genes()

- pairwise:

  whether tests are pairwise (TRUE) or reference-based (FALSE)

- pval_adj:

  p-value adjustment method (default: 'fdr')

## Value

data frame with columns:

- sample: sample name

- gene: gene name

- cluster_1: first cluster (original annotation name)

- cluster_2: second cluster (only if pairwise=TRUE)

- exp_p_val: raw spatial p-value

- exp_adj_p_val: adjusted spatial p-value

- exp_estimate: coefficient estimate for spatial effect

- convergence_status: model convergence status

- comments_spatial: warning/error message if any

## Examples

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
