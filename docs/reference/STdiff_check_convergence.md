<div id="main" class="col-md-9" role="main">

# STdiff_check_convergence: Check model convergence and extract summary statistics

<div class="ref-description section level2">

Checks the convergence status of fitted spatial models and extracts
summary statistics including p-values, coefficient estimates, and
warnings.

This function handles multiple convergence states:

-   Successful fits (HLfit objects)

-   Timeouts ('time_out')

-   Non-convergence ('no_convergence')

-   Unknown errors ('unknown_error')

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STdiff_check_convergence(
  sp_models = NULL,
  meta_dict = NULL,
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

    annotation dictionary mapping coded to original cluster names

-   pairwise:

    whether tests are pairwise (TRUE) or reference-based (FALSE)

-   pval_adj:

    p-value adjustment method for spatial p-values

</div>

<div class="section level2">

## Value

list containing:

-   sp_de: data frame with spatial model results (p-values,
    coefficients, convergence status)

-   convergence_summary: summary of convergence status across all models

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
# After running STdiff_run_spatial, check convergence
conv_result = STdiff_check_convergence(sp_models=spatial_res$sp_models,
                                       meta_dict=prep$meta_dict,
                                       pairwise=FALSE)
print(conv_result$convergence_summary)
} # }
```

</div>

</div>

</div>
