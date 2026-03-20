# STdiff_fit_spatial_models: Fit spaMM spatial models with Matern covariance

Fits spatial mixed models using spaMM::fitme with Matern covariance
structure to test for spatial autocorrelation in gene expression. The
model formula is: exprval ~ meta + Matern(1\|xpos + ypos) with fixed
nu=0.5.

This function handles both pairwise and reference-based comparisons,
applying appropriate data filtering and recoding for each model fit.

## Usage

``` r
STdiff_fit_spatial_models(
  expr_dat = NULL,
  non_sp_mods = NULL,
  annot_dict = NULL,
  verb = NULL,
  prep_pairwise = NULL
)
```

## Arguments

- expr_dat:

  a data frame with expression data and spatial coordinates (columns:
  meta, group, ypos, xpos, and gene expression columns)

- non_sp_mods:

  a list of non-spatial models containing metadata for each gene-cluster
  combination to test

- annot_dict:

  a data frame mapping coded annotations to original annotation names

- verb:

  verbosity level (0, 1, or 2)

- prep_pairwise:

  whether tests are pairwise (TRUE) or reference-based (FALSE)

## Value

list of spatial model results with the following structure:

- For successful fits: 'spmod' contains HLfit object

- For timeouts: 'spmod' = 'time_out'

- For convergence failures: 'spmod' = 'no_conv'

- Each element contains: spmod, sample, gene
