# SThet_invdist_test: Computes spatial autocorrelation with statistical tests

Alternative implementation using spdep::moran.test and spdep::geary.test

## Usage

``` r
SThet_invdist_test(
  x = NULL,
  genes = NULL,
  samples = NULL,
  method = "moran",
  k = NULL,
  overwrite = T,
  cores = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  an STlist

- genes:

  a vector of gene names to compute statistics

- samples:

  the samples to compute statistics

- method:

  The spatial statistic(s) to estimate. It can be set to 'moran',
  'geary' or both. Default is 'moran'

- k:

  the number of neighbors to estimate weights

- overwrite:

  logical indicating if previous statistics should be overwritten

- cores:

  the number of cores to use during computations

## Value

an STlist containing spatial statistics
