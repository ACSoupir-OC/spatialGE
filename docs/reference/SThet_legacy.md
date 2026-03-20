# SThet_legacy: Computes global spatial autocorrelation statistics on gene expression (legacy)

**\[superseded\]**

Computes the global spatial autocorrelation statistics Moran's I and/or
Geary's C for a set of genes

## Usage

``` r
SThet_legacy(
  x = NULL,
  genes = NULL,
  samples = NULL,
  method = "moran",
  k = NULL,
  overwrite = TRUE,
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

  the number of neighbors to estimate weights. By default NULL, meaning
  that spatial weights will be estimated from Euclidean distances. If an
  positive integer is entered, then the faster k nearest-neighbors
  approach is used. Please keep in mind that estimates are not as
  accurate as when using the default distance-based method.

- overwrite:

  logical indicating if previous statistics should be overwritten.
  Default to FALSE (do not overwrite)

- cores:

  integer indicating the number of cores to use during parallelization.
  If NULL, the function uses half of the available cores at a maximum.
  The parallelization uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) and
  works only in Unix systems

- verbose:

  logical, whether to print text to console

## Value

an STlist containing spatial statistics

## Details

**This is the legacy implementation from version 1.x.** Use `SThet` for
new projects.

The function computes global spatial autocorrelation statistics (Moran's
I and/or Geary's C) for the requested genes and samples. Then
computation uses the package `spdep`. The calculated statistics are
stored in the STlist, which can be accessed with the `get_gene_meta`
function. For visual comparative analysis, the function `compare_SThet`
can be used afterwards.
