<div id="main" class="col-md-9" role="main">

# STgradient_core: Core algorithm for spatial gradient analysis

<div class="ref-description section level2">

Core implementation of STgradient - calculates Spearman correlations
between gene expression and distances to reference domain

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STgradient_core(
  x,
  samples,
  annot,
  ref,
  exclude,
  out_rm,
  limit,
  distsumm,
  min_nb,
  robust,
  nb_dist_thr,
  log_dist,
  topgenes,
  cores,
  verbose
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    STlist object with transformed gene expression

-   samples:

    Vector of sample names to process

-   annot:

    Name of annotation column in @spatial_meta

-   ref:

    Reference tissue domain

-   exclude:

    Optional domain to exclude

-   out_rm:

    Remove gene expression outliers (IQR method)

-   limit:

    Limit analysis to spots within this distance threshold

-   distsumm:

    Distance summary metric: "min" or "avg"

-   min_nb:

    Minimum number of neighbors required

-   robust:

    Use robust regression

-   nb_dist_thr:

    Neighborhood distance threshold

-   log_dist:

    Apply log transform to distances

-   topgenes:

    Number of high-variance genes to test

-   cores:

    Number of cores for parallelization

-   verbose:

    Print progress messages

</div>

<div class="section level2">

## Value

List of data frames with correlation results

</div>

<div class="section level2">

## Details

Internal function - use STgradient() for public interface

</div>

</div>
