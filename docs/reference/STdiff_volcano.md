<div id="main" class="col-md-9" role="main">

# STdiff_volcano: Generates volcano plots from STdiff results

<div class="ref-description section level2">

Generates volcano plots of differential expression results from STdiff

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STdiff_volcano(
  x = NULL,
  samples = NULL,
  clusters = NULL,
  pval_thr = 0.05,
  color_pal = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    the output of `STdiff`

-   samples:

    samples to create plots

-   clusters:

    names of the clusters to generate comparisons

-   pval_thr:

    the p-value threshold to color genes with differential expression

-   color_pal:

    the palette to color genes by significance

</div>

<div class="section level2">

## Value

a list of ggplot objects

</div>

<div class="section level2">

## Details

The function generated volcano plots (p-value vs. log-fold change) for
genes tested with `STdiff`. Colors can be customized to show
significance from spatial and non-spatial models

</div>

</div>
