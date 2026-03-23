<div id="main" class="col-md-9" role="main">

# compare_SThet: Compares spatial autocorrelation statistics across samples

<div class="ref-description section level2">

Plots the spatial autocorrelation statistics of genes across samples and
colors samples acording to sample metadata.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
compare_SThet(
  x = NULL,
  samplemeta = NULL,
  genes = NULL,
  color_by = NULL,
  categorical = TRUE,
  color_pal = "muted",
  ptsize = 1
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    an STlist.

-   samplemeta:

    a string indicating the name of the variable in the clinical data
    frame. If NULL, uses sample names

-   genes:

    the name(s) of the gene(s) to plot.

-   color_by:

    the variable in `x@spatial_meta` used to color points in the plot.
    If NULL, each sample is assigned a different color

-   categorical:

    logical indicating whether or not to treat `color_by` as a
    categorical variable. Default is TRUE

-   color_pal:

    a string of a color palette from khroma or RColorBrewer, or a vector
    with colors with enough elements to plot categories.

-   ptsize:

    a number specifying the size of the points. Passed to the `size`
    aesthetic.

</div>

<div class="section level2">

## Value

a list of plots

</div>

<div class="section level2">

## Details

This function takes the names of genes and their Moran's I or Geary's C
computed for multiple samples and to provide a muti-sample comparison.
Samples in the plot can be colored according to sample metadata to
explore potential associations between spatial distribution of gene
expression and sample-level data.

</div>

</div>
