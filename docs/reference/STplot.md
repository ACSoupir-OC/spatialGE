<div id="main" class="col-md-9" role="main">

# STplot: Plots of gene expression, cluster memberships, and metadata in spatial context

<div class="ref-description section level2">

**\[stable\]**

Generates a plot of the location of spots/cells within an spatial
sample, and colors them according to gene expression levels or
spot/cell-level metadata

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STplot(
  x,
  samples = NULL,
  genes = NULL,
  plot_meta = NULL,
  ks = "dtc",
  ws = NULL,
  deepSplit = NULL,
  color_pal = NULL,
  data_type = "tr",
  ptsize = NULL,
  txsize = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    an STlist

-   samples:

    a vector of numbers indicating the ST samples to plot, or their
    sample names. If vector of numbers, it follow the order of samples
    in `names(x@counts)`. If NULL, the function plots all samples

-   genes:

    a vector of gene names or a named list of gene sets. In the latter
    case, the averaged expression of genes within the sets is plotted

-   plot_meta:

    a column name in `x@spatial_meta` to plot

-   ks:

    the k values to plot or 'dtc' to plot results from `dynamicTreeCut`
    clustering solutions. Requires previous analysis with `STclust`

-   ws:

    the spatial weights to plot samples if `STclust` was used

-   deepSplit:

    a logical or positive number indicating the `deepSplit`, if samples
    were analyzed with `STclust`

-   color_pal:

    a string of a color palette from `khroma` or `RColorBrewer`, or a
    vector with enough color names or HEX values

-   data_type:

    one of 'tr' or 'raw', to plot transformed or raw counts respectively

-   ptsize:

    a number specifying the size of the points. Passed to the `size`

-   txsize:

    a number controlling the size of the text in the plot title and
    legend title. Passed to the `element_text` aesthetic.

</div>

<div class="section level2">

## Value

a list of plots

</div>

<div class="section level2">

## Details

The function takes an STlist and plots the cells or spots in their
spatial context. The users can color the spots/cells according to the
expression of selected genes, cluster memberships, or any spot/cell
level metadata included in `x@spatial_meta`. The function also can
average expression of gene sets.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# Load TNBC dataset (from spatialGE_Data package)
# library(spatialGE_Data)
# data(tnbc_bassiouni)
#
# Basic gene expression plot
# STplot(tnbc, genes='CD3E', samples=1, ptsize=1)
#
# Multiple genes
# STplot(tnbc, genes=c('CD3E', 'MS4A1'), samples=1, ptsize=1)
#
# Cluster membership (requires STclust first)
# tnbc_clust <- STclust(tnbc, samples=1)
# STplot(tnbc, samples=1, ks='dtc', ptsize=1)
```

</div>

</div>

</div>
