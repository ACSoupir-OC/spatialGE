<div id="main" class="col-md-9" role="main">

# STclust_hierarchical: Core hierarchical clustering logic

<div class="ref-description section level2">

Performs hierarchical clustering on weighted distance matrices using
either DTC or fixed k

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STclust_hierarchical(weighted_dists, ws, ks, linkage, deepSplit, verbose)
```

</div>

</div>

<div class="section level2">

## Arguments

-   weighted_dists:

    list of weighted distance matrices for one sample (one per ws value)

-   ws:

    vector of spatial weights

-   ks:

    clustering method: 'dtc' for DynamicTreeCut or numeric vector for
    fixed k values

-   linkage:

    linkage method for hclust (e.g., 'ward.D2')

-   deepSplit:

    deepSplit parameter for cutreeDynamic (if ks='dtc')

-   verbose:

    verbosity level

</div>

<div class="section level2">

## Value

list of data frames with cluster assignments for each ws value

</div>

</div>
