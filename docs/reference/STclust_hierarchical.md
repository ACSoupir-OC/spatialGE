# STclust_hierarchical: Core hierarchical clustering logic

Performs hierarchical clustering on weighted distance matrices using
either DTC or fixed k

## Usage

``` r
STclust_hierarchical(weighted_dists, ws, ks, linkage, deepSplit, verbose)
```

## Arguments

- weighted_dists:

  list of weighted distance matrices for one sample (one per ws value)

- ws:

  vector of spatial weights

- ks:

  clustering method: 'dtc' for DynamicTreeCut or numeric vector for
  fixed k values

- linkage:

  linkage method for hclust (e.g., 'ward.D2')

- deepSplit:

  deepSplit parameter for cutreeDynamic (if ks='dtc')

- verbose:

  verbosity level

## Value

list of data frames with cluster assignments for each ws value
