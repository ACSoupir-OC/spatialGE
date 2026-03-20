# get_hier_clusters_dtc: Hierarchical clustering with DynamicTreeCut

Perform hierarchical clustering followed by DynamicTreeCut for adaptive
cluster detection

## Usage

``` r
get_hier_clusters_dtc(
  weighted_dists = NULL,
  ws = NULL,
  deepSplit = NULL,
  linkage = NULL
)
```

## Arguments

- weighted_dists:

  a list of distance matrices (NOT dist objects) for each spatial weight

- ws:

  a vector with spatial weights

- deepSplit:

  a logical or integer (1-4) for cutreeDynamic deepSplit parameter

- linkage:

  a string with the linkage method for hclust (e.g., 'ward.D2')

## Value

a list of data frames with spot/cell cluster assignments for each weight
