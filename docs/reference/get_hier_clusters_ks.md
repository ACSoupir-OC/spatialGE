# get_hier_clusters_ks: Hierarchical clustering with fixed k values

Perform hierarchical clustering followed by cutree for fixed number of
clusters

## Usage

``` r
get_hier_clusters_ks(
  weighted_dists = NULL,
  ws = NULL,
  ks = NULL,
  linkage = NULL
)
```

## Arguments

- weighted_dists:

  a list of distance matrices (NOT dist objects) for each spatial weight

- ws:

  a vector with spatial weights

- ks:

  a vector with k values for cluster detection

- linkage:

  a string with the linkage method for hclust (e.g., 'ward.D2')

## Value

a list of data frames with spot/cell cluster assignments for each weight
