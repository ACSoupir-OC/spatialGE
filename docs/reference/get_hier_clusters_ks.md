<div id="main" class="col-md-9" role="main">

# get_hier_clusters_ks: Hierarchical clustering with fixed k values

<div class="ref-description section level2">

Perform hierarchical clustering followed by cutree for fixed number of
clusters

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
get_hier_clusters_ks(
  weighted_dists = NULL,
  ws = NULL,
  ks = NULL,
  linkage = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   weighted_dists:

    a list of distance matrices (NOT dist objects) for each spatial
    weight

-   ws:

    a vector with spatial weights

-   ks:

    a vector with k values for cluster detection

-   linkage:

    a string with the linkage method for hclust (e.g., 'ward.D2')

</div>

<div class="section level2">

## Value

a list of data frames with spot/cell cluster assignments for each weight

</div>

</div>
