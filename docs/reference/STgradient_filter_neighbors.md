# STgradient_filter_neighbors: Filter reference spots by neighbor count

Identify reference spots with sufficient neighbors

## Usage

``` r
STgradient_filter_neighbors(dist_tmp, ref_tmp, min_nb, nb_dist_thr)
```

## Arguments

- dist_tmp:

  distance matrix

- ref_tmp:

  reference barcodes

- min_nb:

  minimum neighbors required

- nb_dist_thr:

  neighborhood distance threshold

## Value

vector of reference barcodes to keep
