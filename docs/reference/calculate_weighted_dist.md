# calculate_weighted_dist: Calculate weighted distance matrices

Combine scaled expression and spatial distance matrices using
user-specified weights

## Usage

``` r
calculate_weighted_dist(scaled_dists = NULL, ws = NULL)
```

## Arguments

- scaled_dists:

  a list with two matrices: scaled expression and spatial distance
  matrices

- ws:

  vector of spatial weights (0-1) to apply to spatial distances

## Value

a list of weighted distance matrices, one for each ws value
