# STclust_calculate_distances: Calculate scaled distance matrices

Calculate and scale gene expression and spatial distance matrices for a
single sample

## Usage

``` r
STclust_calculate_distances(trcounts_df, coord_dat, dist_metric)
```

## Arguments

- trcounts_df:

  a filtered expression matrix for one sample

- coord_dat:

  a data frame with spot coordinates for one sample

- dist_metric:

  distance metric to use

## Value

a list with scaled distance matrices (expression + spatial)
