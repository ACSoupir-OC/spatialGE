# calculate_dist_matrices: Calculate and scale distance matrices

Calculate Euclidean distance matrices for gene expression and spatial
coordinates, then scale both to 0,1 range

## Usage

``` r
calculate_dist_matrices(expr_dat = NULL, coord_dat = NULL, dist_metric = NULL)
```

## Arguments

- expr_dat:

  a sparse matrix with gene expression (genes in rows, spots/cells in
  columns)

- coord_dat:

  a data frame with columns: 'libname', 'xpos', 'ypos'

- dist_metric:

  character string indicating distance metric (e.g., 'euclidean')

## Value

a list with two matrices: `scale_exp` (scaled expression distances) and
`scale_coord` (scaled coordinate distances)
