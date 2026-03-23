<div id="main" class="col-md-9" role="main">

# calculate_dist_matrices: Calculate and scale distance matrices

<div class="ref-description section level2">

Calculate Euclidean distance matrices for gene expression and spatial
coordinates, then scale both to 0,1 range

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
calculate_dist_matrices(expr_dat = NULL, coord_dat = NULL, dist_metric = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   expr_dat:

    a sparse matrix with gene expression (genes in rows, spots/cells in
    columns)

-   coord_dat:

    a data frame with columns: 'libname', 'xpos', 'ypos'

-   dist_metric:

    character string indicating distance metric (e.g., 'euclidean')

</div>

<div class="section level2">

## Value

a list with two matrices: `scale_exp` (scaled expression distances) and
`scale_coord` (scaled coordinate distances)

</div>

</div>
