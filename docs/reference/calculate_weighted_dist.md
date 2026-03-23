<div id="main" class="col-md-9" role="main">

# calculate_weighted_dist: Calculate weighted distance matrices

<div class="ref-description section level2">

Combine scaled expression and spatial distance matrices using
user-specified weights

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
calculate_weighted_dist(scaled_dists = NULL, ws = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   scaled_dists:

    a list with two matrices: scaled expression and spatial distance
    matrices

-   ws:

    vector of spatial weights (0-1) to apply to spatial distances

</div>

<div class="section level2">

## Value

a list of weighted distance matrices, one for each ws value

</div>

</div>
