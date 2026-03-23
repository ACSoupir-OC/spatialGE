<div id="main" class="col-md-9" role="main">

# STgradient_filter_neighbors: Filter reference spots by neighbor count

<div class="ref-description section level2">

Identify reference spots with sufficient neighbors

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STgradient_filter_neighbors(dist_tmp, ref_tmp, min_nb, nb_dist_thr)
```

</div>

</div>

<div class="section level2">

## Arguments

-   dist_tmp:

    distance matrix

-   ref_tmp:

    reference barcodes

-   min_nb:

    minimum neighbors required

-   nb_dist_thr:

    neighborhood distance threshold

</div>

<div class="section level2">

## Value

vector of reference barcodes to keep

</div>

</div>
