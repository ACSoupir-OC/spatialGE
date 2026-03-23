<div id="main" class="col-md-9" role="main">

# STgradient_summarize_distances: Summarize distances from reference

<div class="ref-description section level2">

Calculate min or average distance from reference for non-reference spots

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STgradient_summarize_distances(
  dist_tmp,
  nonref_tmp,
  ref_tmp,
  nbs_keep,
  distsumm,
  limit
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   dist_tmp:

    distance matrix

-   nonref_tmp:

    non-reference barcodes

-   ref_tmp:

    reference barcodes

-   nbs_keep:

    reference barcodes to keep

-   distsumm:

    distance summary metric ('min' or 'avg')

-   limit:

    optional distance limit

</div>

<div class="section level2">

## Value

data frame with barcode and dist2ref

</div>

</div>
