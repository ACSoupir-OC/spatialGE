# STgradient_summarize_distances: Summarize distances from reference

Calculate min or average distance from reference for non-reference spots

## Usage

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

## Arguments

- dist_tmp:

  distance matrix

- nonref_tmp:

  non-reference barcodes

- ref_tmp:

  reference barcodes

- nbs_keep:

  reference barcodes to keep

- distsumm:

  distance summary metric ('min' or 'avg')

- limit:

  optional distance limit

## Value

data frame with barcode and dist2ref
