# STlist: Creation of STlist objects for spatial transcriptomics analysis

**\[stable\]**

Creates an STlist object from one or multiple spatial transcriptomic
samples.

## Usage

``` r
STlist(
  rnacounts = NULL,
  spotcoords = NULL,
  samples = NULL,
  cores = NULL,
  verbose = TRUE
)
```

## Arguments

- rnacounts:

  the count data which can be provided in various formats (see
  [`STList_legacy`](https://acsoupir-oc.github.io/spatialGE/reference/STList_legacy.md)
  for detailed legacy format descriptions, or new modular
  documentation).

- spotcoords:

  the cell/spot coordinates. Not required if inputs are Visium or
  Xenium.

- samples:

  the sample names/IDs and (optionally) metadata.

- cores:

  integer indicating the number of cores to use during parallelization.

- verbose:

  logical, whether to print text to console.

## Value

an STlist object containing the counts and coordinates, and optionally
sample metadata.

## Details

Objects of the S4 class STlist are the starting point of analyses in
**`spatialGE`**.

**Note:** This is the new modular implementation of `STList` created in
version 2.0.0. For reproducibility, the previous version is available as
`STList_legacy`.

The STlist contains data from one or multiple samples (i.e., tissue
slices), and results from most `spatialGE`'s functions are stored within
the object.

## See also

[`STList_legacy`](https://acsoupir-oc.github.io/spatialGE/reference/STList_legacy.md)
