<div id="main" class="col-md-9" role="main">

# STlist: Creation of STlist objects for spatial transcriptomics analysis

<div class="ref-description section level2">

**\[stable\]**

Creates an STlist object from one or multiple spatial transcriptomic
samples.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STlist(
  rnacounts = NULL,
  spotcoords = NULL,
  samples = NULL,
  cores = NULL,
  verbose = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   rnacounts:

    the count data which can be provided in various formats (see
    `STList_legacy` for detailed legacy format descriptions, or new
    modular documentation).

-   spotcoords:

    the cell/spot coordinates. Not required if inputs are Visium or
    Xenium.

-   samples:

    the sample names/IDs and (optionally) metadata.

-   cores:

    integer indicating the number of cores to use during
    parallelization.

-   verbose:

    logical, whether to print text to console.

</div>

<div class="section level2">

## Value

an STlist object containing the counts and coordinates, and optionally
sample metadata.

</div>

<div class="section level2">

## Details

Objects of the S4 class STlist are the starting point of analyses in
**`spatialGE`**.

**Note:** This is the new modular implementation of `STList` created in
version 2.0.0. For reproducibility, the previous version is available as
`STList_legacy`.

The STlist contains data from one or multiple samples (i.e., tissue
slices), and results from most `spatialGE`'s functions are stored within
the object.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`STList_legacy`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# Example using TNBC test data included with spatialGE
if (requireNamespace("spatialGE", quietly = TRUE)) {
  data_dir <- system.file("tests/testthat/data/tnbc_bassiouni", package="spatialGE")
  if (data_dir != "") {
    # Get file paths for first sample
    sample_dir <- list.dirs(data_dir, recursive = FALSE)[1]
    count_files <- list.files(sample_dir, pattern="\\.h5$", full.names = TRUE)
    coord_dir <- file.path(sample_dir, "spatial")
    coord_files <- list.files(coord_dir, pattern="tissue_positions", full.names = TRUE)
    clin_file <- list.files(data_dir, pattern="clinical", full.names = TRUE)
    
    # Create STlist object
    tnbc <- STlist(rnacounts = count_files, 
                   spotcoords = coord_files, 
                   samples = clin_file,
                   verbose = TRUE)
    
    # Inspect object
    tnbc
    summary(tnbc)
  }
}
```

</div>

</div>

</div>
