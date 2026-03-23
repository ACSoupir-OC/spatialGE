<div id="main" class="col-md-9" role="main">

# STlist (Legacy): Creation of STlist objects for spatial transcriptomics analysis

<div class="ref-description section level2">

**\[superseded\]**

Creates an STlist object from one or multiple spatial transcriptomic
samples.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STList_legacy(
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

    the count data.

-   spotcoords:

    the cell/spot coordinates.

-   samples:

    the sample names/IDs.

-   cores:

    integer indicating the number of cores to use.

-   verbose:

    logical, whether to print text to console.

</div>

<div class="section level2">

## Value

an STlist object.

</div>

<div class="section level2">

## Details

**This function is superseded by the new `STList` implementation in
version 2.0.0.** It is preserved here for reproducibility of older
workflows. Please use `STList` for new projects.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`STList`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# \donttest{
# Using included melanoma example (Thrane et al.)
# Download example data set from spatialGE_Data
thrane_tmp = tempdir()
unlink(thrane_tmp, recursive=TRUE)
dir.create(thrane_tmp)
lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
tryCatch({ # In case data is not available from network
  download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
  #' zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
  unzip(zipfile=zip_tmp, exdir=thrane_tmp)
  # Generate the file paths to be passed to the STlist function
  count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
                            full.names=TRUE, pattern='counts')
  coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
                            full.names=TRUE, pattern='mapping')
  clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
                          full.names=TRUE, pattern='clinical')
  # Create STlist
  library('spatialGE')
  melanoma <- STlist(rnacounts=count_files,
                     spotcoords=coord_files,
                     samples=clin_file)
  melanoma
}, error = function(e) {
  message("Could not run example. Are you connected to the internet?")
  return(NULL)
})
#> Could not run example. Are you connected to the internet?
#> NULL
# }
```

</div>

</div>

</div>
