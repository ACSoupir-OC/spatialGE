<div id="main" class="col-md-9" role="main">

# summarize_STlist: Generates a data frame with summary statistics

<div class="ref-description section level2">

Produces a data frame with counts per gene and counts per ROI/spot/cell

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
summarize_STlist(x = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    an STlist

</div>

<div class="section level2">

## Value

a data frame

</div>

<div class="section level2">

## Details

The function creates a table with counts per gene and counts per region
of interest (ROI), spot, or cell in the samples stored in the STlist

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
                     samples=clin_file) # Only first two samples
  summarize_STlist(melanoma)
}, error = function(e) {
  message("Could not run example. Are you connected to the internet?.")
  return(NULL)
})
#> Could not run example. Are you connected to the internet?.
#> NULL
# }
```

</div>

</div>

</div>
