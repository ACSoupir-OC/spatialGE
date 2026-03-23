<div id="main" class="col-md-9" role="main">

# STclust: Detect clusters of spots/cells (legacy)

<div class="ref-description section level2">

**\[superseded\]**

Perform unsupervised spatially-informed clustering on the spots/cells of
a ST sample

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STclust_legacy(
  x = NULL,
  samples = NULL,
  ws = 0.025,
  dist_metric = "euclidean",
  linkage = "ward.D2",
  ks = "dtc",
  topgenes = 2000,
  deepSplit = FALSE,
  cores = NULL,
  verbose = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    an STlist with normalized expression data

-   samples:

    a vector with strings or a vector with integers indicating the
    samples to run STclust

-   ws:

    a double (0-1) indicating the weight to be applied to spatial
    distances. Defaults to 0.025

-   dist_metric:

    the distance metric to be used. Defaults to 'euclidean'. Other
    options are the same as in `wordspace::dist.matrix`

-   linkage:

    the linkage method applied to hierarchical clustering. Passed to
    `hclust` and defaults to 'ward.D'

-   ks:

    the range of k values to assess. Defaults to `dtc`, meaning
    `cutreeDynamic` is applied

-   topgenes:

    the number of genes with highest spot-to-spot expression variation.
    The variance is calculated via `Seurat::FindVariableFeatures`.

-   deepSplit:

    a logical or integer (1-4), to be passed to `cutreeDynamic` and
    control cluster resolution

-   cores:

    an integer indicating the number of cores to use in parallelization
    (Unix only)

-   verbose:

    either logical or an integer (0, 1, or 2) to increase verbosity

</div>

<div class="section level2">

## Value

an STlist with cluster assignments

</div>

<div class="section level2">

## Details

**This is the legacy implementation from version 1.x.** Use `STclust`
for new projects.

The function takes an STlist and calculates euclidean distances between
cells or spots based on the x,y spatial locations, and the expression of
the top variable genes (`Seurat::FindVariableFeatures`). The resulting
distances are weighted by applying 1-`ws` to the gene expression
distances and `ws` to the spatial distances. Hierarchical clustering is
performed on the sum of the weighted distance matrices. The `STclust`
method allows for identification of tissue niches/domains that are
spatially cohesive.

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
  melanoma <- transform_data(melanoma)
  melanoma <- STclust(melanoma, ws=c(0, 0.025))
  STplot(melanoma, ws=0.025, samples='ST_mel1_rep2', ptsize=1)
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
