<div id="main" class="col-md-9" role="main">

# STclust: Spatial clustering (modular implementation)

<div class="ref-description section level2">

**\[stable\]**

Perform unsupervised spatially-informed clustering on the spots/cells of
an ST sample

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STclust(
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

    a vector with strings or integers indicating samples to cluster

-   ws:

    a double (0-1) indicating weight for spatial distances (default:
    0.025)

-   dist_metric:

    distance metric to use (default: 'euclidean')

-   linkage:

    linkage method for hclust (default: 'ward.D2')

-   ks:

    clustering method: 'dtc' for DynamicTreeCut, or numeric vector for
    fixed k values

-   topgenes:

    number of top variable genes to use (default: 2000)

-   deepSplit:

    deepSplit parameter for cutreeDynamic (if ks='dtc')

-   cores:

    number of cores for parallelization

-   verbose:

    verbosity level (0, 1, or 2)

</div>

<div class="section level2">

## Value

an STlist with cluster assignments added to @spatial_meta

</div>

<div class="section level2">

## Details

The function calculates Euclidean distances between cells or spots based
on:

1.  Expression of top variable genes (Seurat's FindVariableFeatures)

2.  Spatial coordinates (x,y)

These distances are weighted and combined, then hierarchical clustering
is performed. Supports two clustering modes:

-   DTC (DynamicTreeCut): Adaptive cluster detection

-   Fixed k: User-specified number of clusters

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# \donttest{
# TNBC dataset (Bassiouni et al.)
tnbc_tmp = tempdir()
unlink(tnbc_tmp, recursive=TRUE)
dir.create(tnbc_tmp)
lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/tnbc_bassiouni.zip?download='
tryCatch({
  download.file(lk, destfile=file.path(tnbc_tmp, 'tnbc.zip'), mode='wb')
  unzip(file.path(tnbc_tmp, 'tnbc.zip'), exdir=tnbc_tmp)
  tnbc = STlist(rnacounts=list.files(file.path(tnbc_tmp, 'tnbc_bassiouni'), pattern='counts', full.names=TRUE)[1],
                spotcoords=list.files(file.path(tnbc_tmp, 'tnbc_bassiouni'), pattern='mapping', full.names=TRUE)[1],
                samples='TNBC')
  tnbc = transform_data(tnbc)
  tnbc = STclust(tnbc, ks=3, ws=0.025)
}, error = function(e) {
  message("Could not run example. Are you connected to the internet?")
})
#> Matching gene expression and coordinate data...
#> Could not run example. Are you connected to the internet?
# }
```

</div>

</div>

</div>
