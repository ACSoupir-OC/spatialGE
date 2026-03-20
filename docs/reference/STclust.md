# STclust: Spatial clustering (modular implementation)

Perform unsupervised spatially-informed clustering on the spots/cells of
an ST sample

## Usage

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

## Arguments

- x:

  an STlist with normalized expression data

- samples:

  a vector with strings or integers indicating samples to cluster

- ws:

  a double (0-1) indicating weight for spatial distances (default:
  0.025)

- dist_metric:

  distance metric to use (default: 'euclidean')

- linkage:

  linkage method for hclust (default: 'ward.D2')

- ks:

  clustering method: 'dtc' for DynamicTreeCut, or numeric vector for
  fixed k values

- topgenes:

  number of top variable genes to use (default: 2000)

- deepSplit:

  deepSplit parameter for cutreeDynamic (if ks='dtc')

- cores:

  number of cores for parallelization

- verbose:

  verbosity level (0, 1, or 2)

## Value

an STlist with cluster assignments added to @spatial_meta

## Details

The function calculates Euclidean distances between cells or spots based
on:

1.  Expression of top variable genes (Seurat's FindVariableFeatures)

2.  Spatial coordinates (x,y)

These distances are weighted and combined, then hierarchical clustering
is performed. Supports two clustering modes:

- DTC (DynamicTreeCut): Adaptive cluster detection

- Fixed k: User-specified number of clusters

## Examples

``` r
# \donttest{
# Using included melanoma example (Thrane et al.)
thrane_tmp = tempdir()
unlink(thrane_tmp, recursive=TRUE)
dir.create(thrane_tmp)
lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
tryCatch({
  download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
  unzip(file=paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir=thrane_tmp)
  count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='counts')
  coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='mapping')
  clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='clinical')
  melanoma <- STlist(rnacounts=count_files[1], spotcoords=coord_files[1], samples=clin_file)
  melanoma <- transform_data(melanoma)
  melanoma <- STclust(melanoma, ks=2:3, ws=0.025)
  STplot(melanoma, ws=0.025, samples='ST_mel1_rep2', ptsize=1)
}, error = function(e) {
  message("Could not run example. Are you connected to the internet?")
  return(NULL)
})
# }
```
