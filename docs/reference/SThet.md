<div id="main" class="col-md-9" role="main">

# SThet: Computes global spatial autocorrelation statistics on gene expression

<div class="ref-description section level2">

**\[stable\]**

Computes the global spatial autocorrelation statistics Moran's I and/or
Geary's C for a set of genes

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
SThet(
  x = NULL,
  genes = NULL,
  samples = NULL,
  method = "moran",
  k = NULL,
  overwrite = TRUE,
  cores = NULL,
  verbose = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    an STlist

-   genes:

    a vector of gene names to compute statistics

-   samples:

    the samples to compute statistics

-   method:

    The spatial statistic(s) to estimate. It can be set to 'moran',
    'geary' or both. Default is 'moran'

-   k:

    the number of neighbors to estimate weights. By default NULL,
    meaning that spatial weights will be estimated from Euclidean
    distances. If an positive integer is entered, then the faster k
    nearest-neighbors approach is used. Please keep in mind that
    estimates are not as accurate as when using the default
    distance-based method.

-   overwrite:

    logical indicating if previous statistics should be overwritten.
    Default to FALSE (do not overwrite)

-   cores:

    integer indicating the number of cores to use during
    parallelization. If NULL, the function uses half of the available
    cores at a maximum. The parallelization uses `parallel::mclapply`
    and works only in Unix systems

-   verbose:

    logical, whether to print text to console

</div>

<div class="section level2">

## Value

an STlist containing spatial statistics

</div>

<div class="section level2">

## Details

The function computes global spatial autocorrelation statistics (Moran's
I and/or Geary's C) for the requested genes and samples. Then
computation uses the package `spdep`. The calculated statistics are
stored in the STlist, which can be accessed with the `get_gene_meta`
function. For visual comparative analysis, the function `compare_SThet`
can be used afterwards.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# Example with TNBC data showing heterogeneity calculation
# Load TNBC test data
data_dir = system.file("tests/testthat/data/tnbc_bassiouni", package="spatialGE")
count_files <- list.files(data_dir, pattern='counts', full.names=TRUE)
coord_files <- list.files(data_dir, pattern='mapping', full.names=TRUE)
clin_file <- list.files(data_dir, pattern='clinical', full.names=TRUE)

# Create STlist and transform
tnb <- STlist(rnacounts=count_files, spotcoords=coord_files, samples=clin_file)
#> Matching gene expression and coordinate data...
#> Error in counts_df_list[[name_i]]: attempt to select less than one element in get1index
tnb <- transform_data(tnb)
#> Error: object 'tnb' not found

# Calculate spatial autocorrelation for genes of interest
tnb <- SThet(tnb, genes=c('ESR1', 'PGR', 'ERBB2'), method='moran')
#> SThet started.
#> Error: object 'tnb' not found

# Extract and view heterogeneity statistics
meta <- get_gene_meta(tnb, sthet_only=TRUE)
#> Error: object 'tnb' not found
print(meta)
#> Error: object 'meta' not found
```

</div>

</div>

</div>
