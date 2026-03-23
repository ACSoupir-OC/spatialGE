<div id="main" class="col-md-9" role="main">

# Internal MEX Directory Reader

<div class="ref-description section level2">

Reads Visium MEX (tar.gz or directory) Matrix Market format files

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
ingest_visium_mex_dir(
  mex_dir,
  coords_path = NULL,
  sample_name = NULL,
  original_root = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   mex_dir:

    Path to filtered_feature_bc_matrix or raw_feature_bc_matrix
    directory

-   coords_path:

    Optional path to tissue positions CSV file

-   sample_name:

    Optional sample identifier

-   original_root:

    Original root path for coordinate lookup (used with tar extraction)

</div>

<div class="section level2">

## Value

list with counts (sparse matrix) and coords (dataframe with barcode,
imagerow, imagecol)

</div>

<div class="section level2">

## Details

Extracts matrix.mtx, features.tsv, barcodes.tsv, and locates tissue
positions CSV. Handles both compressed and uncompressed formats.

</div>

</div>
