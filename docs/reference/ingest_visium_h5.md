<div id="main" class="col-md-9" role="main">

# Internal H5 Reader

<div class="ref-description section level2">

Reads 10X Visium H5 files and extracts count matrix and coordinates

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
ingest_visium_h5(h5_path, coords_path = NULL, sample_name = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   h5_path:

    Path to filtered_feature_bc_matrix.h5 or raw_feature_bc_matrix.h5

-   coords_path:

    Optional path to tissue positions CSV file

-   sample_name:

    Optional sample identifier

</div>

<div class="section level2">

## Value

list with counts (sparse matrix) and coords (dataframe with barcode,
imagerow, imagecol)

</div>

<div class="section level2">

## Details

Extracts sparse matrix components from H5 structure, locates tissue
positions CSV files, and aligns barcodes between counts and coordinates.

</div>

</div>
