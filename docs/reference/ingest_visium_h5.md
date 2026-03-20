# Internal H5 Reader

Reads 10X Visium H5 files and extracts count matrix and coordinates

## Usage

``` r
ingest_visium_h5(h5_path, coords_path = NULL, sample_name = NULL)
```

## Arguments

- h5_path:

  Path to filtered_feature_bc_matrix.h5 or raw_feature_bc_matrix.h5

- coords_path:

  Optional path to tissue positions CSV file

- sample_name:

  Optional sample identifier

## Value

list with counts (sparse matrix) and coords (dataframe with barcode,
imagerow, imagecol)

## Details

Extracts sparse matrix components from H5 structure, locates tissue
positions CSV files, and aligns barcodes between counts and coordinates.
