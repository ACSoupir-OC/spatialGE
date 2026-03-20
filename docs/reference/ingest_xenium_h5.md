# Internal Xenium H5 Reader

Reads 10X Xenium H5 files and extracts count matrix and coordinates

## Usage

``` r
ingest_xenium_h5(h5_path, coords_path = NULL, sample_name = NULL)
```

## Arguments

- h5_path:

  Path to cell_feature_matrix.h5

- coords_path:

  Optional path to cells.parquet or cells.csv file

- sample_name:

  Optional sample identifier

## Value

list with counts (sparse matrix) and coords (dataframe with barcode,
ypos, xpos)

## Details

Extracts sparse matrix components from H5 structure, locates
cells.parquet or cells.csv coordinate files, and maps cell_id to barcode
format.
