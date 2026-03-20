# Internal Xenium MEX Directory Reader

Reads Xenium MEX (tar.gz or directory) Matrix Market format files

## Usage

``` r
ingest_xenium_mex_dir(
  mex_dir,
  coords_path = NULL,
  sample_name = NULL,
  original_root = NULL
)
```

## Arguments

- mex_dir:

  Path to cell_feature_matrix directory

- coords_path:

  Optional path to cells.parquet or cells.csv file

- sample_name:

  Optional sample identifier

- original_root:

  Original root path for coordinate lookup (used with tar extraction)

## Value

list with counts (sparse matrix) and coords (dataframe with barcode,
ypos, xpos)

## Details

Extracts matrix.mtx, features.tsv, barcodes.tsv, and locates
cells.parquet or cells.csv coordinate files. Handles both compressed and
uncompressed formats.
