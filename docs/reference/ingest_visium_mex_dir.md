# Internal MEX Directory Reader

Reads Visium MEX (tar.gz or directory) Matrix Market format files

## Usage

``` r
ingest_visium_mex_dir(
  mex_dir,
  coords_path = NULL,
  sample_name = NULL,
  original_root = NULL
)
```

## Arguments

- mex_dir:

  Path to filtered_feature_bc_matrix or raw_feature_bc_matrix directory

- coords_path:

  Optional path to tissue positions CSV file

- sample_name:

  Optional sample identifier

- original_root:

  Original root path for coordinate lookup (used with tar extraction)

## Value

list with counts (sparse matrix) and coords (dataframe with barcode,
imagerow, imagecol)

## Details

Extracts matrix.mtx, features.tsv, barcodes.tsv, and locates tissue
positions CSV. Handles both compressed and uncompressed formats.
