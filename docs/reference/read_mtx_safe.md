# Safe Matrix Reader

Robustly reads Matrix market files, handling .gz via temp expansion

## Usage

``` r
read_mtx_safe(path)
```

## Arguments

- path:

  Path to matrix.mtx or matrix.mtx.gz file

## Value

CsparseMatrix from Matrix package

## Details

Decompresses gzip files to temporary location before reading, ensuring
proper handling of compressed Matrix Market format files.
