# Ingest Generic H5 Details (Visium or Xenium or 10x)

Dispatches to appropriate H5 reader based on file format heuristics

## Usage

``` r
# S3 method for class 'source_h5_10x'
dispatch_ingest(source)
```

## Arguments

- source:

  InputSource object of type 'h5_10x'

## Value

list containing counts sparse matrix and coords dataframe

## Details

Uses coordinate file patterns and filename heuristics to distinguish
between Visium and Xenium H5 files. Falls back to Visium as default.
