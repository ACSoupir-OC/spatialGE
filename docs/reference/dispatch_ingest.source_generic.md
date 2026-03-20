# Ingest Generic Data (DataFrames, Lists, Text Files)

Dispatches to appropriate reader for generic tabular data (CSV/TSV)

## Usage

``` r
# S3 method for class 'source_generic'
dispatch_ingest(source)
```

## Arguments

- source:

  InputSource object of type 'generic' or 'list'

## Value

list containing counts sparse matrix and coords dataframe

## Details

Reads delimited text files, converts to sparse matrices, and handles
basic barcode matching with coordinate files. Supports both single
samples and batch processing.
