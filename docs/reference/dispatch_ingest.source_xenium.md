# Ingest Xenium Data

Dispatches to appropriate Xenium reader based on file format (H5 or MEX)

## Usage

``` r
# S3 method for class 'source_xenium'
dispatch_ingest(source)
```

## Arguments

- source:

  InputSource object of type 'xenium'

## Value

list containing counts sparse matrix and coords dataframe
