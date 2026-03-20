# Ingest Visium Data

Dispatches to appropriate Visium reader based on file format (H5 or MEX)

## Usage

``` r
# S3 method for class 'source_visium'
dispatch_ingest(source)
```

## Arguments

- source:

  InputSource object of type 'visium'

## Value

list containing counts sparse matrix and coords dataframe
