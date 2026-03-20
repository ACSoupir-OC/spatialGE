# Dispatch Ingest

Generic function to dispatch parsing logic based on InputSource class

## Usage

``` r
dispatch_ingest(source)
```

## Arguments

- source:

  InputSource object with type, format, rna, coords, and samples fields

## Value

list containing counts matrix and coords dataframe for the input

## Details

Uses R's S3 method dispatch to route to platform-specific ingestors
based on the InputSource object type.
