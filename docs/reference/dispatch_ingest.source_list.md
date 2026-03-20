# Ingest List of DataFrames

Handles named list of dataframes with counts and coordinates

## Usage

``` r
# S3 method for class 'source_list'
dispatch_ingest(source)
```

## Arguments

- source:

  InputSource object of type 'list'

## Value

list containing counts list and coords list (one per sample)

## Details

Processes named list inputs where each element contains gene expression
dataframes. Sorts and aligns count and coordinate dataframes.
