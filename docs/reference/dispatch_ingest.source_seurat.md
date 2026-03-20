# Ingest Seurat Object

Converts Seurat spatial object to STlist format

## Usage

``` r
# S3 method for class 'source_seurat'
dispatch_ingest(source)
```

## Arguments

- source:

  InputSource object of type 'seurat'

## Value

list containing counts list and coords list (one per slice)

## Details

Extracts counts from Spatial assay and coordinates from images slot,
supporting multiple slices from a single Seurat object.
