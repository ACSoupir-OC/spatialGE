<div id="main" class="col-md-9" role="main">

# Ingest Seurat Object

<div class="ref-description section level2">

Converts Seurat spatial object to STlist format

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
# S3 method for class 'source_seurat'
dispatch_ingest(source)
```

</div>

</div>

<div class="section level2">

## Arguments

-   source:

    InputSource object of type 'seurat'

</div>

<div class="section level2">

## Value

list containing counts list and coords list (one per slice)

</div>

<div class="section level2">

## Details

Extracts counts from Spatial assay and coordinates from images slot,
supporting multiple slices from a single Seurat object.

</div>

</div>
