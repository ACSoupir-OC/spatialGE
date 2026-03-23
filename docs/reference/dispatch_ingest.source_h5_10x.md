<div id="main" class="col-md-9" role="main">

# Ingest Generic H5 Details (Visium or Xenium or 10x)

<div class="ref-description section level2">

Dispatches to appropriate H5 reader based on file format heuristics

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
# S3 method for class 'source_h5_10x'
dispatch_ingest(source)
```

</div>

</div>

<div class="section level2">

## Arguments

-   source:

    InputSource object of type 'h5_10x'

</div>

<div class="section level2">

## Value

list containing counts sparse matrix and coords dataframe

</div>

<div class="section level2">

## Details

Uses coordinate file patterns and filename heuristics to distinguish
between Visium and Xenium H5 files. Falls back to Visium as default.

</div>

</div>
