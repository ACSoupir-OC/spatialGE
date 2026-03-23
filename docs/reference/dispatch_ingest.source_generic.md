<div id="main" class="col-md-9" role="main">

# Ingest Generic Data (DataFrames, Lists, Text Files)

<div class="ref-description section level2">

Dispatches to appropriate reader for generic tabular data (CSV/TSV)

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
# S3 method for class 'source_generic'
dispatch_ingest(source)
```

</div>

</div>

<div class="section level2">

## Arguments

-   source:

    InputSource object of type 'generic' or 'list'

</div>

<div class="section level2">

## Value

list containing counts sparse matrix and coords dataframe

</div>

<div class="section level2">

## Details

Reads delimited text files, converts to sparse matrices, and handles
basic barcode matching with coordinate files. Supports both single
samples and batch processing.

</div>

</div>
