<div id="main" class="col-md-9" role="main">

# Dispatch Ingest

<div class="ref-description section level2">

Generic function to dispatch parsing logic based on InputSource class

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
dispatch_ingest(source)
```

</div>

</div>

<div class="section level2">

## Arguments

-   source:

    InputSource object with type, format, rna, coords, and samples
    fields

</div>

<div class="section level2">

## Value

list containing counts matrix and coords dataframe for the input

</div>

<div class="section level2">

## Details

Uses R's S3 method dispatch to route to platform-specific ingestors
based on the InputSource object type.

</div>

</div>
