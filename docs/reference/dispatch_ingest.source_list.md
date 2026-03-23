<div id="main" class="col-md-9" role="main">

# Ingest List of DataFrames

<div class="ref-description section level2">

Handles named list of dataframes with counts and coordinates

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
# S3 method for class 'source_list'
dispatch_ingest(source)
```

</div>

</div>

<div class="section level2">

## Arguments

-   source:

    InputSource object of type 'list'

</div>

<div class="section level2">

## Value

list containing counts list and coords list (one per sample)

</div>

<div class="section level2">

## Details

Processes named list inputs where each element contains gene expression
dataframes. Sorts and aligns count and coordinate dataframes.

</div>

</div>
