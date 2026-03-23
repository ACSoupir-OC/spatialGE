<div id="main" class="col-md-9" role="main">

# Detect Input Source

<div class="ref-description section level2">

Analyzes arguments to determine the Input Source strategy

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
detect_input_source(rnacounts = NULL, spotcoords = NULL, samples = NULL)
```

</div>

</div>

<div class="section level2">

## Arguments

-   rnacounts:

    RNA count data (object, file paths, or directories)

-   spotcoords:

    Spatial coordinates (object or file paths)

-   samples:

    Sample names or metadata

</div>

<div class="section level2">

## Value

list of InputSource objects, one per sample/source to process

</div>

<div class="section level2">

## Details

Converts old-style detect_input() logic into modern modular InputSource
objects that can be dispatched to platform-specific ingestors.

</div>

</div>
