# Detect Input Source

Analyzes arguments to determine the Input Source strategy

## Usage

``` r
detect_input_source(rnacounts = NULL, spotcoords = NULL, samples = NULL)
```

## Arguments

- rnacounts:

  RNA count data (object, file paths, or directories)

- spotcoords:

  Spatial coordinates (object or file paths)

- samples:

  Sample names or metadata

## Value

list of InputSource objects, one per sample/source to process

## Details

Converts old-style detect_input() logic into modern modular InputSource
objects that can be dispatched to platform-specific ingestors.
