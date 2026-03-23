<div id="main" class="col-md-9" role="main">

# Constructor for InputSource

<div class="ref-description section level2">

Constructor for InputSource

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
new_input_source(
  type,
  format = NULL,
  rna = NULL,
  coords = NULL,
  samples = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   type:

    character string of input type (visium, xenium, generic, seurat,
    list)

-   format:

    character string detailed format (h5, mex, csv, etc)

-   rna:

    object or path for counts

-   coords:

    object or path for coordinates

-   samples:

    object or path for samples metadata

</div>

<div class="section level2">

## Value

list with class attribute set

</div>

</div>
