# Constructor for InputSource

Constructor for InputSource

## Usage

``` r
new_input_source(
  type,
  format = NULL,
  rna = NULL,
  coords = NULL,
  samples = NULL
)
```

## Arguments

- type:

  character string of input type (visium, xenium, generic, seurat, list)

- format:

  character string detailed format (h5, mex, csv, etc)

- rna:

  object or path for counts

- coords:

  object or path for coordinates

- samples:

  object or path for samples metadata

## Value

list with class attribute set
