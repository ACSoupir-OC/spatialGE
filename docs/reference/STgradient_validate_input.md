# STgradient_validate_input: Validate STgradient inputs

Validate inputs for STgradient function

## Usage

``` r
STgradient_validate_input(
  x,
  samples,
  topgenes,
  annot,
  ref,
  exclude,
  min_nb,
  nb_dist_thr,
  cores,
  verbose
)
```

## Arguments

- x:

  STlist object

- samples:

  samples to process

- topgenes:

  number of top genes

- annot:

  annotation column

- ref:

  reference cluster

- exclude:

  optional exclude cluster

- min_nb:

  minimum neighbors

- nb_dist_thr:

  neighborhood distance threshold

- cores:

  number of cores

- verbose:

  verbosity

## Value

validated parameters list
