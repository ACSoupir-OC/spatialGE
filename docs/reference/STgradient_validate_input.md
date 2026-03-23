<div id="main" class="col-md-9" role="main">

# STgradient_validate_input: Validate STgradient inputs

<div class="ref-description section level2">

Validate inputs for STgradient function

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

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

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    STlist object

-   samples:

    samples to process

-   topgenes:

    number of top genes

-   annot:

    annotation column

-   ref:

    reference cluster

-   exclude:

    optional exclude cluster

-   min_nb:

    minimum neighbors

-   nb_dist_thr:

    neighborhood distance threshold

-   cores:

    number of cores

-   verbose:

    verbosity

</div>

<div class="section level2">

## Value

validated parameters list

</div>

</div>
