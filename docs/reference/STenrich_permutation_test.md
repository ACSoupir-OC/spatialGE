<div id="main" class="col-md-9" role="main">

# STenrich_permutation_test: Perform permutation testing

<div class="ref-description section level2">

Calculate observed distances and perform random permutations

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STenrich_permutation_test(
  result_df,
  coords_df,
  combo,
  pw_genes,
  samples,
  gene_sets,
  num_sds,
  min_units,
  reps,
  seed,
  cores,
  verbose
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   result_df:

    list of data frames with expression scores

-   coords_df:

    list of coordinate matrices for each sample

-   combo:

    data frame with combinations of samples and gene sets

-   pw_genes:

    list of available genes per gene set

-   samples:

    a vector with sample names to run analysis

-   gene_sets:

    a named list of gene sets to test

-   num_sds:

    number of standard deviations to set minimum gene set expression
    threshold (default: 1)

-   min_units:

    Minimum number of spots with high expression (default: 20)

-   reps:

    the number of random samples to be extracted (default: 1000)

-   seed:

    the seed number for random sampling (default: 12345)

-   cores:

    the number of cores for parallelization

-   verbose:

    verbosity level

</div>

<div class="section level2">

## Value

list of data frames with p-values per sample

</div>

</div>
