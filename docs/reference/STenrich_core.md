<div id="main" class="col-md-9" role="main">

# STenrich_core: Core workflow for spatial enrichment analysis

<div class="ref-description section level2">

Main workflow function that orchestrates the enrichment analysis

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STenrich_core(
  x,
  samples,
  gene_sets,
  score_type = "avg",
  reps = 1000,
  annot = NULL,
  domain = NULL,
  num_sds = 1,
  min_units = 20,
  min_genes = 5,
  pval_adj_method = "BH",
  seed = 12345,
  cores = NULL,
  verbose = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    an STlist with transformed gene expression

-   samples:

    a vector with sample names or indexes to run analysis

-   gene_sets:

    a named list of gene sets to test

-   score_type:

    Controls how gene set expression is calculated. Options: 'avg' or
    'gsva'

-   reps:

    the number of random samples to be extracted. Default is 1000
    replicates

-   annot:

    name of the annotation within `x@spatial_meta` containing spot/cell
    categories

-   domain:

    the domain to restrict the analysis

-   num_sds:

    number of standard deviations to set minimum gene set expression
    threshold (default: 1)

-   min_units:

    Minimum number of spots with high expression (default: 20)

-   min_genes:

    the minimum number of genes of a gene set present in data (default:
    5)

-   pval_adj_method:

    the method for multiple comparison adjustment (default: 'BH')

-   seed:

    the seed number for random sampling (default: 12345)

-   cores:

    the number of cores for parallelization. If NULL, auto-detected

-   verbose:

    logical, whether to print text to console

</div>

<div class="section level2">

## Value

list of data frames with p-values and adjusted p-values per sample

</div>

<div class="section level2">

## Details

This function performs the core workflow:

1.  Prepare data (extract tissue spots, coordinates, match gene sets)

2.  Calculate gene set expression (average or GSVA)

3.  Perform permutation testing

4.  Format results with adjusted p-values

</div>

</div>
