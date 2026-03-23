<div id="main" class="col-md-9" role="main">

# STenrich_prepare_data: Prepare data for enrichment analysis

<div class="ref-description section level2">

Extract tissue spots, coordinates, and gene sets from STlist

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STenrich_prepare_data(x, samples, gene_sets, annot, domain, min_units)
```

</div>

</div>

<div class="section level2">

## Arguments

-   x:

    an STlist with transformed gene expression

-   samples:

    a vector with sample names to run analysis

-   gene_sets:

    a named list of gene sets to test

-   annot:

    name of the annotation within `x@spatial_meta` containing spot/cell
    categories

-   domain:

    the domain to restrict the analysis

-   min_units:

    Minimum number of spots with high expression (default: 20)

</div>

<div class="section level2">

## Value

list with prepared data including tissue_spots, coords_df, combo,
pw_genes

</div>

</div>
