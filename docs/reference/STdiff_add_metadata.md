<div id="main" class="col-md-9" role="main">

# STdiff_add_metadata: Add cluster and sample metadata to results

<div class="ref-description section level2">

Enriches differential expression results with metadata from meta_dict
including original annotation names, sample information, and cluster
mappings

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
STdiff_add_metadata(results = NULL, meta_dict = NULL, pairwise = FALSE)
```

</div>

</div>

<div class="section level2">

## Arguments

-   results:

    a list of data frames containing DE results per sample

-   meta_dict:

    a data frame with columns (orig_annot, coded_annot) mapping coded
    annotations to original annotation names

-   pairwise:

    logical indicating whether tests are pairwise or reference-based

</div>

<div class="section level2">

## Value

enriched list of data frames with complete metadata columns

</div>

</div>
