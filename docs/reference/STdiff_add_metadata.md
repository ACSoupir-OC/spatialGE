# STdiff_add_metadata: Add cluster and sample metadata to results

Enriches differential expression results with metadata from meta_dict
including original annotation names, sample information, and cluster
mappings

## Usage

``` r
STdiff_add_metadata(results = NULL, meta_dict = NULL, pairwise = FALSE)
```

## Arguments

- results:

  a list of data frames containing DE results per sample

- meta_dict:

  a data frame with columns (orig_annot, coded_annot) mapping coded
  annotations to original annotation names

- pairwise:

  logical indicating whether tests are pairwise or reference-based

## Value

enriched list of data frames with complete metadata columns
