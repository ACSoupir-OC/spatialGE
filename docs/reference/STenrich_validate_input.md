# STenrich_validate_input: Validate input parameters

Check all input parameters and return validated values

## Usage

``` r
STenrich_validate_input(
  x,
  samples,
  gene_sets,
  score_type,
  annot,
  domain,
  num_sds,
  min_units,
  min_genes,
  pval_adj_method
)
```

## Arguments

- x:

  an STlist with transformed gene expression

- samples:

  a vector with sample names or indexes to run analysis

- gene_sets:

  a named list of gene sets to test

- score_type:

  Controls how gene set expression is calculated. Options: 'avg' or
  'gsva'

- annot:

  name of the annotation within `x@spatial_meta` containing spot/cell
  categories

- domain:

  the domain to restrict the analysis

- num_sds:

  number of standard deviations to set minimum gene set expression
  threshold (default: 1)

- min_units:

  Minimum number of spots with high expression (default: 20)

- min_genes:

  the minimum number of genes of a gene set present in data (default: 5)

- pval_adj_method:

  the method for multiple comparison adjustment (default: 'BH')

## Value

list with validated parameters
