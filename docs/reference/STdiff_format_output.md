# STdiff_format_output: Format results for user consumption

Formats DE results into a user-friendly structure with consistent column
naming, proper ordering, and optional filtering

## Usage

``` r
STdiff_format_output(
  results = NULL,
  select_cols = NULL,
  filter_sig = FALSE,
  adj_pval_thr = 0.05,
  order_by = "exp_adj_p_val",
  ascending = TRUE
)
```

## Arguments

- results:

  a list of data frames containing DE results per sample

- select_cols:

  optional vector of columns to include in output

- filter_sig:

  logical whether to filter for significant results only

- adj_pval_thr:

  p-value threshold for filtering (if filter_sig=TRUE)

- order_by:

  column name to order results by

- ascending:

  logical indicating sort order

## Value

formatted list of data frames ready for user consumption
