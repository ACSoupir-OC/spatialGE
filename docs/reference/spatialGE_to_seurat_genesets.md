# Create Seurat-compatible gene sets from spatialGE results

**\[stable\]**

Converts STenrich or STgradient results to gene set format for Seurat
GSEA.

## Usage

``` r
spatialGE_to_seurat_genesets(st_result, pval.thr = 0.05, return.type = "list")
```

## Arguments

- st_result:

  results from STenrich or STgradient analysis

- pval.thr:

  p-value threshold for significance (default: 0.05)

- return.type:

  'list' (gene sets) or 'data.frame' (summary table)

## Value

gene sets in Seurat-compatible format

## Details

This function extracts significant gene sets from spatialGE enrichment
analysis and formats them for use with Seurat's AddModuleScore or
similar functions.

## See also

[`STenrich`](https://acsoupir-oc.github.io/spatialGE/reference/STenrich.md),
[`STgradient`](https://acsoupir-oc.github.io/spatialGE/reference/STgradient.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Run enrichment analysis
enrich_result <- STenrich(st_obj, gene_sets = pathway_list)

# Convert to Seurat gene sets
gene_sets <- spatialGE_to_seurat_genesets(enrich_result, pval.thr = 0.05)

# Add module scores
seurat_obj <- Seurat::AddModuleScore(seurat_obj, features = gene_sets)
} # }
```
