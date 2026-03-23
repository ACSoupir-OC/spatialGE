<div id="main" class="col-md-9" role="main">

# Create Seurat-compatible gene sets from spatialGE results

<div class="ref-description section level2">

**\[stable\]**

Converts STenrich or STgradient results to gene set format for Seurat
GSEA.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
spatialGE_to_seurat_genesets(st_result, pval.thr = 0.05, return.type = "list")
```

</div>

</div>

<div class="section level2">

## Arguments

-   st_result:

    results from STenrich or STgradient analysis

-   pval.thr:

    p-value threshold for significance (default: 0.05)

-   return.type:

    'list' (gene sets) or 'data.frame' (summary table)

</div>

<div class="section level2">

## Value

gene sets in Seurat-compatible format

</div>

<div class="section level2">

## Details

This function extracts significant gene sets from spatialGE enrichment
analysis and formats them for use with Seurat's AddModuleScore or
similar functions.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`STenrich`, `STgradient`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

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

</div>

</div>

</div>
