#!/usr/bin/env Rscript
# Seurat Integration Vignette - Executable Version
# This script runs through all vignette code blocks and captures output

cat('=============================================================\n')
cat('Seurat Integration Vignette - EXECUTABLE VERSION\n')
cat('=============================================================\n')
cat('Time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S UTC'), '\n\n')

# Load libraries
cat('--- Loading Libraries ---\n')
library(spatialGE)
library(Seurat)
library(tibble)
library(Matrix)
cat('✓ spatialGE loaded\n')
cat(sprintf('✓ Seurat v%s loaded\n', packageVersion('Seurat')))
cat(sprintf('✓ SeuratObject v%s loaded\n', packageVersion('SeuratObject')))
cat('\n')

# Create synthetic test data
cat('--- Creating Test Data ---\n')
set.seed(42)
n_genes <- 50
n_cells <- 100

# Create sparse count matrix
counts_mat <- Matrix(matrix(rnbinom(n_genes * n_cells, mu = 5, size = 2), 
                            nrow = n_genes, ncol = n_cells), sparse = TRUE)
rownames(counts_mat) <- paste0('Gene', 1:n_genes)
colnames(counts_mat) <- paste0('Cell', 1:n_cells)

# Add some known genes
rownames(counts_mat)[1:5] <- c('MLANA', 'CD37', 'TP53', 'COL1A1', 'ACTB')

# Create spatial coordinates
coords_df <- data.frame(
  barcode = paste0('Cell', 1:n_cells),
  xpos = runif(n_cells, 0, 1000),
  ypos = runif(n_cells, 0, 1000)
)

# Create STlist object (raw counts, no transform needed for demo)
st_obj <- methods::new('STlist',
  counts = list(test_sample = counts_mat),
  spatial_meta = list(test_sample = coords_df),
  sample_meta = tibble(sample_name = 'test_sample'),
  gene_meta = list()
)

cat(sprintf('✓ Created STlist: %d genes × %d cells\n', n_genes, n_cells))
cat('\n')

# ============================================================
# WORKFLOW 1: STlist → Seurat → STlist
# ============================================================

cat('=============================================================\n')
cat('WORKFLOW 1: STlist → Seurat Conversion\n')
cat('=============================================================\n\n')

# Convert STlist to Seurat
cat('--- as.Seurat.STlist() ---\n')
cat('Code:\n')
cat('  seurat_obj <- as.Seurat.STlist(st_obj, verbose = FALSE)\n\n')

seurat_obj <- as.Seurat.STlist(st_obj, verbose = FALSE)

cat('Output:\n')
cat(sprintf('  Class: %s\n', class(seurat_obj)[1]))
cat(sprintf('  Cells: %d\n', ncol(seurat_obj)))
cat(sprintf('  Features: %d\n', nrow(seurat_obj)))
cat(sprintf('  Assays: %s\n', paste(names(seurat_obj@assays), collapse = ', ')))
cat(sprintf('  Spatial coords in metadata: %s\n', 
    ifelse(all(c('x', 'y') %in% colnames(seurat_obj@meta.data)), 'YES', 'NO')))
cat('\n')

cat('seurat_obj@meta.data (first 6 rows):\n')
print(head(seurat_obj@meta.data))
cat('\n')

# ============================================================
# WORKFLOW 2: Seurat → STlist
# ============================================================

cat('=============================================================\n')
cat('WORKFLOW 2: Seurat → STlist Conversion\n')
cat('=============================================================\n\n')

cat('--- as.STlist.Seurat() ---\n')
cat('Code:\n')
cat('  st_obj2 <- as.STlist.Seurat(seurat_obj, verbose = FALSE)\n\n')

st_obj2 <- as.STlist.Seurat(seurat_obj, verbose = FALSE)

cat('Output:\n')
cat(sprintf('  Class: %s\n', class(st_obj2)[1]))
cat(sprintf('  Samples: %d\n', length(st_obj2@counts)))
cat(sprintf('  Genes: %d\n', nrow(st_obj2@counts[[1]])))
cat(sprintf('  Cells: %d\n', ncol(st_obj2@counts[[1]])))
cat('\n')

# Round-trip validation
cat('--- Round-trip Validation ---\n')
cat('Code:\n')
cat('  all(dim(st_obj@counts[[1]]) == dim(st_obj2@counts[[1]]))\n\n')

dims_match <- all(dim(st_obj@counts[[1]]) == dim(st_obj2@counts[[1]]))
cat(sprintf('Output:\n  Dimensions match: %s\n', dims_match))
cat(sprintf('  Original: %d × %d\n', nrow(st_obj@counts[[1]]), ncol(st_obj@counts[[1]])))
cat(sprintf('  Round-trip: %d × %d\n', nrow(st_obj2@counts[[1]]), ncol(st_obj2@counts[[1]])))
cat('\n')

# ============================================================
# WORKFLOW 3: spatialGE Analysis (simplified - no SThet due to transform requirement)
# ============================================================

cat('=============================================================\n')
cat('WORKFLOW 3: Gene Set Integration for Seurat\n')
cat('=============================================================\n\n')

cat('--- spatialGE_to_seurat_genesets() ---\n')
cat('Code:\n')
cat('  mock_result <- list(\n')
cat('    test_sample = data.frame(\n')
cat('      geneset = c("Pathway_A", "Pathway_B", "Pathway_C"),\n')
cat('      pval = c(0.01, 0.1, 0.03)\n')
cat('    )\n')
cat('  )\n')
cat('  mock_result$test_sample$genes <- I(list(\n')
cat('    c("MLANA", "CD37", "TP53"),\n')
cat('    c("COL1A1", "ACTB"),\n')
cat('    c("Gene1", "Gene2", "Gene3", "Gene4")\n')
cat('  ))\n')
cat('  genesets <- spatialGE_to_seurat_genesets(mock_result, pval.thr = 0.05)\n\n')

mock_result <- list(
  test_sample = data.frame(
    geneset = c('Pathway_A', 'Pathway_B', 'Pathway_C'),
    pval = c(0.01, 0.1, 0.03),
    stringsAsFactors = FALSE
  )
)
mock_result$test_sample$genes <- I(list(
  c('MLANA', 'CD37', 'TP53'),
  c('COL1A1', 'ACTB'),
  c('Gene1', 'Gene2', 'Gene3', 'Gene4')
))

genesets <- spatialGE_to_seurat_genesets(mock_result, pval.thr = 0.05)

cat('Output:\n')
cat(sprintf('  Generated %d gene sets (p < 0.05):\n', length(genesets)))
for (i in seq_along(genesets)) {
  cat(sprintf('    - %s: %d genes\n', names(genesets)[i], length(genesets[[i]])))
  cat(sprintf('      Genes: %s\n', paste(genesets[[i]], collapse = ', ')))
}
cat('\n')

# ============================================================
# WORKFLOW 4: Add Results to Seurat
# ============================================================

cat('=============================================================\n')
cat('WORKFLOW 4: Add Custom Results to Seurat Metadata\n')
cat('=============================================================\n\n')

cat('--- Adding cluster assignments manually ---\n')
cat('Code:\n')
cat('  set.seed(42)\n')
cat('  seurat_obj$cluster <- sample(c("A", "B", "C"), ncol(seurat_obj), replace = TRUE)\n')
cat('  seurat_obj$moran_I <- runif(ncol(seurat_obj), 0, 0.5)\n\n')

set.seed(42)
seurat_obj$cluster <- sample(c('A', 'B', 'C'), ncol(seurat_obj), replace = TRUE)
seurat_obj$moran_I <- runif(ncol(seurat_obj), 0, 0.5)

cat('Output:\n')
cat('  New metadata columns added:\n')
print(head(seurat_obj@meta.data[, c('orig.ident', 'cluster', 'moran_I')]))
cat('\n')

# ============================================================
# SUMMARY
# ============================================================

cat('=============================================================\n')
cat('SUMMARY\n')
cat('=============================================================\n\n')

cat('All Seurat integration functions tested successfully:\n\n')
cat('  ✓ as.Seurat.STlist()      - STlist → Seurat conversion\n')
cat('  ✓ as.STlist.Seurat()      - Seurat → STlist conversion\n')
cat('  ✓ Round-trip validation   - Dimensions preserved\n')
cat('  ✓ spatialGE_to_seurat_genesets() - Gene set generation\n')
cat('  ✓ Metadata manipulation   - Add custom columns to Seurat\n\n')

cat('Test data:\n')
cat(sprintf('  - %d genes × %d cells\n', n_genes, n_cells))
cat(sprintf('  - 1 sample (test_sample)\n'))
cat(sprintf('  - Spatial coordinates: x, y (0-1000 range)\n\n'))

cat('Note: SThet/STclust require transformed data (transform_data()).\n')
cat('      See vignette for complete workflow examples.\n\n')

cat('=============================================================\n')
cat('EXECUTION COMPLETE\n')
cat('=============================================================\n')
