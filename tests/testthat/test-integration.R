##
# Integration Tests - End-to-End Workflows
#
# Tests that all components work together. Uses pre-computed results
# where possible to avoid long runtimes. Intermediate outputs saved
# to archive/ (not tracked by git).
#

library(testthat)

# ============================================================================
# Setup
# ============================================================================

# Create archive directory for intermediate outputs
archive_dir = file.path(getwd(), "archive")
if (!dir.exists(archive_dir)) {
  dir.create(archive_dir, recursive = TRUE)
}

# Note: setup.R is automatically loaded by testthat before tests run

# ============================================================================
# TEST 1: Basic Spatial Workflow - Component Integration
# ============================================================================
# Verifies: STlist → STclust → STdiff → STenrich → STplot can be called in sequence

test_that("Basic spatial workflow: All functions callable", {
  skip_if_no_data("integration")
  
  cat("\n=== Integration Test 1: Basic Spatial Workflow ===\n")
  
  # Step 1: Load test data
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  cat("Step 1: STlist loaded -", length(st_obj@tr_counts), "samples\n")
  
  expect_s4_class(st_obj, "STlist")
  expect_true(length(st_obj@tr_counts) > 0)
  
  # Step 2: Verify STclust results exist (from setup)
  expect_true(cluster_col %in% colnames(st_obj@spatial_meta[[1]]))
  cat("Step 2: STclust verified - column:", cluster_col, "\n")
  
  # Step 3: STdiff - Verify function is available and callable
  # (Full STdiff is tested in test-STdiff-complete.R)
  cat("Step 3: STdiff function available\n")
  expect_true(exists("STdiff"))
  expect_type(STdiff, "closure")
  
  # Step 4: STenrich - Run with minimal parameters
  cat("Step 4: Running STenrich (100 permutations)...\n")
  gene_sets = create_gene_sets(st_obj)
  stenrich_result = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 100,
    cores = 1,
    verbose = FALSE
  )
  
  expect_true(is.list(stenrich_result))
  expect_true(length(stenrich_result) > 0)
  expect_true(nrow(stenrich_result[[1]]) > 0)
  cat("  STenrich complete -", nrow(stenrich_result[[1]]), "gene sets\n")
  
  # Save STenrich output to archive
  stenrich_file = file.path(archive_dir, "integration_stenrich_result.rds")
  saveRDS(stenrich_result, stenrich_file)
  cat("  Saved to:", stenrich_file, "\n")
  
  # Step 5: STplot - Verify plotting functions work
  cat("Step 5: Testing STplot functions...\n")
  
  # Plot spatial expression
  p1 = plot_spatial_expression(
    x = st_obj,
    sample = 1,
    genes = rownames(st_obj@tr_counts[[1]])[1]
  )
  expect_s3_class(p1, "gg")
  cat("  plot_spatial_expression OK\n")
  
  # Plot gene set
  p2 = plot_spatial_geneset(
    x = st_obj,
    sample = 1,
    geneset = gene_sets[[1]]
  )
  expect_s3_class(p2, "gg")
  cat("  plot_spatial_geneset OK\n")
  
  # Save plots to archive
  suppressMessages({
    ggsave(file.path(archive_dir, "integration_plot_spatial_expression.pdf"), p1, dpi = 150)
    ggsave(file.path(archive_dir, "integration_plot_geneset.pdf"), p2, dpi = 150)
  })
  cat("  Plots saved to archive/\n")
  
  cat("\n=== Integration Test 1: PASSED ===\n\n")
})

# ============================================================================
# TEST 2: Gradient Analysis Workflow
# ============================================================================
# Verifies: STlist → STgradient → visualization

test_that("Gradient analysis workflow: STgradient callable", {
  skip_if_no_data("integration")
  
  cat("\n=== Integration Test 2: Gradient Analysis Workflow ===\n")
  
  # Step 1: Load test data
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  cat("Step 1: STlist loaded -", length(st_obj@tr_counts), "samples\n")
  
  # Step 2: STgradient - Run with minimal genes for speed
  cat("Step 2: Running STgradient (5 genes)...\n")
  stgradient_result = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 5,
    annot = cluster_col,
    ref = "1",
    distsumm = "min",
    cores = 1,
    verbose = FALSE
  )
  
  expect_true(is.list(stgradient_result))
  expect_true(length(stgradient_result) > 0)
  expect_true(nrow(stgradient_result[[1]]) > 0)
  expect_true("gene" %in% colnames(stgradient_result[[1]]))
  cat("  STgradient complete -", nrow(stgradient_result[[1]]), "genes\n")
  
  # Save STgradient output to archive
  stgradient_file = file.path(archive_dir, "integration_stgradient_result.rds")
  saveRDS(stgradient_result, stgradient_file)
  cat("  Saved to:", stgradient_file, "\n")
  
  cat("\n=== Integration Test 2: PASSED ===\n\n")
})

# ============================================================================
# TEST 3: Heterogeneity Workflow
# ============================================================================
# Verifies: STlist → SThet → results

test_that("Heterogeneity workflow: SThet callable", {
  skip_if_no_data("integration")
  
  cat("\n=== Integration Test 3: Heterogeneity Workflow ===\n")
  
  # Step 1: Load test data
  st_obj = get_test_stlist()
  cat("Step 1: STlist loaded -", length(st_obj@tr_counts), "samples\n")
  
  # Step 2: SThet - Calculate heterogeneity metrics
  cat("Step 2: Running SThet (5 genes)...\n")
  sthet_result = SThet(
    x = st_obj,
    samples = 1,
    genes = rownames(st_obj@tr_counts[[1]])[1:5],
    cores = 1,
    verbose = FALSE
  )
  
  # Verify SThet output (returns modified STlist with results in @gene_meta)
  expect_s4_class(sthet_result, "STlist")
  
  # Check that results are stored in gene_meta
  expect_true("moran_i" %in% colnames(sthet_result@gene_meta[[1]]))
  
  # Extract moran values for tested genes
  test_genes = rownames(st_obj@tr_counts[[1]])[1:5]
  moran_vals = sthet_result@gene_meta[[1]]$moran_i[
    sthet_result@gene_meta[[1]]$gene %in% test_genes
  ]
  
  expect_true(length(moran_vals) > 0)
  cat("  SThet complete -", length(moran_vals), "genes analyzed\n")
  
  # Verify statistical values are in expected range (Moran's I: -1 to 1)
  expect_true(all(moran_vals >= -1 & moran_vals <= 1, na.rm = TRUE))
  cat("  Statistical values in expected range\n")
  
  # Save SThet output to archive
  sthet_file = file.path(archive_dir, "integration_sthet_result.rds")
  saveRDS(sthet_result, sthet_file)
  cat("  Saved to:", sthet_file, "\n")
  
  cat("\n=== Integration Test 3: PASSED ===\n\n")
})

# ============================================================================
# TEST 4: Error Handling Across Pipeline
# ============================================================================

test_that("Error handling: Invalid inputs caught", {
  skip_if_no_data("integration")
  
  cat("\n=== Integration Test 4: Error Handling ===\n")
  
  st_obj = get_test_stlist()
  
  # Test 1: Invalid annotation in STgradient
  cat("Test 1: Invalid annotation in STgradient...\n")
  expect_message(
    STgradient(x = st_obj, samples = 1, annot = "INVALID", verbose = FALSE),
    "annotation.*not present"
  )
  cat("  STgradient error handling OK\n")
  
  # Test 2: Invalid samples in STenrich
  cat("Test 2: Invalid samples in STenrich...\n")
  gene_sets = create_gene_sets(st_obj)
  expect_error(
    STenrich(x = st_obj, gene_sets = gene_sets, samples = 999),
    "sample"
  )
  cat("  STenrich error handling OK\n")
  
  cat("\n=== Integration Test 4: PASSED ===\n\n")
})

# ============================================================================
# TEST 5: Reproducibility
# ============================================================================

test_that("Reproducibility: Same inputs produce same outputs", {
  skip_if_no_data("integration")
  
  cat("\n=== Integration Test 5: Reproducibility ===\n")
  
  st_obj = get_test_stlist()
  test_genes = rownames(st_obj@tr_counts[[1]])[1:3]
  
  # Run SThet twice with same seed
  set.seed(42)
  result1 = SThet(
    x = st_obj,
    samples = 1,
    genes = test_genes,
    cores = 1,
    verbose = FALSE
  )
  
  set.seed(42)
  result2 = SThet(
    x = st_obj,
    samples = 1,
    genes = test_genes,
    cores = 1,
    verbose = FALSE
  )
  
  # Results should be identical (extract from @gene_meta)
  moran1 = result1@gene_meta[[1]]$moran_i[result1@gene_meta[[1]]$gene %in% test_genes]
  moran2 = result2@gene_meta[[1]]$moran_i[result2@gene_meta[[1]]$gene %in% test_genes]
  expect_equal(moran1, moran2)
  cat("  SThet reproducibility OK\n")
  
  # Run STenrich twice with same seed
  gene_sets = create_gene_sets(st_obj)
  set.seed(42)
  enrich1 = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 50,
    cores = 1,
    verbose = FALSE
  )
  
  set.seed(42)
  enrich2 = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 50,
    cores = 1,
    verbose = FALSE
  )
  
  # P-values should be identical (permutation test with same seed)
  expect_equal(enrich1[[1]]$p_value, enrich2[[1]]$p_value)
  cat("  STenrich reproducibility OK\n")
  
  cat("\n=== Integration Test 5: PASSED ===\n\n")
})

# ============================================================================
# Cleanup
# ============================================================================

cat("\n=== Integration Tests Complete ===\n")
cat("Archive directory:", archive_dir, "\n")
cat("Files saved:\n")
if (dir.exists(archive_dir)) {
  files = list.files(archive_dir, full.names = TRUE)
  for (f in files) {
    cat("  -", f, "\n")
  }
}
cat("\n")
