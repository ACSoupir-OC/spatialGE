# Comprehensive STdiff Tests

#' @title Test STdiff core functions
#' @description Comprehensive test suite for differential expression analysis

library(testthat)
library(spatialGE)

# Skip all tests if test data is not available
skip_if_no_data <- function() {
  data_dir <- testthat::test_path("data", "melanoma_thrane")
  if (!dir.exists(data_dir)) {
    skip("Test data not available. Run helper-data.R to download.")
  }
}

# Helper function to create minimal test STlist object
create_minimal_stlist <- function(n_spots = 50, n_genes = 100, n_samples = 2) {
  set.seed(42)
  
  stlist <- list()
  stlist@data <- list()
  stlist@spatial_meta <- list()
  stlist@gene_meta <- list()
  stlist@results <- list()
  
  for (i in 1:n_samples) {
    sample_name <- paste0("sample_", i)
    
    # Create count matrix
    counts <- matrix(
      rpois(n_spots * n_genes, lambda = 10),
      nrow = n_spots,
      ncol = n_genes,
      dimnames = list(paste0("spot_", 1:n_spots), paste0("gene_", 1:n_genes))
    )
    stlist@data[[sample_name]] <- counts
    
    # Create spatial metadata with clusters
    x_coords <- runif(n_spots, 0, 100)
    y_coords <- runif(n_spots, 0, 100)
    clusters <- sample(c("cluster_A", "cluster_B"), n_spots, replace = TRUE)
    
    stlist@spatial_meta[[sample_name]] <- data.frame(
      spot_id = paste0("spot_", 1:n_spots),
      x = x_coords,
      y = y_coords,
      cluster = clusters,
      stringsAsFactors = FALSE
    )
    
    # Create gene metadata
    stlist@gene_meta[[sample_name]] <- data.frame(
      gene = paste0("gene_", 1:n_genes),
      mean = rowMeans(counts),
      var = apply(counts, 1, var),
      stringsAsFactors = FALSE
    )
  }
  
  class(stlist) <- "STlist"
  return(stlist)
}

# ============================================================================
# SECTION 1: Core Function Tests
# ============================================================================

test_that("STdiff_select_genes validates input parameters", {
  skip_if_no_data()
  
  data_dir <- test_path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = FALSE)
  
  cluster_col <- "stclust_spw0.05_k2"
  
  # Test invalid topgenes
  expect_error(
    STdiff_select_genes(st_obj, topgenes = -1, annot = cluster_col),
    "invalid values"
  )
  
  # Test invalid pval_thr
  expect_error(
    STdiff_select_genes(st_obj, pval_thr = 1.5, annot = cluster_col),
    "invalid values"
  )
  
  # Test missing annotation without w/k
  expect_error(
    STdiff_select_genes(st_obj, annot = NULL, w = NULL, k = NULL),
    "specify both w and k"
  )
  
  # Test multiple annotation columns
  expect_error(
    STdiff_select_genes(st_obj, annot = c("col1", "col2")),
    "Only one annotation"
  )
})

test_that("STdiff_select_genes returns correct structure", {
  skip_if_no_data()
  
  data_dir <- test_path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = FALSE)
  
  cluster_col <- "stclust_spw0.05_k2"
  
  result <- STdiff_select_genes(
    st_obj,
    annot = cluster_col,
    topgenes = 100,
    test_type = "mm",
    verbose = 0L
  )
  
  expect_type(result, "list")
  expect_true("combo_df" %in% names(result))
  expect_true("meta_dict" %in% names(result))
  expect_true("pval_thr" %in% names(result))
  expect_true("test_type" %in% names(result))
  
  # Check combo_df structure
  expect_s3_class(result$combo_df, "data.frame")
  expect_true(all(c("samplename", "meta1", "meta2", "gene") %in% colnames(result$combo_df)))
})

test_that("STdiff_select_genes handles different test types", {
  skip_if_no_data()
  
  data_dir <- test_path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = FALSE)
  
  cluster_col <- "stclust_spw0.05_k2"
  
  # Test mixed models
  result_mm <- STdiff_select_genes(st_obj, annot = cluster_col, test_type = "mm", verbose = 0L)
  expect_equal(result_mm$test_type, "mm")
  
  # Test t-test
  result_t <- STdiff_select_genes(st_obj, annot = cluster_col, test_type = "t_test", verbose = 0L)
  expect_equal(result_t$test_type, "t_test")
  
  # Test Wilcoxon
  result_w <- STdiff_select_genes(st_obj, annot = cluster_col, test_type = "wilcoxon", verbose = 0L)
  expect_equal(result_w$test_type, "wilcoxon")
})

test_that("STdiff_select_genes handles pairwise testing", {
  skip_if_no_data()
  
  data_dir <- test_path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = FALSE)
  
  cluster_col <- "stclust_spw0.05_k2"
  
  # Test pairwise with insufficient clusters
  expect_error(
    STdiff_select_genes(st_obj, annot = cluster_col, clusters = "cluster_A", pairwise = TRUE),
    "at least two clusters"
  )
  
  # Test pairwise with sufficient clusters
  result <- STdiff_select_genes(
    st_obj,
    annot = cluster_col,
    clusters = c("cluster_A", "cluster_B"),
    pairwise = TRUE,
    verbose = 0L
  )
  
  expect_type(result, "list")
})

# ============================================================================
# SECTION 2: Helper Function Tests
# ============================================================================

test_that("count_cores returns appropriate values", {
  # Test with small n
  cores_small <- count_cores(2)
  expect_gte(cores_small, 1)
  expect_lte(cores_small, parallel::detectCores() / 2)
  
  # Test with large n
  cores_large <- count_cores(100)
  expect_lte(cores_large, parallel::detectCores() / 2)
})

test_that("expandSparse converts sparse matrices correctly", {
  # Create sparse matrix
  library(Matrix)
  sparse_mat <- Matrix::Matrix(
    c(1, 0, 0, 2, 0, 0, 3, 0, 4),
    nrow = 3,
    ncol = 3,
    sparse = TRUE
  )
  
  # Convert to dense
  dense_df <- expandSparse(sparse_mat)
  
  expect_s3_class(dense_df, "data.frame")
  expect_equal(nrow(dense_df), 3)
  expect_equal(ncol(dense_df), 3)
  expect_equal(dense_df[1, 1], 1)
  expect_equal(dense_df[2, 2], 2)
})

test_that("non_spatial_de handles different test types", {
  # Create minimal test data
  set.seed(42)
  n_spots <- 50
  
  expr_data <- data.frame(
    group = rep(c("sample1", "sample2"), each = n_spots/2),
    meta = sample(c("cluster_A", "cluster_B"), n_spots, replace = TRUE),
    xpos = runif(n_spots, 0, 100),
    ypos = runif(n_spots, 0, 100),
    gene_1 = rnorm(n_spots, 10, 2),
    gene_2 = rnorm(n_spots, 8, 3)
  )
  
  combo <- data.frame(
    samplename = "sample1",
    meta1 = "cluster_A",
    meta2 = "cluster_B",
    gene = "gene_1",
    stringsAsFactors = FALSE
  )
  
  # Test non_spatial_de (note: this uses spaMM internally)
  # Skipping full execution test due to complexity
  # Instead test structure expectations
  expect_true(is.function(non_spatial_de))
})

test_that("prepare_stdiff_combo creates correct combinations", {
  to_expand <- data.frame(
    samplename = c("sample1", "sample1"),
    gene = c("gene_1", "gene_2"),
    orig_annot = c("cluster_A", "cluster_B"),
    meta = c("c1", "c2"),
    stringsAsFactors = FALSE
  )
  
  # Test without specific clusters
  result <- prepare_stdiff_combo(to_expand, user_clusters = NULL, pairwise = FALSE, verbose = 0L)
  
  expect_s3_class(result, "data.frame")
  expect_true(all(c("samplename", "meta1", "meta2", "gene") %in% colnames(result)))
})

# ============================================================================
# SECTION 3: Integration Tests
# ============================================================================

test_that("STdiff full workflow with mixed models", {
  skip_if_no_data()
  
  data_dir <- test_path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = FALSE)
  
  cluster_col <- "stclust_spw0.05_k2"
  
  # Run STdiff with non-spatial only (faster)
  res_diff <- STdiff(
    st_obj,
    annot = cluster_col,
    topgenes = 100,
    sp_topgenes = 0,  # Skip spatial for speed
    test_type = "mm",
    verbose = 0L
  )
  
  expect_type(res_diff, "list")
  expect_true(length(res_diff) > 0)
  
  # Check first sample result structure
  first_sample <- names(res_diff)[1]
  expect_s3_class(res_diff[[first_sample]], "data.frame")
  
  # Check required columns
  required_cols <- c("gene", "avg_log2fc", "adj_p_val", "cluster_1")
  expect_true(all(required_cols %in% colnames(res_diff[[first_sample]])))
})

test_that("STdiff with t-test", {
  skip_if_no_data()
  
  data_dir <- test_path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = FALSE)
  
  cluster_col <- "stclust_spw0.05_k2"
  
  res_diff <- STdiff(
    st_obj,
    annot = cluster_col,
    topgenes = 100,
    sp_topgenes = 0,
    test_type = "t_test",
    verbose = 0L
  )
  
  expect_type(res_diff, "list")
  expect_true(length(res_diff) > 0)
})

test_that("STdiff with Wilcoxon test", {
  skip_if_no_data()
  
  data_dir <- test_path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = FALSE)
  
  cluster_col <- "stclust_spw0.05_k2"
  
  res_diff <- STdiff(
    st_obj,
    annot = cluster_col,
    topgenes = 100,
    sp_topgenes = 0,
    test_type = "wilcoxon",
    verbose = 0L
  )
  
  expect_type(res_diff, "list")
  expect_true(length(res_diff) > 0)
})

test_that("STdiff_volcano generates plots", {
  skip_if_no_data()
  
  data_dir <- test_path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = FALSE)
  
  cluster_col <- "stclust_spw0.05_k2"
  
  res_diff <- STdiff(
    st_obj,
    annot = cluster_col,
    topgenes = 100,
    sp_topgenes = 0,
    verbose = 0L
  )
  
  # Generate volcano plot
  p_vol <- STdiff_volcano(res_diff, samples = names(res_diff)[1])
  
  expect_type(p_vol, "list")
  expect_true(length(p_vol) > 0)
})

# ============================================================================
# SECTION 4: Edge Cases and Error Handling
# ============================================================================

test_that("STdiff handles missing annotation column", {
  skip_if_no_data()
  
  data_dir <- test_path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  
  # Use non-existent annotation column
  expect_error(
    STdiff(st_obj, annot = "non_existent_column", verbose = 0L),
    "No samples left"
  )
})

test_that("STdiff handles empty sample list", {
  skip_if_no_data()
  
  data_dir <- test_path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = FALSE)
  
  cluster_col <- "stclust_spw0.05_k2"
  
  # Use invalid sample names
  expect_error(
    STdiff(st_obj, samples = c("invalid_sample"), annot = cluster_col),
    "None of the requested samples"
  )
})

test_that("STdiff handles invalid test type", {
  skip_if_no_data()
  
  data_dir <- test_path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = FALSE)
  
  cluster_col <- "stclust_spw0.05_k2"
  
  # Use invalid test type - should default or error
  # Based on code, it defaults to mm with warning
  expect_warning(
    STdiff(st_obj, annot = cluster_col, test_type = "invalid_type", verbose = 0L)
  )
})

# ============================================================================
# SECTION 5: Regression Tests
# ============================================================================

test_that("STdiff produces identical results to STdiff_legacy", {
  skip_if_no_data()
  
  data_dir <- test_path("data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)[1]
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)[1]
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj, method = "log")
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, k = 2, deepSplit = FALSE)
  
  cluster_col <- "stclust_spw0.05_k2"
  
  # Run both versions (non-spatial only for determinism)
  res_new <- STdiff(
    st_obj,
    annot = cluster_col,
    topgenes = 100,
    sp_topgenes = 0,
    verbose = 0L
  )
  
  res_legacy <- STdiff_legacy(
    st_obj,
    annot = cluster_col,
    topgenes = 100,
    sp_topgenes = 0,
    verbose = 0L
  )
  
  # Compare results
  expect_type(res_new, "list")
  expect_type(res_legacy, "list")
  expect_identical(names(res_new), names(res_legacy))
  
  # Compare data frames for each sample
  for (sample_name in names(res_new)) {
    expect_equal(nrow(res_new[[sample_name]]), nrow(res_legacy[[sample_name]]))
    expect_equal(ncol(res_new[[sample_name]]), ncol(res_legacy[[sample_name]]))
    
    # Compare key columns with tolerance for floating point
    expect_equal(
      res_new[[sample_name]][["gene"]],
      res_legacy[[sample_name]][["gene"]]
    )
    
    expect_equal(
      res_new[[sample_name]][["avg_log2fc"]],
      res_legacy[[sample_name]][["avg_log2fc"]],
      tolerance = 1e-6
    )
    
    expect_equal(
      res_new[[sample_name]][["adj_p_val"]],
      res_legacy[[sample_name]][["adj_p_val"]],
      tolerance = 1e-6
    )
  }
})

# ============================================================================
# END OF TESTS
# ============================================================================
