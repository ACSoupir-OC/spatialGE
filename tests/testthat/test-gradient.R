test_that("STgradient (new implementation) works", {
library(spatialGE)
  
  # Load data from test data directory
  data_dir <- file.path(getwd(), "data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj)
  
  # Run STclust first (required for annot column)
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, ks = 'dtc', deepSplit = FALSE)
  
  # Identify a reference cluster
  ref_cluster <- "1" # Assuming clusters are 1, 2...
  
  # Run new STgradient - use correct column name from STclust output
  res_gradient <- STgradient(st_obj, samples = 1, annot = "stclust_spw0.05_dsplFalse", 
                              ref = ref_cluster, topgenes = 50, verbose = FALSE)
  
  expect_type(res_gradient, "list")
  expect_true(length(res_gradient) > 0)
  expect_s3_class(res_gradient[[1]], "data.frame")
  expect_true("min_lm_coef" %in% colnames(res_gradient[[1]]))
  expect_true("min_spearman_r_pval_adj" %in% colnames(res_gradient[[1]]))
})


test_that("STgradient_legacy still works for reproducibility", {
library(spatialGE)
  
  # Load data from test data directory
  data_dir <- file.path(getwd(), "data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj)
  
  # Run STclust first (required for annot column)
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, ks = 'dtc', deepSplit = FALSE)
  
  # Identify a reference cluster
  ref_cluster <- "1"
  
  # Run legacy
  res_legacy <- STgradient_legacy(st_obj, samples = 1, annot = "stclust_spw0.05_dsplFalse", 
                                   ref = ref_cluster, topgenes = 50, verbose = FALSE)
  
  expect_type(res_legacy, "list")
  expect_true(length(res_legacy) > 0)
  expect_s3_class(res_legacy[[1]], "data.frame")
  expect_true("min_lm_coef" %in% colnames(res_legacy[[1]]))
  expect_true("min_spearman_r_pval_adj" %in% colnames(res_legacy[[1]]))
})


test_that("STgradient produces identical results to STgradient_legacy", {
library(spatialGE)
  
  # Load data from test data directory
  data_dir <- file.path(getwd(), "data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj)
  
  # Run STclust first (required for annot column)
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, ks = 'dtc', deepSplit = FALSE)
  
  ref_cluster <- "1"
  
  # Run both implementations
  res_new <- STgradient(st_obj, samples = 1, annot = "stclust_spw0.05_dsplFalse", 
                         ref = ref_cluster, topgenes = 50, verbose = FALSE)
  res_legacy <- STgradient_legacy(st_obj, samples = 1, annot = "stclust_spw0.05_dsplFalse", 
                                   ref = ref_cluster, topgenes = 50, verbose = FALSE)
  
  # Compare results
  expect_equal(nrow(res_new[[1]]), nrow(res_legacy[[1]]))
  expect_equal(ncol(res_new[[1]]), ncol(res_legacy[[1]]))
  expect_identical(colnames(res_new[[1]]), colnames(res_legacy[[1]]))
  
  # Compare key columns (with tolerance for floating point)
  expect_equal(res_new[[1]]$min_lm_coef, res_legacy[[1]]$min_lm_coef, 
               tolerance = 1e-10)
  expect_equal(res_new[[1]]$min_lm_pval, res_legacy[[1]]$min_lm_pval, 
               tolerance = 1e-10)
  expect_equal(res_new[[1]]$min_spearman_r, res_legacy[[1]]$min_spearman_r, 
               tolerance = 1e-10)
  expect_equal(res_new[[1]]$min_spearman_r_pval, res_legacy[[1]]$min_spearman_r_pval, 
               tolerance = 1e-10)
  expect_equal(res_new[[1]]$min_spearman_r_pval_adj, res_legacy[[1]]$min_spearman_r_pval_adj, 
               tolerance = 1e-10)
})


test_that("STgradient handles different distance summaries", {
library(spatialGE)
  
  # Load data from test data directory
  data_dir <- file.path(getwd(), "data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj)
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, ks = 'dtc', deepSplit = FALSE)
  
  # Test 'min' distance summary
  res_min <- STgradient(st_obj, samples = 1, annot = "stclust_spw0.05_dsplFalse", 
                         ref = "1", topgenes = 50, distsumm = 'min', verbose = FALSE)
  expect_true("min_lm_coef" %in% colnames(res_min[[1]]))
  
  # Test 'avg' distance summary
  res_avg <- STgradient(st_obj, samples = 1, annot = "stclust_spw0.05_dsplFalse", 
                         ref = "1", topgenes = 50, distsumm = 'avg', verbose = FALSE)
  expect_true("avg_lm_coef" %in% colnames(res_avg[[1]]))
})


test_that("STgradient handles edge cases (no gradients found)", {
library(spatialGE)
  
  # Load data from test data directory
  data_dir <- file.path(getwd(), "data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj)
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, ks = 'dtc', deepSplit = FALSE)
  
  # Run with parameters that might not find gradients
  res <- STgradient(st_obj, samples = 1, annot = "stclust_spw0.05_dsplFalse", 
                     ref = "1", topgenes = 10, min_nb = 100, verbose = FALSE)
  
  # Should return empty list or list with empty data frames
  expect_type(res, "list")
})


test_that("STgradient handles robust=FALSE to avoid convergence warnings", {
library(spatialGE)
  
  # Load data from test data directory
  data_dir <- file.path(getwd(), "data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj)
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, ks = 'dtc', deepSplit = FALSE)
  
  # Run with robust=FALSE to bypass MASS::rlm convergence issues
  # This should NOT produce the 'rlm failed to converge' warning
  res_no_robust <- STgradient(st_obj, samples = 1, annot = "stclust_spw0.05_dsplFalse", 
                               ref = "1", topgenes = 50, robust = FALSE, verbose = FALSE)
  
  expect_type(res_no_robust, "list")
  expect_true(length(res_no_robust) > 0)
  expect_s3_class(res_no_robust[[1]], "data.frame")
})


test_that("STgradient produces results with robust=TRUE (may warn)", {
library(spatialGE)
  
  # Load data from test data directory
  data_dir <- file.path(getwd(), "data", "melanoma_thrane")
  counts <- list.files(data_dir, pattern = "counts", full.names = TRUE)
  coords <- list.files(data_dir, pattern = "mapping", full.names = TRUE)
  clinical <- file.path(data_dir, "thrane_clinical.csv")
  
  st_obj <- STlist(rnacounts = counts, spotcoords = coords, samples = clinical)
  st_obj <- transform_data(st_obj)
  st_obj <- STclust(st_obj, topgenes = 100, ws = 0.05, ks = 'dtc', deepSplit = FALSE)
  
  # Run with robust=TRUE (default) - expect convergence warning for some genes
  # The warning is expected and documented behavior - users should be aware of it
  res_robust <- suppressWarnings(
    STgradient(st_obj, samples = 1, annot = "stclust_spw0.05_dsplFalse", 
              ref = "1", topgenes = 50, robust = TRUE, verbose = FALSE)
  )
  
  expect_type(res_robust, "list")
  expect_true(length(res_robust) > 0)
  expect_s3_class(res_robust[[1]], "data.frame")
})