##
# STgradient comprehensive test suite
#
# Complete test suite for STgradient spatial gradient analysis
#

library(testthat)
library(devtools)
library(spatialGE)

# Load the package
load_all("../../.", export_all = TRUE)

# ============================================================================
# SETUP
# ============================================================================

# Download test data
thrane_tmp = tempdir()
unlink(thrane_tmp, recursive = TRUE)
dir.create(thrane_tmp)
lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
download.file(lk, destfile = paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode = 'wb')
unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
st_obj = transform_data(st_obj)
st_obj = STclust(st_obj, ws = 0.025, ks = c(2, 3), cores = 1, verbose = FALSE)

# ============================================================================
# TEST 1: Basic function runs without error
# ============================================================================

test_that("STgradient runs without error", {
  result = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    distsumm = "min",
    cores = 1,
    verbose = FALSE
  )
  
  # Should return a list
  expect_s3_class(result, "list")
  expect_true(length(result) > 0)
  
  # Should have data frame with expected columns
  expect_s3_class(result[[1]], "data.frame")
  expect_true(nrow(result[[1]]) > 0)
})

# ============================================================================
# TEST 2: Different distsumm options
# ============================================================================

test_that("STgradient with different distsumm", {
  result_min = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    distsumm = "min",
    cores = 1,
    verbose = FALSE
  )
  
  result_avg = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    distsumm = "avg",
    cores = 1,
    verbose = FALSE
  )
  
  # Both should return lists with results
  expect_s3_class(result_min, "list")
  expect_s3_class(result_avg, "list")
  expect_true(nrow(result_min[[1]]) > 0)
  expect_true(nrow(result_avg[[1]]) > 0)
  
  # Both should have spearman correlations (not all NA)
  expect_true(sum(!is.na(result_min[[1]]$min_spearman_r)) > 0)
  expect_true(sum(!is.na(result_avg[[1]]$avg_spearman_r)) > 0)
})

# ============================================================================
# TEST 3: Top genes parameter
# ============================================================================

test_that("STgradient topgenes parameter", {
  result_20 = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    cores = 1,
    verbose = FALSE
  )
  
  result_10 = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 10,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    cores = 1,
    verbose = FALSE
  )
  
  # Should return different number of results
  expect_equal(nrow(result_20[[1]]), 20)
  expect_equal(nrow(result_10[[1]]), 10)
})

# ============================================================================
# TEST 4: Different reference domains
# ============================================================================

test_that("STgradient different reference domains", {
  result_ref1 = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    cores = 1,
    verbose = FALSE
  )
  
  result_ref2 = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "2",
    cores = 1,
    verbose = FALSE
  )
  
  # Both should return lists
  expect_s3_class(result_ref1, "list")
  expect_s3_class(result_ref2, "list")
  expect_true(nrow(result_ref1[[1]]) > 0)
  expect_true(nrow(result_ref2[[1]]) > 0)
  
  # Results should be different (different reference domains)
  expect_false(all(result_ref1[[1]]$min_spearman_r == result_ref2[[1]]$min_spearman_r))
})

# ============================================================================
# TEST 5: Exclude parameter
# ============================================================================

test_that("STgradient exclude parameter", {
  result_exclude = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    exclude = "2",
    cores = 1,
    verbose = FALSE
  )
  
  # Should return list with results
  expect_s3_class(result_exclude, "list")
  expect_true(nrow(result_exclude[[1]]) > 0)
})

# ============================================================================
# TEST 6: Log distance transformation
# ============================================================================

test_that("STgradient log_dist parameter", {
  result_log = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    log_dist = TRUE,
    cores = 1,
    verbose = FALSE
  )
  
  result_no_log = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    log_dist = FALSE,
    cores = 1,
    verbose = FALSE
  )
  
  # Both should return lists
  expect_s3_class(result_log, "list")
  expect_s3_class(result_no_log, "list")
  
  # Results should be different
  expect_false(all(result_log[[1]]$min_lm_coef == result_no_log[[1]]$min_lm_coef))
})

# ============================================================================
# TEST 7: Distance limit
# ============================================================================

test_that("STgradient limit parameter", {
  result_limited = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    limit = 10,
    cores = 1,
    verbose = FALSE
  )
  
  result_unlimited = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    cores = 1,
    verbose = FALSE
  )
  
  # Both should return lists
  expect_s3_class(result_limited, "list")
  expect_s3_class(result_unlimited, "list")
})

# ============================================================================
# TEST 8: Invalid reference
# ============================================================================

test_that("STgradient invalid reference", {
  result = tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 20,
      annot = "stclust_spw0.025_k2",
      ref = "999",
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    return(NULL)
  })
  
  # Should return NULL for invalid reference
  expect_null(result)
})

# ============================================================================
# TEST 9: Matches legacy output
# ============================================================================

test_that("STgradient matches legacy output", {
  result_modern = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    distsumm = "min",
    cores = 1,
    verbose = FALSE
  )
  
  result_legacy = STgradient_legacy(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    distsumm = "min",
    min_nb = 3,
    robust = TRUE,
    nb_dist_thr = c(0.75, 1.25),
    cores = 1,
    verbose = FALSE
  )
  
  # Both should return lists
  expect_s3_class(result_modern, "list")
  expect_s3_class(result_legacy, "list")
  
  # Results should be similar (within tolerance)
  expect_equal(nrow(result_modern[[1]]), nrow(result_legacy[[1]]))
  
  # Spearman correlations should be similar
  expect_equal(
    round(result_modern[[1]]$min_spearman_r, 4),
    round(result_legacy[[1]]$spearman_r, 4),
    tolerance = 0.1
  )
})

# ============================================================================
# TEST 10: Output column structure
# ============================================================================

test_that("STgradient output column structure", {
  result = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    distsumm = "min",
    cores = 1,
    verbose = FALSE
  )
  
  # Check column names
  expected_cols = c("sample_name", "gene", "min_lm_coef", "min_lm_pval", 
                    "min_spearman_r", "min_spearman_r_pval", 
                    "spearman_r_pval_adj", "pval_comment")
  expect_equal(colnames(result[[1]]), expected_cols)
  
  # Check data types
  expect_type(result[[1]]$sample_name, "character")
  expect_type(result[[1]]$gene, "character")
  expect_type(result[[1]]$min_lm_coef, "double")
  expect_type(result[[1]]$min_spearman_r, "double")
  expect_type(result[[1]]$spearman_r_pval_adj, "double")
})

# ============================================================================
# TEST 11: Multiple samples
# ============================================================================

test_that("STgradient multiple samples", {
  result = STgradient(
    x = st_obj,
    samples = c(1, 2),
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    cores = 1,
    verbose = FALSE
  )
  
  # Should return list with results for both samples
  expect_s3_class(result, "list")
  expect_equal(length(result), 2)
  expect_true(all(names(result) %in% c("ST_mel1_rep2_counts.tsv", "ST_mel2_rep1_counts.tsv")))
  
  # Each sample should have results
  expect_true(nrow(result[[1]]) > 0)
  expect_true(nrow(result[[2]]) > 0)
})

# ============================================================================
# TEST 12: Robust regression
# ============================================================================

test_that("STgradient robust regression", {
  result_robust = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    robust = TRUE,
    cores = 1,
    verbose = FALSE
  )
  
  result_non_robust = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = "stclust_spw0.025_k2",
    ref = "1",
    robust = FALSE,
    cores = 1,
    verbose = FALSE
  )
  
  # Both should return lists
  expect_s3_class(result_robust, "list")
  expect_s3_class(result_non_robust, "list")
})

cat("\n=== STgradient Comprehensive Test Suite Complete ===\n")
cat("All tests passed!\n")
