##
# STenrich comprehensive test suite
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

# Define gene sets
genes = rownames(st_obj@tr_counts[[1]])
gene_sets = list(
  GS1 = genes[1:10],
  GS2 = genes[11:20],
  GS3 = genes[21:30]
)

# ============================================================================
# TEST 1: Basic function runs without error
# ============================================================================

test_that("STenrich runs without error", {
  result = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 100,
    min_units = 10,
    cores = 1,
    verbose = FALSE
  )
  
  # Should return a list
  expect_s3_class(result, "list")
  expect_true(length(result) > 0)
  
  # Should have data frame with expected columns
  expect_s3_class(result[[1]], "data.frame")
  expect_true(nrow(result[[1]]) > 0)
  
  # Check column names
  expected_cols = c("sample_name", "gene_set", "size_test", "size_gene_set", 
                    "prop_size_test", "p_value", "adj_p_value")
  expect_equal(colnames(result[[1]]), expected_cols)
})

# ============================================================================
# TEST 2: Multiple samples
# ============================================================================

test_that("STenrich with multiple samples", {
  result = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = c(1, 2),
    reps = 100,
    min_units = 10,
    cores = 1,
    verbose = FALSE
  )
  
  # Should return list with results for both samples
  expect_s3_class(result, "list")
  expect_equal(length(result), 2)
  
  # Each sample should have results
  expect_true(nrow(result[[1]]) > 0)
  expect_true(nrow(result[[2]]) > 0)
})

# ============================================================================
# TEST 3: Different score types
# ============================================================================

test_that("STenrich with avg score type", {
  result_avg = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    score_type = "avg",
    reps = 100,
    min_units = 10,
    cores = 1,
    verbose = FALSE
  )
  
  # Should return list with results
  expect_s3_class(result_avg, "list")
  expect_true(nrow(result_avg[[1]]) > 0)
})

# ============================================================================
# TEST 4: Different min_units
# ============================================================================

test_that("STenrich with different min_units", {
  result_10 = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 100,
    min_units = 10,
    cores = 1,
    verbose = FALSE
  )
  
  result_20 = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 100,
    min_units = 20,
    cores = 1,
    verbose = FALSE
  )
  
  # Both should return lists
  expect_s3_class(result_10, "list")
  expect_s3_class(result_20, "list")
})

# ============================================================================
# TEST 5: Different num_sds
# ============================================================================

test_that("STenrich with different num_sds", {
  result_1sd = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 100,
    num_sds = 1,
    min_units = 10,
    cores = 1,
    verbose = FALSE
  )
  
  result_2sd = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 100,
    num_sds = 2,
    min_units = 10,
    cores = 1,
    verbose = FALSE
  )
  
  # Both should return lists
  expect_s3_class(result_1sd, "list")
  expect_s3_class(result_2sd, "list")
})

# ============================================================================
# TEST 6: Different reps
# ============================================================================

test_that("STenrich with different reps", {
  result_100 = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 100,
    min_units = 10,
    cores = 1,
    verbose = FALSE
  )
  
  result_500 = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 500,
    min_units = 10,
    cores = 1,
    verbose = FALSE
  )
  
  # Both should return lists
  expect_s3_class(result_100, "list")
  expect_s3_class(result_500, "list")
  
  # Results should have same number of rows
  expect_equal(nrow(result_100[[1]]), nrow(result_500[[1]]))
})

# ============================================================================
# TEST 7: P-value adjustment
# ============================================================================

test_that("STenrich p-value adjustment", {
  result = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 100,
    min_units = 10,
    pval_adj_method = "BH",
    cores = 1,
    verbose = FALSE
  )
  
  # Should have adjusted p-values
  expect_true("adj_p_value" %in% colnames(result[[1]]))
  
  # Adjusted p-values should be between 0 and 1
  expect_true(all(result[[1]]$adj_p_value >= 0, na.rm = TRUE))
  expect_true(all(result[[1]]$adj_p_value <= 1, na.rm = TRUE))
})

# ============================================================================
# TEST 8: Matches legacy output
# ============================================================================

test_that("STenrich matches legacy output", {
  result_modern = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 100,
    min_units = 10,
    cores = 1,
    verbose = FALSE
  )
  
  result_legacy = STenrich_legacy(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 100,
    min_units = 10,
    cores = 1,
    verbose = FALSE
  )
  
  # Both should return lists
  expect_s3_class(result_modern, "list")
  expect_s3_class(result_legacy, "list")
  
  # Results should be identical
  expect_equal(nrow(result_modern[[1]]), nrow(result_legacy[[1]]))
  expect_equal(
    round(result_modern[[1]]$p_value, 4),
    round(result_legacy[[1]]$p_value, 4)
  )
})

# ============================================================================
# TEST 9: Invalid inputs
# ============================================================================

test_that("STenrich invalid STlist", {
  expect_error({
    STenrich(
      x = NULL,
      gene_sets = gene_sets,
      samples = 1
    )
  })
})

# ============================================================================
# TEST 10: Domain filtering
# ============================================================================

test_that("STenrich with domain filtering", {
  # First run clustering to get domains
  st_obj_clustered = STclust(st_obj, ws = 0.025, ks = c(2, 3), cores = 1, verbose = FALSE)
  
  result = STenrich(
    x = st_obj_clustered,
    gene_sets = gene_sets,
    samples = 1,
    annot = "stclust_spw0.025_k2",
    domain = "1",
    reps = 100,
    min_units = 10,
    cores = 1,
    verbose = FALSE
  )
  
  # Should return list with results
  expect_s3_class(result, "list")
  expect_true(length(result) > 0)
})

cat("\n=== STenrich Comprehensive Test Suite Complete ===\n")
