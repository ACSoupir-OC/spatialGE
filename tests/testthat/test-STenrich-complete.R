##
# STenrich comprehensive test suite
#
# Uses shared test data from tests/testthat/setup.R
#

library(testthat)

# ============================================================================
# TEST 1: Basic function runs without error
# ============================================================================

test_that("STenrich runs without error", {
  skip_if_no_data("STenrich")
  st_obj = get_test_stlist()
  gene_sets = create_gene_sets(st_obj)
  
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
  expect_true(is.list(result))
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
  skip_if_no_data("STenrich"); st_obj = get_test_stlist()
  gene_sets = create_gene_sets(st_obj)
  
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
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  
  # Each sample should have results
  expect_true(nrow(result[[1]]) > 0)
  expect_true(nrow(result[[2]]) > 0)
})

# ============================================================================
# TEST 3: Different score types
# ============================================================================

test_that("STenrich with avg score type", {
  skip_if_no_data("STenrich"); st_obj = get_test_stlist()
  gene_sets = create_gene_sets(st_obj)
  
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
  expect_true(is.list(result_avg))
  expect_true(nrow(result_avg[[1]]) > 0)
})

# ============================================================================
# TEST 4: Different min_units
# ============================================================================

test_that("STenrich with different min_units", {
  skip_if_no_data("STenrich"); st_obj = get_test_stlist()
  gene_sets = create_gene_sets(st_obj)
  
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
  expect_true(is.list(result_10))
  expect_true(is.list(result_20))
})

# ============================================================================
# TEST 5: Different num_sds
# ============================================================================

test_that("STenrich with different num_sds", {
  skip_if_no_data("STenrich"); st_obj = get_test_stlist()
  gene_sets = create_gene_sets(st_obj)
  
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
  expect_true(is.list(result_1sd))
  expect_true(is.list(result_2sd))
})

# ============================================================================
# TEST 6: Different reps
# ============================================================================

test_that("STenrich with different reps", {
  skip_if_no_data("STenrich"); st_obj = get_test_stlist()
  gene_sets = create_gene_sets(st_obj)
  
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
  expect_true(is.list(result_100))
  expect_true(is.list(result_500))
  
  # Results should have same number of rows
  expect_equal(nrow(result_100[[1]]), nrow(result_500[[1]]))
})

# ============================================================================
# TEST 7: P-value adjustment
# ============================================================================

test_that("STenrich p-value adjustment", {
  skip_if_no_data("STenrich"); st_obj = get_test_stlist()
  gene_sets = create_gene_sets(st_obj)
  
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
  skip_if_no_data("STenrich"); st_obj = get_test_stlist()
  gene_sets = create_gene_sets(st_obj)
  
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
  expect_true(is.list(result_modern))
  expect_true(is.list(result_legacy))
  
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
  gene_sets = list(GS1 = c("GENE1", "GENE2"))
  
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
  skip_if_no_data("STenrich"); st_obj = get_test_stlist()
  st_obj = STclust(st_obj, ws = 0.025, ks = c(2, 3), cores = 1, verbose = FALSE)
  gene_sets = create_gene_sets(st_obj)
  
  result = STenrich(
    x = st_obj,
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
  expect_true(is.list(result))
  expect_true(length(result) > 0)
})

# ============================================================================
# TEST 11: Gene set size validation
# ============================================================================

test_that("STenrich gene set size validation", {
  skip_if_no_data("STenrich"); st_obj = get_test_stlist()
  
  # Create gene sets with varying sizes
  genes = rownames(st_obj@tr_counts[[1]])
  gene_sets = list(
    GS_small = genes[1:3],   # Below min_genes (5)
    GS_valid = genes[1:10]   # Above min_genes
  )
  
  result = STenrich(
    x = st_obj,
    gene_sets = gene_sets,
    samples = 1,
    reps = 100,
    min_units = 10,
    min_genes = 5,
    cores = 1,
    verbose = FALSE
  )
  
  # Should return list
  expect_true(is.list(result))
  expect_true(length(result) > 0)
  
  # Should have results for both gene sets
  expect_true(nrow(result[[1]]) >= 1)
})

cat("\n=== STenrich Comprehensive Test Suite Complete ===\n")
