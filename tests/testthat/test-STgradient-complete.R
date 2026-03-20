##
# STgradient comprehensive test suite
#
# Uses shared test data from tests/testthat/setup.R
#

library(testthat)

# ============================================================================
# TEST 1: Basic function runs without error
# ============================================================================

test_that("STgradient runs without error", {
  skip_if_no_data("STgradient")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  result = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 20,
    annot = cluster_col,
    ref = "1",
    distsumm = "min",
    cores = 1,
    verbose = FALSE
  )
  
  # Should return a list
  expect_true(is.list(result))
  expect_true(length(result) > 0)
  expect_s3_class(result[[1]], "data.frame")
  expect_true(nrow(result[[1]]) > 0)
})

# ============================================================================
# TEST 2: Different distsumm options
# ============================================================================

test_that("STgradient with different distsumm", {
  skip_if_no_data("STgradient")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  # Test with mean
  result_mean = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 10,
    annot = cluster_col,
    ref = "1",
    distsumm = "mean",
    cores = 1,
    verbose = FALSE
  )
  
  # Test with min
  result_min = STgradient(
    x = st_obj,
    samples = 1,
    topgenes = 10,
    annot = cluster_col,
    ref = "1",
    distsumm = "min",
    cores = 1,
    verbose = FALSE
  )
  
  expect_true(length(result_mean) > 0)
  expect_true(length(result_min) > 0)
})

# ============================================================================
# TEST 3: Invalid annot message
# ============================================================================

test_that("STgradient - Invalid annot message", {
  skip_if_no_data("STgradient")
  st_obj = get_test_stlist()
  
  # Both modern and legacy use message() not warning() for invalid annot
  expect_message(
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 10,
      annot = "INVALID_COLUMN",
      verbose = FALSE
    ),
    "annotation.*not present"
  )
})

# ============================================================================
# TEST 4: Invalid sample error
# ============================================================================

test_that("STgradient - Invalid sample error", {
  skip_if_no_data("STgradient")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  expect_error(
    STgradient(
      x = st_obj,
      samples = 999,
      topgenes = 10,
      annot = cluster_col,
      verbose = FALSE
    )
  )
})

# ============================================================================
# TEST 5: Multiple samples
# ============================================================================

test_that("STgradient - Multiple samples", {
  skip_if_no_data("STgradient")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  result = STgradient(
    x = st_obj,
    samples = c(1, 2),
    topgenes = 5,
    annot = cluster_col,
    ref = "1",
    distsumm = "min",
    cores = 1,
    verbose = FALSE
  )
  
  expect_true(length(result) >= 1)
})

cat('\n========================================\n')
cat('All STgradient tests complete!\n')
cat('========================================\n')
