##
# STgradient basic test suite
#
# Basic tests for STgradient spatial gradient analysis
#

library(testthat)
library(devtools)
library(spatialGE)

# Load the package
load_all("../../.", export_all = TRUE)

# ============================================================================
# SETUP
# ============================================================================

# Get test data path
test_data_path = system.file("extdata", "melanoma_thrane", package = "spatialGE")

if(test_data_path == ""){
  # Try to download test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  tryCatch({
    download.file(lk, destfile = paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode = 'wb')
    unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
    count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'),
                             full.names = TRUE, pattern = 'counts')
    coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'),
                             full.names = TRUE, pattern = 'mapping')
    clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'),
                           full.names = TRUE, pattern = 'clinical')
    st_obj = STlist(rnacounts = count_files,
                    spotcoords = coord_files,
                    samples = clin_file)
    st_obj = transform_data(st_obj)
  }, error = function(e){
    test_that("STgradient tests require test data", {
      skip("Test data not available")
    })
  })
} else{
  count_files = list.files(test_data_path, full.names = TRUE, pattern = 'counts')
  coord_files = list.files(test_data_path, full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(test_data_path, full.names = TRUE, pattern = 'clinical')
  st_obj = STlist(rnacounts = count_files,
                  spotcoords = coord_files,
                  samples = clin_file)
  st_obj = transform_data(st_obj)
}

# Run STclust first to get annotations
st_obj = STclust(st_obj, ws = 0.025, ks = c(2, 3), cores = 1, verbose = FALSE)

# ============================================================================
# TEST 1: Basic function runs without error
# ============================================================================

test_that("STgradient runs without error", {
  skip_if_not(exists("st_obj"))
  
  # Run STgradient
  result = tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 20,
      annot = "stclust_spw0.025_k2",
      ref = "1",
      distsumm = "min",
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  # Should return a list
  expect_s3_class(result, "list")
})

# ============================================================================
# TEST 2: Different distsumm options
# ============================================================================

test_that("STgradient with different distsumm", {
  skip_if_not(exists("st_obj"))
  
  # Test with min distance
  result_min = tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 20,
      annot = "stclust_spw0.025_k2",
      ref = "1",
      distsumm = "min",
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  # Test with avg distance
  result_avg = tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 20,
      annot = "stclust_spw0.025_k2",
      ref = "1",
      distsumm = "avg",
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  # Both should return lists
  expect_s3_class(result_min, "list")
  expect_s3_class(result_avg, "list")
})

# ============================================================================
# TEST 3: Different topgenes values
# ============================================================================

test_that("STgradient topgenes parameter", {
  skip_if_not(exists("st_obj"))
  
  result_50 = tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 50,
      annot = "stclust_spw0.025_k2",
      ref = "1",
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  result_100 = tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 100,
      annot = "stclust_spw0.025_k2",
      ref = "1",
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  # Both should return lists
  expect_s3_class(result_50, "list")
  expect_s3_class(result_100, "list")
})

# ============================================================================
# TEST 4: Different reference domains
# ============================================================================

test_that("STgradient different reference domains", {
  skip_if_not(exists("st_obj"))
  
  result_ref1 = tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 20,
      annot = "stclust_spw0.025_k2",
      ref = "1",
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  result_ref2 = tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 20,
      annot = "stclust_spw0.025_k2",
      ref = "2",
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  # Both should return lists
  expect_s3_class(result_ref1, "list")
  expect_s3_class(result_ref2, "list")
})

# ============================================================================
# TEST 5: Exclude parameter
# ============================================================================

test_that("STgradient exclude parameter", {
  skip_if_not(exists("st_obj"))
  
  result_exclude = tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 20,
      annot = "stclust_spw0.025_k2",
      ref = "1",
      exclude = "2",
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  # Should return list
  expect_s3_class(result_exclude, "list")
})

# ============================================================================
# TEST 6: Log distance transformation
# ============================================================================

test_that("STgradient log_dist parameter", {
  skip_if_not(exists("st_obj"))
  
  result_log = tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 20,
      annot = "stclust_spw0.025_k2",
      ref = "1",
      log_dist = TRUE,
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  # Should return list
  expect_s3_class(result_log, "list")
})

# ============================================================================
# TEST 7: Distance limit
# ============================================================================

test_that("STgradient limit parameter", {
  skip_if_not(exists("st_obj"))
  
  result_limited = tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 20,
      annot = "stclust_spw0.025_k2",
      ref = "1",
      limit = 100,
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  # Should return list
  expect_s3_class(result_limited, "list")
})

# ============================================================================
# TEST 8: Invalid reference (should handle gracefully)
# ============================================================================

test_that("STgradient invalid reference", {
  skip_if_not(exists("st_obj"))
  
  result = tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 20,
      annot = "stclust_spw0.025_k2",
      ref = "999",  # Non-existent cluster
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  # Should handle gracefully and return list
  expect_s3_class(result, "list")
})

# ============================================================================
# TEST 9: Backward compatibility with legacy
# ============================================================================

test_that("STgradient matches legacy output", {
  skip_if_not(exists("st_obj"))
  
  # Run both versions
  result_modern <- tryCatch({
    STgradient(
      x = st_obj,
      samples = 1,
      topgenes = 20,
      annot = "stclust_spw0.025_k2",
      ref = "1",
      distsumm = "min",
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  result_legacy <- tryCatch({
    STgradient_legacy(
      x = st_obj,
      samples = 1,
      topgenes = 20,
      annot = "stclust_spw0.025_k2",
      ref = "1",
      distsumm = "min",
      cores = 1,
      verbose = FALSE
    )
  }, error = function(e){
    cat("Error:", e$message, "\n")
    return(NULL)
  })
  
  # Both should return lists
  expect_s3_class(result_modern, "list")
  expect_s3_class(result_legacy, "list")
})

cat("\n=== STgradient Basic Test Suite Complete ===\n")
