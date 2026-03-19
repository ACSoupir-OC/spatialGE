##
# STgradient comprehensive test suite
#

library(testthat)
library(devtools)
library(spatialGE)

# Load the package
load_all("../../.", export_all = TRUE)

# ============================================================================
# TEST 1: Basic function runs without error
# ============================================================================

test_that("STgradient runs without error", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb')
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  st_obj = STclust(st_obj, ws = 0.025, ks = c(2, 3), cores = 1, verbose = FALSE)
  
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
  expect_s3_class(result[[1]], "data.frame")
  expect_true(nrow(result[[1]]) > 0)
})

# ============================================================================
# TEST 2: Different distsumm options
# ============================================================================

test_that("STgradient with different distsumm", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb')
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  st_obj = STclust(st_obj, ws = 0.025, ks = c(2, 3), cores = 1, verbose = FALSE)
  
  result_min = STgradient(x = st_obj, samples = 1, topgenes = 20, annot = "stclust_spw0.025_k2", ref = "1", distsumm = "min", cores = 1, verbose = FALSE)
  result_avg = STgradient(x = st_obj, samples = 1, topgenes = 20, annot = "stclust_spw0.025_k2", ref = "1", distsumm = "avg", cores = 1, verbose = FALSE)
  
  expect_s3_class(result_min, "list")
  expect_s3_class(result_avg, "list")
  expect_true(sum(!is.na(result_min[[1]]$min_spearman_r)) > 0)
  expect_true(sum(!is.na(result_avg[[1]]$avg_spearman_r)) > 0)
})

# ============================================================================
# TEST 3: Top genes parameter
# ============================================================================

test_that("STgradient topgenes parameter", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb')
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  st_obj = STclust(st_obj, ws = 0.025, ks = c(2, 3), cores = 1, verbose = FALSE)
  
  result_20 = STgradient(x = st_obj, samples = 1, topgenes = 20, annot = "stclust_spw0.025_k2", ref = "1", cores = 1, verbose = FALSE)
  result_10 = STgradient(x = st_obj, samples = 1, topgenes = 10, annot = "stclust_spw0.025_k2", ref = "1", cores = 1, verbose = FALSE)
  
  expect_equal(nrow(result_20[[1]]), 20)
  expect_equal(nrow(result_10[[1]]), 10)
})

# ============================================================================
# TEST 4: Different reference domains
# ============================================================================

test_that("STgradient different reference domains", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb')
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  st_obj = STclust(st_obj, ws = 0.025, ks = c(2, 3), cores = 1, verbose = FALSE)
  
  result_ref1 = STgradient(x = st_obj, samples = 1, topgenes = 20, annot = "stclust_spw0.025_k2", ref = "1", cores = 1, verbose = FALSE)
  result_ref2 = STgradient(x = st_obj, samples = 1, topgenes = 20, annot = "stclust_spw0.025_k2", ref = "2", cores = 1, verbose = FALSE)
  
  expect_s3_class(result_ref1, "list")
  expect_s3_class(result_ref2, "list")
  expect_false(all(result_ref1[[1]]$min_spearman_r == result_ref2[[1]]$min_spearman_r))
})

# ============================================================================
# TEST 5: Matches legacy output
# ============================================================================

test_that("STgradient matches legacy output", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb')
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  st_obj = STclust(st_obj, ws = 0.025, ks = c(2, 3), cores = 1, verbose = FALSE)
  
  result_modern = STgradient(x = st_obj, samples = 1, topgenes = 20, annot = "stclust_spw0.025_k2", ref = "1", distsumm = "min", cores = 1, verbose = FALSE)
  result_legacy = STgradient_legacy(x = st_obj, samples = 1, topgenes = 20, annot = "stclust_spw0.025_k2", ref = "1", distsumm = "min", min_nb = 3, robust = TRUE, nb_dist_thr = c(0.75, 1.25), cores = 1, verbose = FALSE)
  
  expect_equal(nrow(result_modern[[1]]), nrow(result_legacy[[1]]))
  expect_equal(round(result_modern[[1]]$min_spearman_r, 4), round(result_legacy[[1]]$spearman_r, 4), tolerance = 0.1)
})

cat("\n=== STgradient Test Suite Complete ===\n")
