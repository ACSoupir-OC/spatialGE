##
# SThet Legacy Test Suite - Baseline Tests
# 
# These tests establish baseline behavior BEFORE refactoring
# Run these to ensure legacy implementation works correctly
#

library(testthat)
library(devtools)
library(spatialGE)

# Load the package
load_all(".", export_all = TRUE)

# ============================================================================
# TEST 1: Basic SThet_legacy runs without error (Moran's I)
# ============================================================================

test_that("SThet_legacy basic Moran's I", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb', timeout = 60)
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  
  # Run SThet with Moran's I
  result = SThet_legacy(
    x = st_obj,
    genes = c("MLANA", "TP53", "CD37"),
    samples = 1,
    method = "moran",
    cores = 1,
    verbose = FALSE
  )
  
  # Should return STlist
  expect_s4_class(result, "STlist")
  
  # Should have moran_i column in gene_meta
  expect_true("moran_i" %in% colnames(result@gene_meta[[1]]))
  
  # Should have non-NA values for requested genes
  moran_values = result@gene_meta[[1]][result@gene_meta[[1]]$gene %in% c("MLANA", "TP53", "CD37"), "moran_i"]
  expect_true(sum(!is.na(moran_values)) > 0)
})

# ============================================================================
# TEST 2: SThet with Geary's C
# ============================================================================

test_that("SThet Geary's C", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb', timeout = 60)
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  
  # Run SThet with Geary's C
  result = SThet_legacy(
    x = st_obj,
    genes = c("MLANA", "TP53"),
    samples = 1,
    method = "geary",
    cores = 1,
    verbose = FALSE
  )
  
  # Should return STlist
  expect_s4_class(result, "STlist")
  
  # Should have geary_c column in gene_meta
  expect_true("geary_c" %in% colnames(result@gene_meta[[1]]))
  
  # Should have non-NA values
  geary_values = result@gene_meta[[1]][result@gene_meta[[1]]$gene %in% c("MLANA", "TP53"), "geary_c"]
  expect_true(sum(!is.na(geary_values)) > 0)
})

# ============================================================================
# TEST 3: SThet with both Moran's I and Geary's C
# ============================================================================

test_that("SThet both methods", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb', timeout = 60)
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  
  # Run SThet with both methods
  result = SThet_legacy(
    x = st_obj,
    genes = c("MLANA", "TP53"),
    samples = 1,
    method = c("moran", "geary"),
    cores = 1,
    verbose = FALSE
  )
  
  # Should have both columns
  expect_true("moran_i" %in% colnames(result@gene_meta[[1]]))
  expect_true("geary_c" %in% colnames(result@gene_meta[[1]]))
})

# ============================================================================
# TEST 4: SThet with k-nearest neighbors
# ============================================================================

test_that("SThet with k-nearest neighbors", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb', timeout = 60)
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  
  # Run SThet with k-nearest neighbors
  result = SThet_legacy(
    x = st_obj,
    genes = c("MLANA", "TP53"),
    samples = 1,
    method = "moran",
    k = 5,
    cores = 1,
    verbose = FALSE
  )
  
  # Should return STlist
  expect_s4_class(result, "STlist")
  
  # Should have moran_i values
  moran_values = result@gene_meta[[1]][result@gene_meta[[1]]$gene %in% c("MLANA", "TP53"), "moran_i"]
  expect_true(sum(!is.na(moran_values)) > 0)
})

# ============================================================================
# TEST 5: SThet with multiple samples
# ============================================================================

test_that("SThet multiple samples", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb', timeout = 60)
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  
  # Run SThet with multiple samples
  result = SThet_legacy(
    x = st_obj,
    genes = c("MLANA", "TP53"),
    samples = c(1, 2),
    method = "moran",
    cores = 1,
    verbose = FALSE
  )
  
  # Should have results for both samples
  expect_true(length(result@gene_meta) >= 2)
  
  # Both samples should have moran_i
  expect_true("moran_i" %in% colnames(result@gene_meta[[1]]))
  expect_true("moran_i" %in% colnames(result@gene_meta[[2]]))
})

# ============================================================================
# TEST 6: SThet invalid gene
# ============================================================================

test_that("SThet invalid gene warning", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb', timeout = 60)
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  
  # Run SThet with invalid gene - should message but not error
  expect_message({
    result = SThet_legacy(
      x = st_obj,
      genes = c("MLANA", "INVALID_GENE_XYZ"),
      samples = 1,
      method = "moran",
      cores = 1,
      verbose = FALSE
    )
  }, "Not present")
})

# ============================================================================
# TEST 7: SThet no genes error
# ============================================================================

test_that("SThet no genes error", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb', timeout = 60)
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  
  # Run SThet with no genes - should error
  expect_error({
    SThet_legacy(
      x = st_obj,
      genes = NULL,
      samples = 1,
      method = "moran",
      cores = 1,
      verbose = FALSE
    )
  }, "Please enter one or more genes")
})

# ============================================================================
# TEST 8: SThet overwrite parameter
# ============================================================================

test_that("SThet overwrite parameter", {
  # Create test data
  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive = TRUE)
  dir.create(thrane_tmp)
  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), mode = 'wb', timeout = 60)
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir = thrane_tmp)
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names = TRUE, pattern = 'clinical')
  
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  st_obj = transform_data(st_obj)
  
  # Run SThet first time
  result1 = SThet_legacy(
    x = st_obj,
    genes = c("MLANA"),
    samples = 1,
    method = "moran",
    overwrite = TRUE,
    cores = 1,
    verbose = FALSE
  )
  
  moran1 = result1@gene_meta[[1]][result1@gene_meta[[1]]$gene == "MLANA", "moran_i"]
  
  # Run SThet second time with overwrite=FALSE
  result2 = SThet_legacy(
    x = result1,
    genes = c("MLANA"),
    samples = 1,
    method = "moran",
    overwrite = FALSE,
    cores = 1,
    verbose = FALSE
  )
  
  moran2 = result2@gene_meta[[1]][result2@gene_meta[[1]]$gene == "MLANA", "moran_i"]
  
  # Values should be the same (not overwritten)
  expect_equal(moran1, moran2)
})

cat("\n=== SThet Legacy Test Suite Complete ===\n")
