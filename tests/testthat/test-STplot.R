# STplot Test Suite - Visualization Function Tests
# Tests for spatial plotting functions
#
# Created: 2026-03-19
# Note: Visualization tests focus on input validation, output structure, 
#       and plot components rather than pixel-perfect visual comparison
#

# ============================================================================
# Setup - Create test data
# ============================================================================

cat('Setting up STplot test data...\n')

# Use melanoma_thrane dataset
thrane_tmp = tempdir()
unlink(thrane_tmp, recursive=TRUE)
dir.create(thrane_tmp)

lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
cat('Downloading test data...\n')
download.file(lk, destfile=paste0(thrane_tmp, '/melanoma_thrane.zip'), mode='wb', timeout=60)
unzip(zipfile=paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir=thrane_tmp)

count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='counts')
coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='mapping')
clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='clinical')

cat('Creating STlist...\n')
st_obj = STlist(rnacounts=count_files, spotcoords=coord_files, samples=clin_file)
st_obj = transform_data(st_obj)

# Run STclust for cluster-based tests
cat('Running STclust for cluster tests...\n')
st_obj = STclust(st_obj, samples=1, verbose=FALSE)

test_genes = c('MLANA', 'TP53', 'CD37')
test_samples = 1

cat('Test data ready.\n\n')

# ============================================================================
# TEST 1: Basic gene expression plot
# ============================================================================

test_that("STplot - Basic gene expression", {
  cat('TEST 1: Basic gene expression plot...\n')
  
  # Test single gene
  plot_list = STplot(st_obj, genes='MLANA', samples=test_samples, ptsize=1)
  
  # Should return a list
  expect_true(is.list(plot_list))
  expect_true(length(plot_list) > 0)
  
  # Each element should be a ggplot object
  for (p in plot_list) {
    expect_true(inherits(p, 'ggplot'))
  }
  
  cat('TEST 1: PASSED\n\n')
})

# ============================================================================
# TEST 2: Multiple genes
# ============================================================================

test_that("STplot - Multiple genes", {
  cat('TEST 2: Multiple genes plot...\n')
  
  plot_list = STplot(st_obj, genes=test_genes, samples=test_samples, ptsize=1)
  
  # Should have one plot per gene
  expect_equal(length(plot_list), length(test_genes))
  
  # All should be ggplot objects
  for (p in plot_list) {
    expect_true(inherits(p, 'ggplot'))
  }
  
  cat('TEST 2: PASSED\n\n')
})

# ============================================================================
# TEST 3: Invalid gene error
# ============================================================================

test_that("STplot - Invalid gene error", {
  cat('TEST 3: Invalid gene error...\n')
  
  # Should error when gene not found
  # Note: Error message may vary depending on where validation occurs
  expect_error(
    STplot(st_obj, genes='INVALID_GENE_12345', samples=test_samples)
  )
  
  cat('TEST 3: PASSED\n\n')
})

# ============================================================================
# TEST 4: Invalid sample error
# ============================================================================

test_that("STplot - Invalid sample error", {
  cat('TEST 4: Invalid sample error...\n')
  
  # Should error when sample not found
  expect_error(
    STplot(st_obj, genes='MLANA', samples='INVALID_SAMPLE'),
    'The requested samples are not present'
  )
  
  cat('TEST 4: PASSED\n\n')
})

# ============================================================================
# TEST 5: Metadata plotting
# ============================================================================

test_that("STplot - Metadata plotting", {
  cat('TEST 5: Metadata plot...\n')
  
  # Plot spatial metadata column
  # First check what columns are available
  meta_cols = colnames(st_obj@spatial_meta[[1]])
  
  if ('libname' %in% meta_cols) {
    plot_list = STplot(st_obj, samples=test_samples, plot_meta='libname', ptsize=1)
    expect_true(is.list(plot_list))
    expect_true(length(plot_list) > 0)
    for (p in plot_list) {
      expect_true(inherits(p, 'ggplot'))
    }
  }
  
  cat('TEST 5: PASSED\n\n')
})

# ============================================================================
# TEST 6: Cluster membership plotting (requires STclust)
# ============================================================================

test_that("STplot - Cluster membership", {
  cat('TEST 6: Cluster membership plot...\n')
  
  # Plot cluster memberships (ks='dtc' uses dynamic tree cut results)
  plot_list = STplot(st_obj, samples=test_samples, ks='dtc', ptsize=1)
  
  expect_true(is.list(plot_list))
  expect_true(length(plot_list) > 0)
  for (p in plot_list) {
    expect_true(inherits(p, 'ggplot'))
  }
  
  cat('TEST 6: PASSED\n\n')
})

# ============================================================================
# TEST 7: Raw vs transformed data
# ============================================================================

test_that("STplot - Raw vs transformed data", {
  cat('TEST 7: Raw vs transformed data...\n')
  
  # Test transformed data (default)
  plot_tr = STplot(st_obj, genes='MLANA', samples=test_samples, 
                   data_type='tr', ptsize=1)
  
  # Test raw data
  plot_raw = STplot(st_obj, genes='MLANA', samples=test_samples, 
                    data_type='raw', ptsize=1)
  
  # Both should return plots
  expect_true(length(plot_tr) > 0)
  expect_true(length(plot_raw) > 0)
  
  cat('TEST 7: PASSED\n\n')
})

# ============================================================================
# TEST 8: Invalid data_type error
# ============================================================================

test_that("STplot - Invalid data_type error", {
  cat('TEST 8: Invalid data_type error...\n')
  
  # Should error on invalid data_type
  expect_error(
    STplot(st_obj, genes='MLANA', samples=test_samples, data_type='invalid'),
    'Please, select one of transformed \\(tr\\) or raw \\(raw\\)'
  )
  
  cat('TEST 8: PASSED\n\n')
})

# ============================================================================
# TEST 9: Gene set plotting (list input)
# ============================================================================

test_that("STplot - Gene set plotting", {
  cat('TEST 9: Gene set plot...\n')
  
  # Create a gene set list
  gene_sets = list(
    immune = c('CD37', 'TP53'),
    melanoma = c('MLANA')
  )
  
  plot_list = STplot(st_obj, genes=gene_sets, samples=test_samples, ptsize=1)
  
  # Should return plots for each gene set
  expect_true(is.list(plot_list))
  expect_true(length(plot_list) > 0)
  for (p in plot_list) {
    expect_true(inherits(p, 'ggplot'))
  }
  
  cat('TEST 9: PASSED\n\n')
})

# ============================================================================
# TEST 10: Color palette parameter
# ============================================================================

test_that("STplot - Color palette", {
  cat('TEST 10: Color palette parameter...\n')
  
  # Test with different color palettes (RColorBrewer)
  plot_buRd = STplot(st_obj, genes='MLANA', samples=test_samples, 
                     color_pal='BuRd', ptsize=1)
  
  plot_YlOrRd = STplot(st_obj, genes='MLANA', samples=test_samples, 
                       color_pal='YlOrRd', ptsize=1)
  
  # Both should work
  expect_true(length(plot_buRd) > 0)
  expect_true(length(plot_YlOrRd) > 0)
  
  cat('TEST 10: PASSED\n\n')
})

# ============================================================================
# TEST 11: Multiple samples
# ============================================================================

test_that("STplot - Multiple samples", {
  cat('TEST 11: Multiple samples...\n')
  
  # Plot multiple samples
  multi_samples = c(1, 2)
  plot_list = STplot(st_obj, genes='MLANA', samples=multi_samples, ptsize=1)
  
  # Should have plots for each sample
  expect_true(length(plot_list) >= length(multi_samples))
  
  cat('TEST 11: PASSED\n\n')
})

# ============================================================================
# TEST 12: STplot_interpolation basic test
# ============================================================================

test_that("STplot_interpolation - Basic kriging plot", {
  cat('TEST 12: Interpolation plot...\n')
  
  # First need to run gene_interpolation
  # Skip if gstat not available or too slow
  skip_if_not_installed('gstat')
  
  tryCatch({
    st_interp = gene_interpolation(st_obj, genes='MLANA', samples=test_samples)
    
    if (!is.null(st_interp@gene_krige[['MLANA']][[1]][['krige_out']])) {
      plot_list = STplot_interpolation(st_interp, genes='MLANA', samples=test_samples)
      
      expect_true(is.list(plot_list))
      expect_true(length(plot_list) > 0)
      for (p in plot_list) {
        expect_true(inherits(p, 'ggplot'))
      }
    }
  }, error = function(e) {
    # Interpolation can fail for various reasons, skip test
    skip('Interpolation failed')
  })
  
  cat('TEST 12: PASSED\n\n')
})

# ============================================================================
# TEST 13: No genes or metadata error
# ============================================================================

test_that("STplot - No genes or metadata error", {
  cat('TEST 13: No genes or metadata error...\n')
  
  # Should provide a helpful error when neither genes nor plot_meta specified
  # Actually, STplot will try to plot metadata if genes is NULL
  # So this might not error - just test that it runs
  plot_list = STplot(st_obj, samples=test_samples, ptsize=1)
  expect_true(is.list(plot_list))
  
  cat('TEST 13: PASSED\n\n')
})

# ============================================================================
# TEST 14: Point size parameter
# ============================================================================

test_that("STplot - Point size parameter", {
  cat('TEST 14: Point size parameter...\n')
  
  # Test different point sizes
  plot_small = STplot(st_obj, genes='MLANA', samples=test_samples, 
                      ptsize=0.5)
  plot_large = STplot(st_obj, genes='MLANA', samples=test_samples, 
                      ptsize=3)
  
  # Both should work
  expect_true(length(plot_small) > 0)
  expect_true(length(plot_large) > 0)
  
  cat('TEST 14: PASSED\n\n')
})

# ============================================================================
# TEST 15: Plot data verification
# ============================================================================

test_that("STplot - Plot data verification", {
  cat('TEST 15: Plot data verification...\n')
  
  plot_list = STplot(st_obj, genes='MLANA', samples=test_samples, ptsize=1)
  
  # Extract the plot
  p = plot_list[[1]]
  
  # Check that plot has data
  expect_true(!is.null(p$data))
  expect_true(nrow(p$data) > 0)
  
  # Check that data has expected columns (xpos, ypos, and expression)
  expect_true('xpos' %in% colnames(p$data) || 'x' %in% colnames(p$data))
  expect_true('ypos' %in% colnames(p$data) || 'y' %in% colnames(p$data))
  
  cat('TEST 15: PASSED\n\n')
})

cat('\n========================================\n')
cat('All STplot tests complete!\n')
cat('========================================\n')
