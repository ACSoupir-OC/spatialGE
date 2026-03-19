# SThet Comprehensive Test Suite - Modern vs Legacy Comparison
# Tests the modular SThet implementation against legacy baseline
#
# Created: 2026-03-19
# Pattern: Create test data → Run both versions → Compare results
#

# ============================================================================
# Setup - Create test data
# ============================================================================

cat('Setting up test data...\n')

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

test_genes = c('MLANA', 'TP53', 'CD37', 'CD8A', 'CD4')
test_samples = 1

cat('Test data ready.\n\n')

# ============================================================================
# TEST 1: Basic Moran's I comparison
# ============================================================================

test_that("SThet modern vs legacy - Moran's I", {
  cat('TEST 1: Moran\'s I comparison...\n')
  
  result_modern = SThet(st_obj, genes=test_genes, samples=test_samples, 
                        method='moran', cores=1, verbose=FALSE)
  result_legacy = SThet_legacy(st_obj, genes=test_genes, samples=test_samples, 
                               method='moran', cores=1, verbose=FALSE)
  
  for (gene in test_genes) {
    val_mod = result_modern@gene_meta[[1]][
      result_modern@gene_meta[[1]]$gene == gene, 'moran_i']
    val_leg = result_legacy@gene_meta[[1]][
      result_legacy@gene_meta[[1]]$gene == gene, 'moran_i']
    
    expect_equal(as.numeric(val_mod), as.numeric(val_leg), tolerance=1e-10)
  }
  
  cat('TEST 1: PASSED\n\n')
})

# ============================================================================
# TEST 2: Geary's C comparison
# ============================================================================

test_that("SThet modern vs legacy - Geary's C", {
  cat('TEST 2: Geary\'s C comparison...\n')
  
  # Reset gene_meta to avoid overwrite issues
  st_obj2 = st_obj
  
  result_modern = SThet(st_obj2, genes=test_genes, samples=test_samples, 
                        method='geary', cores=1, verbose=FALSE)
  result_legacy = SThet_legacy(st_obj2, genes=test_genes, samples=test_samples, 
                               method='geary', cores=1, verbose=FALSE)
  
  for (gene in test_genes) {
    val_mod = result_modern@gene_meta[[1]][
      result_modern@gene_meta[[1]]$gene == gene, 'geary_c']
    val_leg = result_legacy@gene_meta[[1]][
      result_legacy@gene_meta[[1]]$gene == gene, 'geary_c']
    
    expect_equal(as.numeric(val_mod), as.numeric(val_leg), tolerance=1e-10)
  }
  
  cat('TEST 2: PASSED\n\n')
})

# ============================================================================
# TEST 3: Both methods together
# ============================================================================

test_that("SThet modern vs legacy - Both methods", {
  cat('TEST 3: Both methods comparison...\n')
  
  st_obj3 = st_obj
  
  result_modern = SThet(st_obj3, genes=test_genes[1:3], samples=test_samples, 
                        method=c('moran', 'geary'), cores=1, verbose=FALSE)
  result_legacy = SThet_legacy(st_obj3, genes=test_genes[1:3], samples=test_samples, 
                               method=c('moran', 'geary'), cores=1, verbose=FALSE)
  
  for (gene in test_genes[1:3]) {
    # Check Moran's I
    val_mod_moran = result_modern@gene_meta[[1]][
      result_modern@gene_meta[[1]]$gene == gene, 'moran_i']
    val_leg_moran = result_legacy@gene_meta[[1]][
      result_legacy@gene_meta[[1]]$gene == gene, 'moran_i']
    expect_equal(as.numeric(val_mod_moran), as.numeric(val_leg_moran), tolerance=1e-10)
    
    # Check Geary's C
    val_mod_geary = result_modern@gene_meta[[1]][
      result_modern@gene_meta[[1]]$gene == gene, 'geary_c']
    val_leg_geary = result_legacy@gene_meta[[1]][
      result_legacy@gene_meta[[1]]$gene == gene, 'geary_c']
    expect_equal(as.numeric(val_mod_geary), as.numeric(val_leg_geary), tolerance=1e-10)
  }
  
  cat('TEST 3: PASSED\n\n')
})

# ============================================================================
# TEST 4: k-nearest neighbors comparison
# ============================================================================

test_that("SThet modern vs legacy - k-nearest neighbors", {
  cat('TEST 4: k-nearest neighbors comparison...\n')
  
  st_obj4 = st_obj
  
  result_modern = SThet(st_obj4, genes=test_genes[1:2], samples=test_samples, 
                        method='moran', k=5, cores=1, verbose=FALSE)
  result_legacy = SThet_legacy(st_obj4, genes=test_genes[1:2], samples=test_samples, 
                               method='moran', k=5, cores=1, verbose=FALSE)
  
  for (gene in test_genes[1:2]) {
    val_mod = result_modern@gene_meta[[1]][
      result_modern@gene_meta[[1]]$gene == gene, 'moran_i']
    val_leg = result_legacy@gene_meta[[1]][
      result_legacy@gene_meta[[1]]$gene == gene, 'moran_i']
    
    expect_equal(as.numeric(val_mod), as.numeric(val_leg), tolerance=1e-10)
  }
  
  cat('TEST 4: PASSED\n\n')
})

# ============================================================================
# TEST 5: Multiple samples
# ============================================================================

test_that("SThet modern vs legacy - Multiple samples", {
  cat('TEST 5: Multiple samples comparison...\n')
  
  st_obj5 = st_obj
  multi_samples = c(1, 2)
  
  result_modern = SThet(st_obj5, genes=test_genes[1:2], samples=multi_samples, 
                        method='moran', cores=1, verbose=FALSE)
  result_legacy = SThet_legacy(st_obj5, genes=test_genes[1:2], samples=multi_samples, 
                               method='moran', cores=1, verbose=FALSE)
  
  for (samp in multi_samples) {
    for (gene in test_genes[1:2]) {
      val_mod = result_modern@gene_meta[[samp]][
        result_modern@gene_meta[[samp]]$gene == gene, 'moran_i']
      val_leg = result_legacy@gene_meta[[samp]][
        result_legacy@gene_meta[[samp]]$gene == gene, 'moran_i']
      
      expect_equal(as.numeric(val_mod), as.numeric(val_leg), tolerance=1e-10)
    }
  }
  
  cat('TEST 5: PASSED\n\n')
})

# ============================================================================
# TEST 6: Invalid gene warning
# ============================================================================

test_that("SThet modern vs legacy - Invalid gene warning", {
  cat('TEST 6: Invalid gene warning...\n')
  
  st_obj6 = st_obj
  genes_with_invalid = c('MLANA', 'INVALID_GENE_12345', 'TP53')
  
  # Capture messages using testthat's capture_messages
  msg_modern = testthat::capture_messages({
    result_modern = SThet(st_obj6, genes=genes_with_invalid, samples=test_samples, 
                          method='moran', cores=1, verbose=FALSE)
  })
  
  msg_legacy = testthat::capture_messages({
    result_legacy = SThet_legacy(st_obj6, genes=genes_with_invalid, samples=test_samples, 
                                 method='moran', cores=1, verbose=FALSE)
  })
  
  # Both should warn about invalid gene
  expect_true(any(grepl('INVALID_GENE_12345', msg_modern)))
  expect_true(any(grepl('INVALID_GENE_12345', msg_legacy)))
  
  # Valid genes should still work
  val_mod = result_modern@gene_meta[[1]][
    result_modern@gene_meta[[1]]$gene == 'MLANA', 'moran_i']
  val_leg = result_legacy@gene_meta[[1]][
    result_legacy@gene_meta[[1]]$gene == 'MLANA', 'moran_i']
  
  expect_equal(as.numeric(val_mod), as.numeric(val_leg), tolerance=1e-10)
  
  cat('TEST 6: PASSED\n\n')
})

# ============================================================================
# TEST 7: No genes error
# ============================================================================

test_that("SThet modern vs legacy - No genes error", {
  cat('TEST 7: No genes error...\n')
  
  expect_error(
    SThet(st_obj, genes=NULL, samples=test_samples, method='moran', cores=1, verbose=FALSE),
    'Please enter one or more genes'
  )
  
  expect_error(
    SThet_legacy(st_obj, genes=NULL, samples=test_samples, method='moran', cores=1, verbose=FALSE),
    'Please enter one or more genes'
  )
  
  cat('TEST 7: PASSED\n\n')
})

# ============================================================================
# TEST 8: Overwrite parameter
# ============================================================================

test_that("SThet modern vs legacy - Overwrite parameter", {
  cat('TEST 8: Overwrite parameter...\n')
  
  st_obj8 = st_obj
  
  # First run
  result1_modern = SThet(st_obj8, genes='MLANA', samples=test_samples, 
                         method='moran', cores=1, verbose=FALSE)
  result1_legacy = SThet_legacy(st_obj8, genes='MLANA', samples=test_samples, 
                                method='moran', cores=1, verbose=FALSE)
  
  # Second run with overwrite=FALSE (should not change values)
  val1_mod = result1_modern@gene_meta[[1]][
    result1_modern@gene_meta[[1]]$gene == 'MLANA', 'moran_i']
  val1_leg = result1_legacy@gene_meta[[1]][
    result1_legacy@gene_meta[[1]]$gene == 'MLANA', 'moran_i']
  
  # Modify the value manually
  result1_modern@gene_meta[[1]][
    result1_modern@gene_meta[[1]]$gene == 'MLANA', 'moran_i'] = 999
  result1_legacy@gene_meta[[1]][
    result1_legacy@gene_meta[[1]]$gene == 'MLANA', 'moran_i'] = 999
  
  # Run again with overwrite=FALSE (should keep 999)
  result2_modern = SThet(result1_modern, genes='MLANA', samples=test_samples, 
                         method='moran', overwrite=FALSE, cores=1, verbose=FALSE)
  result2_legacy = SThet_legacy(result1_legacy, genes='MLANA', samples=test_samples, 
                                method='moran', overwrite=FALSE, cores=1, verbose=FALSE)
  
  val2_mod = result2_modern@gene_meta[[1]][
    result2_modern@gene_meta[[1]]$gene == 'MLANA', 'moran_i']
  val2_leg = result2_legacy@gene_meta[[1]][
    result2_legacy@gene_meta[[1]]$gene == 'MLANA', 'moran_i']
  
  # Should still be 999 (not overwritten)
  expect_equal(as.numeric(val2_mod), 999, tolerance=1e-10)
  expect_equal(as.numeric(val2_leg), 999, tolerance=1e-10)
  
  cat('TEST 8: PASSED\n\n')
})

# ============================================================================
# TEST 9: SThet_invdist_test comparison with legacy
# ============================================================================

test_that("SThet_invdist_test modern vs legacy", {
  cat('TEST 9: invdist_test comparison...\n')
  
  st_obj9 = st_obj
  
  result_modern = SThet_invdist_test(st_obj9, genes=test_genes[1:2], samples=test_samples, 
                                     method='moran', cores=1)
  result_legacy = SThet_invdist_test_legacy(st_obj9, genes=test_genes[1:2], samples=test_samples, 
                                            method='moran', cores=1)
  
  for (gene in test_genes[1:2]) {
    val_mod = result_modern@gene_meta[[1]][
      result_modern@gene_meta[[1]]$gene == gene, 'moran_i']
    val_leg = result_legacy@gene_meta[[1]][
      result_legacy@gene_meta[[1]]$gene == gene, 'moran_i']
    
    expect_equal(as.numeric(val_mod), as.numeric(val_leg), tolerance=1e-10)
  }
  
  cat('TEST 9: PASSED\n\n')
})

cat('\n========================================\n')
cat('All SThet comprehensive tests complete!\n')
cat('========================================\n')
