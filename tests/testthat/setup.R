# ============================================================================
# Test Setup - Create shared test data for all testthat tests
# ============================================================================
#
# This file is sourced once before all tests run. It creates test data that
# is shared across all test files, avoiding redundant downloads and processing.
#
# Created: 2026-03-19
#

# Load required packages
library(testthat)
library(devtools)
library(spatialGE)

# Load the package (from parent directory)
# Note: When run via devtools::test(), package is already loaded
# This is for direct source() execution
if (!"spatialGE" %in% search()) {
  tryCatch({
    load_all("..", export_all = TRUE, quiet = TRUE)
  }, error = function(e) {
    cat('Note: Package already loaded or load failed:', conditionMessage(e), '\n')
  })
}

# Global environment for test data
.test_env = new.env(parent = emptyenv())

# ============================================================================
# Create Test Data
# ============================================================================

cat('Setting up shared test data for STdiff/STenrich/STgradient tests...\n')

# Create temp directory for data
thrane_tmp = tempdir()
unlink(thrane_tmp, recursive = TRUE)
dir.create(thrane_tmp, showWarnings = FALSE)

# Download melanoma_thrane dataset
cat('Downloading test data (melanoma_thrane)...\n')
lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='

tryCatch({
  download.file(lk, destfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), 
                mode = 'wb', timeout = 120, quiet = TRUE)
  
  cat('Extracting test data...\n')
  unzip(zipfile = paste0(thrane_tmp, '/melanoma_thrane.zip'), 
        exdir = thrane_tmp, overwrite = TRUE)
  
  count_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), 
                           full.names = TRUE, pattern = 'counts')
  coord_files = list.files(paste0(thrane_tmp, '/melanoma_thrane'), 
                           full.names = TRUE, pattern = 'mapping')
  clin_file = list.files(paste0(thrane_tmp, '/melanoma_thrane'), 
                         full.names = TRUE, pattern = 'clinical')
  
  cat('Creating STlist object...\n')
  st_obj = STlist(rnacounts = count_files, spotcoords = coord_files, samples = clin_file)
  
  cat('Transforming data...\n')
  st_obj = transform_data(st_obj)
  
  # Run STclust for cluster annotations
  cat('Running STclust for cluster annotations...\n')
  st_obj = STclust(st_obj, samples = 1, w = 0.025, deepSplit = FALSE, verbose = FALSE)
  
  # Get cluster column name (STclust creates column with params in name)
  cluster_col = grep('stclust', colnames(st_obj@spatial_meta[[1]]), value = TRUE)[1]
  
  # Store in test environment
  .test_env$st_obj = st_obj
  .test_env$cluster_col = cluster_col
  .test_env$thrane_tmp = thrane_tmp
  
  cat('Test data ready. Cluster column:', cluster_col, '\n')
  cat('Test samples available:', paste(names(st_obj@spatial_meta), collapse = ', '), '\n\n')
  
}, error = function(e) {
  cat('ERROR setting up test data:', conditionMessage(e), '\n')
  cat('Tests will be skipped.\n\n')
})

# ============================================================================
# Helper Functions
# ============================================================================

#' Get shared test STlist object
#' @return STlist object or NULL if setup failed
get_test_stlist = function() {
  if (exists("st_obj", envir = .test_env)) {
    return(.test_env$st_obj)
  }
  return(NULL)
}

#' Get cluster column name
#' @return column name string or NULL
get_cluster_col = function() {
  if (exists("cluster_col", envir = .test_env)) {
    return(.test_env$cluster_col)
  }
  return(NULL)
}

#' Skip test if test data not available
#' @param test_name Name of test for error message
skip_if_no_data = function(test_name) {
  if (is.null(get_test_stlist())) {
    testthat::skip(paste('No test data available for', test_name))
  }
}

#' Create gene sets for STenrich tests
#' @param st_obj STlist object
#' @return list of gene sets
create_gene_sets = function(st_obj) {
  genes = rownames(st_obj@tr_counts[[1]])
  list(
    GS1 = genes[1:10],
    GS2 = genes[11:20],
    GS3 = genes[21:30]
  )
}

cat('Setup complete. Helper functions available:\n')
cat('  - get_test_stlist()\n')
cat('  - get_cluster_col()\n')
cat('  - skip_if_no_data(test_name)\n')
cat('  - create_gene_sets(st_obj)\n\n')
