# STdiff Test Suite - Fast Version
# Tests for spatial differential expression analysis
# Optimized for speed while maintaining coverage
#
# Uses shared test data from tests/testthat/setup.R
# Created: 2026-03-19
#

library(testthat)

# ============================================================================
# TEST 1: Basic STdiff - Non-spatial tests
# ============================================================================

test_that("STdiff - Basic non-spatial tests", {
  skip_if_no_data("STdiff")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  cat('TEST 1: Basic non-spatial tests (10 genes)...\n')
  
  result = STdiff(
    x = st_obj,
    samples = 1,
    annot = cluster_col,
    test_type = 't_test',
    sp_topgenes = 0,
    topgenes = 10,
    verbose = 0
  )
  
  expect_true(is.list(result))
  expect_true(length(result) > 0)
  expect_true(is.data.frame(result[[1]]))
  expect_true(nrow(result[[1]]) > 0)
  
  # Check required columns
  expect_true('sample' %in% colnames(result[[1]]))
  expect_true('gene' %in% colnames(result[[1]]))
  expect_true('avg_log2fc' %in% colnames(result[[1]]))
  expect_true('ttest_p_val' %in% colnames(result[[1]]))
  expect_true('adj_p_val' %in% colnames(result[[1]]))
  
  cat('TEST 1: PASSED\n\n')
})

# ============================================================================
# TEST 2: Wilcoxon test type
# ============================================================================

test_that("STdiff - Wilcoxon test", {
  skip_if_no_data("STdiff")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  cat('TEST 2: Wilcoxon test (10 genes)...\n')
  
  result = STdiff(
    x = st_obj,
    samples = 1,
    annot = cluster_col,
    test_type = 'wilcoxon',
    sp_topgenes = 0,
    topgenes = 10,
    verbose = 0
  )
  
  expect_true(length(result) > 0)
  expect_true('wilcox_p_val' %in% colnames(result[[1]]))
  
  cat('TEST 2: PASSED\n\n')
})

# ============================================================================
# TEST 3: Invalid test type error
# ============================================================================

test_that("STdiff - Invalid test type error", {
  skip_if_no_data("STdiff")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  cat('TEST 3: Invalid test type error...\n')
  
  expect_error(
    STdiff(
      x = st_obj,
      samples = 1,
      annot = cluster_col,
      test_type = 'invalid_test',
      sp_topgenes = 0,
      verbose = 0
    ),
    'Test type is one of'
  )
  
  cat('TEST 3: PASSED\n\n')
})

# ============================================================================
# TEST 4: Invalid sp_topgenes error
# ============================================================================

test_that("STdiff - Invalid sp_topgenes error", {
  skip_if_no_data("STdiff")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  cat('TEST 4: Invalid sp_topgenes error...\n')
  
  expect_error(
    STdiff(x = st_obj, samples = 1, annot = cluster_col, sp_topgenes = 1.5, verbose = 0)
  )
  
  expect_error(
    STdiff(x = st_obj, samples = 1, annot = cluster_col, sp_topgenes = -0.1, verbose = 0)
  )
  
  cat('TEST 4: PASSED\n\n')
})

# ============================================================================
# TEST 5: Multiple samples
# ============================================================================

test_that("STdiff - Multiple samples", {
  skip_if_no_data("STdiff")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  cat('TEST 5: Multiple samples (5 genes)...\n')
  
  result = STdiff(
    x = st_obj,
    samples = c(1, 2),
    annot = cluster_col,
    test_type = 't_test',
    sp_topgenes = 0,
    topgenes = 5,
    verbose = 0
  )
  
  expect_true(length(result) >= 1)
  
  cat('TEST 5: PASSED\n\n')
})

# ============================================================================
# TEST 6: Pairwise tests
# ============================================================================

test_that("STdiff - Pairwise tests", {
  skip_if_no_data("STdiff")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  cat('TEST 6: Pairwise tests (5 genes)...\n')
  
  result = STdiff(
    x = st_obj,
    samples = 1,
    annot = cluster_col,
    pairwise = TRUE,
    test_type = 't_test',
    sp_topgenes = 0,
    topgenes = 5,
    verbose = 0
  )
  
  expect_true(length(result) > 0)
  
  cat('TEST 6: PASSED\n\n')
})

# ============================================================================
# TEST 7: P-value adjustment methods
# ============================================================================

test_that("STdiff - P-value adjustment", {
  skip_if_no_data("STdiff")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  cat('TEST 7: P-value adjustment (5 genes)...\n')
  
  result_fdr = STdiff(
    x = st_obj,
    samples = 1,
    annot = cluster_col,
    pval_adj = 'fdr',
    test_type = 't_test',
    sp_topgenes = 0,
    topgenes = 5,
    verbose = 0
  )
  
  result_bonf = STdiff(
    x = st_obj,
    samples = 1,
    annot = cluster_col,
    pval_adj = 'bonferroni',
    test_type = 't_test',
    sp_topgenes = 0,
    topgenes = 5,
    verbose = 0
  )
  
  expect_true(length(result_fdr) > 0)
  expect_true(length(result_bonf) > 0)
  
  cat('TEST 7: PASSED\n\n')
})

# ============================================================================
# TEST 8: Result structure verification
# ============================================================================

test_that("STdiff - Result structure", {
  skip_if_no_data("STdiff")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  cat('TEST 8: Result structure verification...\n')
  
  result = STdiff(
    x = st_obj,
    samples = 1,
    annot = cluster_col,
    test_type = 't_test',
    sp_topgenes = 0,
    topgenes = 10,
    verbose = 0
  )
  
  res_df = result[[1]]
  
  # P-values between 0 and 1
  expect_true(all(res_df$ttest_p_val >= 0, na.rm = TRUE))
  expect_true(all(res_df$ttest_p_val <= 1, na.rm = TRUE))
  
  # Log2FC is numeric
  expect_true(is.numeric(res_df$avg_log2fc))
  
  cat('TEST 8: PASSED\n\n')
})

# ============================================================================
# TEST 9: Invalid sample error
# ============================================================================

test_that("STdiff - Invalid sample error", {
  skip_if_no_data("STdiff")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  cat('TEST 9: Invalid sample error...\n')
  
  expect_error(
    STdiff(
      x = st_obj,
      samples = 'INVALID_SAMPLE',
      annot = cluster_col,
      sp_topgenes = 0,
      verbose = 0
    )
  )
  
  cat('TEST 9: PASSED\n\n')
})

# ============================================================================
# TEST 10: Modern vs Legacy comparison
# ============================================================================

test_that("STdiff - Modern vs legacy", {
  skip_if_no_data("STdiff")
  st_obj = get_test_stlist()
  cluster_col = get_cluster_col()
  
  cat('TEST 10: Modern vs legacy comparison (5 genes)...\n')
  
  result_modern = STdiff(
    x = st_obj,
    samples = 1,
    annot = cluster_col,
    test_type = 't_test',
    sp_topgenes = 0,
    topgenes = 5,
    verbose = 0
  )
  
  result_legacy = STdiff_legacy(
    x = st_obj,
    samples = 1,
    annot = cluster_col,
    test_type = 't_test',
    sp_topgenes = 0,
    topgenes = 5,
    verbose = 0
  )
  
  expect_true(is.list(result_modern))
  expect_true(is.list(result_legacy))
  expect_true(length(result_modern) > 0)
  expect_true(length(result_legacy) > 0)
  
  # Check gene overlap
  genes_modern = sort(result_modern[[1]]$gene)
  genes_legacy = sort(result_legacy[[1]]$gene)
  overlap = length(intersect(genes_modern, genes_legacy))
  expect_true(overlap >= min(length(genes_modern), length(genes_legacy)) * 0.8)
  
  cat('TEST 10: PASSED\n\n')
})

cat('\n========================================\n')
cat('All STdiff tests complete!\n')
cat('========================================\n')
