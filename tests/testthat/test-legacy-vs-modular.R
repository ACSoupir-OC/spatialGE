#' Test Legacy vs New Modular Function Equivalence
#'
#' This test suite verifies that the refactored modular functions
#' produce identical (or numerically equivalent) outputs to the legacy
#' monolithic implementations.
#'
#' @import testthat

test_that("STclust legacy vs modular produce identical cluster assignments", {
  skip_if_not_installed("testthat")
  
  # Load test data - TNBC Bassiouni dataset (Visium)
  data_dir <- file.path("data", "tnbc_bassiouni")
  skip_if_not(dir.exists(data_dir), "TNBC data not found")
  
  # Get sample directories
  sample_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
  sample_dirs <- grep("sample_", sample_dirs, value = TRUE)[1:2]  # Use 2 samples for speed
  skip_if(length(sample_dirs) < 2, "Not enough samples")
  
  clinical <- file.path(data_dir, "bassiouni_clinical.csv")
  
  # Create STlist object (Visium format)
  tnbc <- STlist(rnacounts = sample_dirs, samples = clinical)
  tnbc <- transform_data(tnbc)
  
  # Legacy approach
  set.seed(123)
  tnbc_legacy <- tnbc
  tnbc_legacy <- STclust_legacy(tnbc_legacy, w = 0.5, k = 5, deepSplit = 2)
  
  # New modular approach
  set.seed(123)
  tnbc_new <- tnbc
  tnbc_new <- STclust(tnbc_new, w = 0.5, k = 5, deepSplit = 2)
  
  # Compare cluster assignments
  expect_equal(
    tnbc_legacy@spatial_meta[[1]]$cluster,
    tnbc_new@spatial_meta[[1]]$cluster,
    tolerance = 0,
    info = "Cluster assignments should be identical"
  )
})

test_that("SThet legacy vs modular produce identical Moran's I", {
  data_dir <- file.path("data", "tnbc_bassiouni")
  skip_if_not(dir.exists(data_dir), "TNBC data not found")
  
  sample_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
  sample_dirs <- grep("sample_", sample_dirs, value = TRUE)[1:2]
  skip_if(length(sample_dirs) < 2, "Not enough samples")
  
  clinical <- file.path(data_dir, "bassiouni_clinical.csv")
  
  tnbc <- STlist(rnacounts = sample_dirs, samples = clinical)
  tnbc <- transform_data(tnbc)
  
  # Legacy approach
  tnbc_legacy <- tnbc
  tnbc_legacy <- SThet_legacy(tnbc_legacy, genes = c('CD3E', 'MS4A1'))
  
  # New modular approach
  tnbc_new <- tnbc
  tnbc_new <- SThet(tnbc_new, genes = c('CD3E', 'MS4A1'))
  
  # Compare Moran's I values
  legacy_moran <- tnbc_legacy@gene_meta[[1]]$moranI
  new_moran <- tnbc_new@gene_meta[[1]]$moranI
  
  expect_equal(
    legacy_moran,
    new_moran,
    tolerance = 1e-6,
    info = "Moran's I values should be identical"
  )
  
  # Compare p-values
  legacy_pval <- tnbc_legacy@gene_meta[[1]]$moranI_pval
  new_pval <- tnbc_new@gene_meta[[1]]$moranI_pval
  
  expect_equal(
    legacy_pval,
    new_pval,
    tolerance = 1e-6,
    info = "Moran's I p-values should be identical"
  )
})

test_that("STdiff legacy vs modular produce equivalent p-values (non-spatial)", {
  data_dir <- file.path("data", "tnbc_bassiouni")
  skip_if_not(dir.exists(data_dir), "TNBC data not found")
  
  sample_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
  sample_dirs <- grep("sample_", sample_dirs, value = TRUE)[1:2]
  skip_if(length(sample_dirs) < 2, "Not enough samples")
  
  clinical <- file.path(data_dir, "bassiouni_clinical.csv")
  
  tnbc <- STlist(rnacounts = sample_dirs, samples = clinical)
  tnbc <- transform_data(tnbc)
  
  # Legacy approach (non-spatial only for speed)
  set.seed(456)
  result_legacy <- STdiff_legacy(
    x = tnbc,
    annot = 'tissue_type',
    test_type = 'wilcoxon',
    sp_topgenes = 0  # Skip spatial models for speed
  )
  
  # New modular approach
  set.seed(456)
  tnbc2 <- tnbc  # Fresh copy
  result_new <- STdiff(
    x = tnbc2,
    annot = 'tissue_type',
    test_type = 'wilcoxon',
    sp_topgenes = 0
  )
  
  # Compare number of genes tested
  expect_equal(
    nrow(result_legacy[[1]]),
    nrow(result_new[[1]]),
    info = "Same number of genes should be tested"
  )
  
  # Compare p-values (allow small numerical differences)
  legacy_pvals <- result_legacy[[1]]$wilcox_p_val
  new_pvals <- result_new[[1]]$wilcox_p_val
  
  # Remove NA values
  legacy_pvals <- legacy_pvals[!is.na(legacy_pvals)]
  new_pvals <- new_pvals[!is.na(new_pvals)]
  
  # Check correlation (should be very high)
  if(length(legacy_pvals) > 2) {
    correlation <- cor(legacy_pvals, new_pvals, method = "spearman")
    expect_gt(correlation, 0.99, info = "P-values should be highly correlated")
    
    # Check mean absolute difference
    mean_diff <- mean(abs(legacy_pvals - new_pvals), na.rm = TRUE)
    expect_lt(mean_diff, 0.01, info = "Mean p-value difference should be small")
  }
})

test_that("STgradient legacy vs modular produce equivalent gradient statistics", {
  data_dir <- file.path("data", "tnbc_bassiouni")
  skip_if_not(dir.exists(data_dir), "TNBC data not found")
  
  sample_dirs <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
  sample_dirs <- grep("sample_", sample_dirs, value = TRUE)[1:2]
  skip_if(length(sample_dirs) < 2, "Not enough samples")
  
  clinical <- file.path(data_dir, "bassiouni_clinical.csv")
  
  tnbc <- STlist(rnacounts = sample_dirs, samples = clinical)
  tnbc <- transform_data(tnbc)
  
  # Legacy approach
  set.seed(321)
  tnbc_legacy <- tnbc
  tnbc_legacy <- STgradient_legacy(
    x = tnbc_legacy,
    annot = 'tissue_type',
    ref = 'Tumor',
    distsumm = 'min',
    min_nb = 3
  )
  
  # New modular approach
  set.seed(321)
  tnbc_new <- tnbc
  tnbc_new <- STgradient(
    x = tnbc_new,
    annot = 'tissue_type',
    ref = 'Tumor',
    distsumm = 'min',
    min_nb = 3
  )
  
  # Compare gradient coefficients
  legacy_coef <- tnbc_legacy@gene_meta[[1]]$gradient_coef
  new_coef <- tnbc_new@gene_meta[[1]]$gradient_coef
  
  # Remove NA values
  legacy_coef <- legacy_coef[!is.na(legacy_coef)]
  new_coef <- new_coef[!is.na(new_coef)]
  
  if(length(legacy_coef) > 2) {
    # Check correlation
    correlation <- cor(legacy_coef, new_coef, method = "pearson")
    expect_gt(correlation, 0.99, info = "Gradient coefficients should be highly correlated")
    
    # Compare p-values
    legacy_pvals <- tnbc_legacy@gene_meta[[1]]$gradient_pval
    new_pvals <- tnbc_new@gene_meta[[1]]$gradient_pval
    
    legacy_pvals <- legacy_pvals[!is.na(legacy_pvals)]
    new_pvals <- new_pvals[!is.na(new_pvals)]
    
    correlation_pval <- cor(legacy_pvals, new_pvals, method = "spearman")
    expect_gt(correlation_pval, 0.95, info = "Gradient p-values should be correlated")
  }
})
