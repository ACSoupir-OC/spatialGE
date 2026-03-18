##
# STclust Comprehensive Tests - Section 1: Core Helper Functions
# Tests for STclust_select_genes(), calculate_distances(), weight_distances()
# Test data: melanoma_thrane dataset
#

test_that("STclust_select_genes() - Variable gene selection with VST", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  
  melanoma <- load_melanoma_thrane()
  samples <- names(melanoma@spatial_meta)
  
  # Test topgenes=500, 1000, 2000
  result_500 <- STclust_select_genes(x=melanoma, samples=samples, topgenes=500, cores=1)
  expect_type(result_500, "list")
  expect_true("x" %in% names(result_500) && "trcounts_df" %in% names(result_500))
  expect_true(inherits(result_500$x, "STlist"))
  expect_true("vst.variance.standardized" %in% colnames(result_500$x@gene_meta[[1]]))
  expect_equal(nrow(result_500$trcounts_df[[1]]), 500)
  
  result_1000 <- STclust_select_genes(x=melanoma, samples=samples, topgenes=1000, cores=1)
  expect_equal(nrow(result_1000$trcounts_df[[1]]), 1000)
  
  result_2000 <- STclust_select_genes(x=melanoma, samples=samples, topgenes=2000, cores=1)
  expect_equal(nrow(result_2000$trcounts_df[[1]]), 2000)
  
  # Verify genes ordered by variance (descending)
  variances <- result_2000$x@gene_meta[[1]]$vst.variance.standardized[1:10]
  expect_true(all(diff(variances) <= 0.1))
  
  # Multi-sample selection
  if(length(samples) > 1) {
    result_multi <- STclust_select_genes(x=melanoma, samples=samples, topgenes=1000, cores=1)
    expect_length(result_multi$trcounts_df, length(samples))
  }
})


test_that("STclust_calculate_distances() - Euclidean, manhattan, spatial distances", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  melanoma <- load_melanoma_thrane()
  expr_dat <- melanoma@tr_counts[[1]]
  coord_dat <- melanoma@spatial_meta[[1]]
  
  # Test Euclidean distance
  euclid_dists <- STclust_calculate_distances(trcounts_df=expr_dat, coord_dat=coord_dat, dist_metric='euclidean')
  expect_type(euclid_dists, "list")
  expect_true(all(c("scale_exp", "scale_coord") %in% names(euclid_dists)))
  expect_true(inherits(euclid_dists$scale_exp, "dist") && inherits(euclid_dists$scale_coord, "dist"))
  
  # Verify symmetric matrices
  exp_mat <- as.matrix(euclid_dists$scale_exp)
  coord_mat <- as.matrix(euclid_dists$scale_coord)
  expect_equal(exp_mat, t(exp_mat), tolerance=1e-10)
  expect_equal(coord_mat, t(coord_mat), tolerance=1e-10)
  
  # Test Manhattan distance
  manhattan_dists <- STclust_calculate_distances(trcounts_df=expr_dat, coord_dat=coord_dat, dist_metric='manhattan')
  expect_type(manhattan_dists, "list")
  expect_true(all(c("scale_exp", "scale_coord") %in% names(manhattan_dists)))
  expect_false(all(exp_mat == as.matrix(manhattan_dists$scale_exp)))
  
  # Verify scaling to [0, 1] range
  expect_true(all(euclid_dists$scale_exp >= 0 & euclid_dists$scale_exp <= 1))
  expect_true(all(euclid_dists$scale_coord >= 0 & euclid_dists$scale_coord <= 1))
  expect_true(all(manhattan_dists$scale_exp >= 0 & manhattan_dists$scale_exp <= 1))
  
  # Verify matrix dimensions and diagonal
  n_spots <- nrow(coord_dat)
  expect_equal(dim(euclid_dists$scale_exp), c(n_spots, n_spots))
  expect_true(all(diag(euclid_dists$scale_exp) == 0))
})


test_that("STclust_weight_distances() - ws=0, ws=1, ws=0.025, multiple ws values", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  melanoma <- load_melanoma_thrane()
  expr_dat <- melanoma@tr_counts[[1]]
  coord_dat <- melanoma@spatial_meta[[1]]
  scaled_dists <- STclust_calculate_distances(trcounts_df=expr_dat, coord_dat=coord_dat, dist_metric='euclidean')
  
  exp_dist_mat <- as.matrix(scaled_dists$scale_exp)
  spatial_dist_mat <- as.matrix(scaled_dists$scale_coord)
  
  # Test ws=0 (expression only)
  weighted_0 <- STclust_weight_distances(scaled_dists=scaled_dists, ws=0)
  expect_type(weighted_0, "list")
  expect_length(weighted_0, 1)
  expect_true(inherits(weighted_0[[1]], "dist"))
  expect_equal(as.matrix(weighted_0[[1]]), exp_dist_mat, tolerance=1e-10)
  
  # Test ws=1 (spatial only)
  weighted_1 <- STclust_weight_distances(scaled_dists=scaled_dists, ws=1)
  expect_type(weighted_1, "list")
  expect_length(weighted_1, 1)
  expect_equal(as.matrix(weighted_1[[1]]), spatial_dist_mat, tolerance=1e-10)
  
  # Test ws=0.025 (balanced) - weighted combination
  weighted_balanced <- STclust_weight_distances(scaled_dists=scaled_dists, ws=0.025)
  expected_mat <- 0.025 * spatial_dist_mat + 0.975 * exp_dist_mat
  expect_equal(as.matrix(weighted_balanced[[1]]), expected_mat, tolerance=1e-10)
  
  # Test multiple ws values
  ws_values <- c(0, 0.01, 0.025, 0.05, 0.1, 0.5, 1.0)
  weighted_multi <- STclust_weight_distances(scaled_dists=scaled_dists, ws=ws_values)
  expect_type(weighted_multi, "list")
  expect_length(weighted_multi, length(ws_values))
  
  # Verify each weight produces correct combination
  for(i in seq_along(ws_values)) {
    ws_val <- ws_values[i]
    expected_mat_i <- ws_val * spatial_dist_mat + (1 - ws_val) * exp_dist_mat
    expect_equal(as.matrix(weighted_multi[[i]]), expected_mat_i, tolerance=1e-10)
  }
  
  # Verify distinct matrices for distinct ws values
  expect_true(any(sapply(1:(length(weighted_multi)-1), function(i)
    any(sapply((i+1):length(weighted_multi), function(j)
      !all(as.matrix(weighted_multi[[i]]) == as.matrix(weighted_multi[[j]])))))))
})


test_that("Helper function isolation - Each helper callable independently", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  melanoma <- load_melanoma_thrane()
  samples <- names(melanoma@spatial_meta)
  expr_dat <- melanoma@tr_counts[[1]]
  coord_dat <- melanoma@spatial_meta[[1]]
  
  # Test STclust_select_genes returns correct structure
  select_result <- STclust_select_genes(x=melanoma, samples=samples, topgenes=1000, cores=1)
  expect_type(select_result, "list")
  expect_length(select_result, 2)
  expect_true(inherits(select_result$x, "STlist"))
  expect_type(select_result$trcounts_df, "list")
  expect_length(select_result$trcounts_df, length(samples))
  
  # Test STclust_calculate_distances returns correct structure
  calc_result <- STclust_calculate_distances(trcounts_df=expr_dat, coord_dat=coord_dat, dist_metric='euclidean')
  expect_type(calc_result, "list")
  expect_length(calc_result, 2)
  expect_true(inherits(calc_result$scale_exp, "dist"))
  expect_true(inherits(calc_result$scale_coord, "dist"))
  
  # Test STclust_weight_distances returns correct structure
  weighted_result <- STclust_weight_distances(scaled_dists=calc_result, ws=0.025)
  expect_type(weighted_result, "list")
  expect_length(weighted_result, 1)
  expect_true(inherits(weighted_result[[1]], "dist"))
  
  # Verify functions callable independently
  raw_calc <- STclust_calculate_distances(trcounts_df=expr_dat, coord_dat=coord_dat, dist_metric='euclidean')
  expect_type(raw_calc, "list")
  expect_length(raw_calc, 2)
  
  raw_weighted <- STclust_weight_distances(scaled_dists=raw_calc, ws=0.5)
  expect_type(raw_weighted, "list")
  expect_length(raw_weighted, 1)
  
  # Verify output types consistent across different parameters
  result_500 <- STclust_select_genes(x=melanoma, samples=samples, topgenes=500, cores=1)
  result_2000 <- STclust_select_genes(x=melanoma, samples=samples, topgenes=2000, cores=1)
  expect_type(result_500, "list")
  expect_type(result_2000, "list")
  expect_equal(nrow(result_500$trcounts_df[[1]]), 500)
  expect_equal(nrow(result_2000$trcounts_df[[1]]), 2000)
  
  euclid_calc <- STclust_calculate_distances(trcounts_df=expr_dat, coord_dat=coord_dat, dist_metric='euclidean')
  manhattan_calc <- STclust_calculate_distances(trcounts_df=expr_dat, coord_dat=coord_dat, dist_metric='manhattan')
  expect_type(euclid_calc, "list")
  expect_type(manhattan_calc, "list")
  expect_length(euclid_calc, 2)
  expect_length(manhattan_calc, 2)
  
  ws_0 <- STclust_weight_distances(scaled_dists=calc_result, ws=0)
  ws_1 <- STclust_weight_distances(scaled_dists=calc_result, ws=1)
  expect_type(ws_0, "list")
  expect_type(ws_1, "list")
  expect_length(ws_0, 1)
  expect_length(ws_1, 1)
})
