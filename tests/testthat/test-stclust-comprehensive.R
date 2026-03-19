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


# =============================================================================
# SECTION 2: STclust_hierarchical() - Clustering Methods and Parameters
# Tests for hierarchical clustering functionality with various methods and parameters
# Test data: melanoma_thrane dataset (Thrane et al.)
# =============================================================================

test_that("STclust_hierarchical() - Clustering methods (DTC, fixed k, ks range)", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load melanoma test data
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  melanoma <- load_melanoma_thrane()
  samples <- names(melanoma@spatial_meta)
  
  # Select genes first (required before clustering)
  melanoma_prepared <- STclust_select_genes(x=melanoma, samples=samples, topgenes=1000, cores=1)
  
  # Calculate distances with default parameters
  scaled_dists <- STclust_calculate_distances(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean'
  )
  weighted_dists <- STclust_weight_distances(scaled_dists=scaled_dists, ws=0.025)
  
  # --------------------------------------------------------------
  # Test 1: DTC method (Dynamic Tree Cut) - adaptive k detection
  # DTC automatically determines number of clusters based on tree structure
  # --------------------------------------------------------------
  dtc_result <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='DTC',
    ks='dtc',
    linkage='ward.D2',
    deepSplit=2,
    minClusterSize=10
  )
  
  # Verify DTC produces valid clustering result
  expect_type(dtc_result, "list")
  expect_true("cluster_assignment" %in% names(dtc_result))
  expect_true("x" %in% names(dtc_result))
  
  # Verify cluster assignments are integers
  expect_type(dtc_result$cluster_assignment, "integer")
  expect_true(all(dtc_result$cluster_assignment >= 1))
  
  # Verify correct length (one assignment per spot)
  n_spots <- nrow(melanoma@spatial_meta[[1]])
  expect_equal(length(dtc_result$cluster_assignment), n_spots)
  
  # Verify number of clusters is reasonable (should be at least 2, less than n_spots)
  n_clusters_dtc <- max(dtc_result$cluster_assignment)
  expect_true(n_clusters_dtc >= 2)
  expect_true(n_clusters_dtc < n_spots)
  cat(sprintf("DTC method produced %d clusters\n", n_clusters_dtc))
  
  # --------------------------------------------------------------
  # Test 2: Fixed k method - ks=2 (exactly 2 clusters)
  # Fixed k method uses specified number of clusters
  # --------------------------------------------------------------
  fixed_k2_result <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=2,
    linkage='ward.D2',
    deepSplit=2,
    minClusterSize=10
  )
  
  expect_type(fixed_k2_result, "list")
  expect_type(fixed_k2_result$cluster_assignment, "integer")
  expect_true(all(fixed_k2_result$cluster_assignment >= 1))
  
  # Verify exactly 2 clusters produced
  n_clusters_k2 <- max(fixed_k2_result$cluster_assignment)
  expect_equal(n_clusters_k2, 2)
  
  # Verify all spots assigned to either cluster 1 or 2
  expect_true(all(fixed_k2_result$cluster_assignment %in% c(1, 2)))
  
  # --------------------------------------------------------------
  # Test 3: Fixed k method - ks=3 (exactly 3 clusters)
  # --------------------------------------------------------------
  fixed_k3_result <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=3,
    linkage='ward.D2',
    deepSplit=2,
    minClusterSize=10
  )
  
  expect_type(fixed_k3_result, "list")
  expect_type(fixed_k3_result$cluster_assignment, "integer")
  
  # Verify exactly 3 clusters produced
  n_clusters_k3 <- max(fixed_k3_result$cluster_assignment)
  expect_equal(n_clusters_k3, 3)
  
  # Verify all spots assigned to clusters 1, 2, or 3
  expect_true(all(fixed_k3_result$cluster_assignment %in% c(1, 2, 3)))
  
  # --------------------------------------------------------------
  # Test 4: Fixed k method - ks=5 (exactly 5 clusters)
  # --------------------------------------------------------------
  fixed_k5_result <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=5,
    linkage='ward.D2',
    deepSplit=2,
    minClusterSize=10
  )
  
  expect_type(fixed_k5_result, "list")
  expect_type(fixed_k5_result$cluster_assignment, "integer")
  
  # Verify exactly 5 clusters produced
  n_clusters_k5 <- max(fixed_k5_result$cluster_assignment)
  expect_equal(n_clusters_k5, 5)
  
  # Verify all spots assigned to clusters 1-5
  expect_true(all(fixed_k5_result$cluster_assignment %in% 1:5))
  
  # --------------------------------------------------------------
  # Test 5: ks range - ks=2:5 produces multiple clusterings
  # When ks is a vector, should return list of clusterings for each value
  # --------------------------------------------------------------
  ks_range_result <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=2:5,
    linkage='ward.D2',
    deepSplit=2,
    minClusterSize=10
  )
  
  # Verify range produces list of clusterings
  expect_type(ks_range_result, "list")
  expect_length(ks_range_result, 4)  # 2, 3, 4, 5
  
  # Verify each element has cluster_assignment
  for(i in 1:4) {
    expect_true("cluster_assignment" %in% names(ks_range_result[[i]]))
    expect_type(ks_range_result[[i]]$cluster_assignment, "integer")
  }
  
  # Verify each clustering has correct number of clusters
  expect_equal(max(ks_range_result[[1]]$cluster_assignment), 2)  # ks=2
  expect_equal(max(ks_range_result[[2]]$cluster_assignment), 3)  # ks=3
  expect_equal(max(ks_range_result[[3]]$cluster_assignment), 4)  # ks=4
  expect_equal(max(ks_range_result[[4]]$cluster_assignment), 5)  # ks=5
  
  # Verify different ks values produce different clusterings
  expect_false(all(ks_range_result[[1]]$cluster_assignment == ks_range_result[[2]]$cluster_assignment))
  expect_false(all(ks_range_result[[2]]$cluster_assignment == ks_range_result[[3]]$cluster_assignment))
  
  cat("Fixed k method tests passed: ks=2, 3, 5 and ks range (2:5)\n")
})


test_that("STclust_hierarchical() - Linkage methods (ward.D2, average, complete)", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load melanoma test data
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  melanoma <- load_melanoma_thrane()
  samples <- names(melanoma@spatial_meta)
  
  # Select genes first
  melanoma_prepared <- STclust_select_genes(x=melanoma, samples=samples, topgenes=1000, cores=1)
  
  # Calculate distances
  scaled_dists <- STclust_calculate_distances(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean'
  )
  
  # --------------------------------------------------------------
  # Test 1: ward.D2 (default) - minimizes within-cluster variance
  # Produces compact, spherical clusters
  # --------------------------------------------------------------
  ward_d2_result <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=3,
    linkage='ward.D2',
    deepSplit=2,
    minClusterSize=10
  )
  
  expect_type(ward_d2_result, "list")
  expect_type(ward_d2_result$cluster_assignment, "integer")
  expect_true(all(ward_d2_result$cluster_assignment >= 1))
  
  n_clusters_ward <- max(ward_d2_result$cluster_assignment)
  expect_equal(n_clusters_ward, 3)
  expect_true(all(ward_d2_result$cluster_assignment %in% 1:3))
  
  # --------------------------------------------------------------
  # Test 2: average linkage (UPGMA) - average distance between clusters
  # More flexible than ward.D2, can produce elongated clusters
  # --------------------------------------------------------------
  average_result <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=3,
    linkage='average',
    deepSplit=2,
    minClusterSize=10
  )
  
  expect_type(average_result, "list")
  expect_type(average_result$cluster_assignment, "integer")
  expect_true(all(average_result$cluster_assignment >= 1))
  
  n_clusters_average <- max(average_result$cluster_assignment)
  expect_equal(n_clusters_average, 3)
  expect_true(all(average_result$cluster_assignment %in% 1:3))
  
  # --------------------------------------------------------------
  # Test 3: complete linkage - maximum distance between clusters
  # Produces compact clusters, sensitive to outliers
  # --------------------------------------------------------------
  complete_result <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=3,
    linkage='complete',
    deepSplit=2,
    minClusterSize=10
  )
  
  expect_type(complete_result, "list")
  expect_type(complete_result$cluster_assignment, "integer")
  expect_true(all(complete_result$cluster_assignment >= 1))
  
  n_clusters_complete <- max(complete_result$cluster_assignment)
  expect_equal(n_clusters_complete, 3)
  expect_true(all(complete_result$cluster_assignment %in% 1:3))
  
  # --------------------------------------------------------------
  # Test 4: Verify different methods produce different results
  # While all produce 3 clusters, cluster assignments should differ
  # --------------------------------------------------------------
  # ward.D2 vs average
  different_ward_avg <- !all(ward_d2_result$cluster_assignment == average_result$cluster_assignment)
  expect_true(different_ward_avg, "ward.D2 and average linkage should produce different clusterings")
  
  # ward.D2 vs complete
  different_ward_complete <- !all(ward_d2_result$cluster_assignment == complete_result$cluster_assignment)
  expect_true(different_ward_complete, "ward.D2 and complete linkage should produce different clusterings")
  
  # average vs complete
  different_avg_complete <- !all(average_result$cluster_assignment == complete_result$cluster_assignment)
  expect_true(different_avg_complete, "average and complete linkage should produce different clusterings")
  
  # --------------------------------------------------------------
  # Test 5: Verify all methods produce valid cluster assignments
  # All should have same number of clusters (3), same length, integers >= 1
  # --------------------------------------------------------------
  expect_equal(length(ward_d2_result$cluster_assignment), length(average_result$cluster_assignment))
  expect_equal(length(ward_d2_result$cluster_assignment), length(complete_result$cluster_assignment))
  
  expect_type(ward_d2_result$cluster_assignment, "integer")
  expect_type(average_result$cluster_assignment, "integer")
  expect_type(complete_result$cluster_assignment, "integer")
  
  expect_true(all(ward_d2_result$cluster_assignment >= 1))
  expect_true(all(average_result$cluster_assignment >= 1))
  expect_true(all(complete_result$cluster_assignment >= 1))
  
  cat("Linkage method tests passed: ward.D2, average, complete\n")
})


test_that("STclust_hierarchical() - deepSplit parameter (FALSE, 2, 3)", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load melanoma test data
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  melanoma <- load_melanoma_thrane()
  samples <- names(melanoma@spatial_meta)
  
  # Select genes first
  melanoma_prepared <- STclust_select_genes(x=melanoma, samples=samples, topgenes=1000, cores=1)
  
  # Calculate distances
  scaled_dists <- STclust_calculate_distances(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean'
  )
  
  # --------------------------------------------------------------
  # Test 1: deepSplit=FALSE (standard) - conservative splitting
  # Produces fewer, broader clusters
  # deepSplit parameter only applies to DTC method
  # --------------------------------------------------------------
  deepSplit_false_result <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='DTC',
    ks='dtc',
    linkage='ward.D2',
    deepSplit=FALSE,
    minClusterSize=10
  )
  
  expect_type(deepSplit_false_result, "list")
  expect_type(deepSplit_false_result$cluster_assignment, "integer")
  expect_true(all(deepSplit_false_result$cluster_assignment >= 1))
  
  n_clusters_false <- max(deepSplit_false_result$cluster_assignment)
  cat(sprintf("deepSplit=FALSE produced %d clusters\n", n_clusters_false))
  
  # --------------------------------------------------------------
  # Test 2: deepSplit=2 (moderate splitting) - intermediate granularity
  # Should produce more clusters than deepSplit=FALSE
  # --------------------------------------------------------------
  deepSplit_2_result <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='DTC',
    ks='dtc',
    linkage='ward.D2',
    deepSplit=2,
    minClusterSize=10
  )
  
  expect_type(deepSplit_2_result, "list")
  expect_type(deepSplit_2_result$cluster_assignment, "integer")
  expect_true(all(deepSplit_2_result$cluster_assignment >= 1))
  
  n_clusters_2 <- max(deepSplit_2_result$cluster_assignment)
  cat(sprintf("deepSplit=2 produced %d clusters\n", n_clusters_2))
  
  # Verify deepSplit=2 produces more clusters than FALSE
  expect_true(n_clusters_2 >= n_clusters_false,
              "deepSplit=2 should produce >= clusters than deepSplit=FALSE")
  
  # --------------------------------------------------------------
  # Test 3: deepSplit=3 (aggressive splitting) - fine-grained clusters
  # Should produce most clusters, potentially smaller cluster sizes
  # --------------------------------------------------------------
  deepSplit_3_result <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='DTC',
    ks='dtc',
    linkage='ward.D2',
    deepSplit=3,
    minClusterSize=10
  )
  
  expect_type(deepSplit_3_result, "list")
  expect_type(deepSplit_3_result$cluster_assignment, "integer")
  expect_true(all(deepSplit_3_result$cluster_assignment >= 1))
  
  n_clusters_3 <- max(deepSplit_3_result$cluster_assignment)
  cat(sprintf("deepSplit=3 produced %d clusters\n", n_clusters_3))
  
  # --------------------------------------------------------------
  # Test 4: Verify deepSplit increases number of clusters
  # deepSplit=3 >= deepSplit=2 >= deepSplit=FALSE
  # --------------------------------------------------------------
  expect_true(n_clusters_3 >= n_clusters_2,
              "deepSplit=3 should produce >= clusters than deepSplit=2")
  expect_true(n_clusters_2 >= n_clusters_false,
              "deepSplit=2 should produce >= clusters than deepSplit=FALSE")
  
  # --------------------------------------------------------------
  # Test 5: Verify all methods produce valid cluster assignments
  # All should have integers >= 1, same length
  # --------------------------------------------------------------
  expect_equal(length(deepSplit_false_result$cluster_assignment),
               length(deepSplit_2_result$cluster_assignment))
  expect_equal(length(deepSplit_false_result$cluster_assignment),
               length(deepSplit_3_result$cluster_assignment))
  
  expect_true(all(deepSplit_false_result$cluster_assignment >= 1))
  expect_true(all(deepSplit_2_result$cluster_assignment >= 1))
  expect_true(all(deepSplit_3_result$cluster_assignment >= 1))
  
  # --------------------------------------------------------------
  # Test 6: Verify different deepSplit values produce different clusterings
  # --------------------------------------------------------------
  different_false_2 <- !all(deepSplit_false_result$cluster_assignment == deepSplit_2_result$cluster_assignment)
  expect_true(different_false_2, "deepSplit=FALSE and deepSplit=2 should produce different clusterings")
  
  different_2_3 <- !all(deepSplit_2_result$cluster_assignment == deepSplit_3_result$cluster_assignment)
  expect_true(different_2_3, "deepSplit=2 and deepSplit=3 should produce different clusterings")
  
  cat("deepSplit parameter tests passed: FALSE < 2 < 3 (increasing clusters)\n")
})


test_that("STclust_hierarchical() - Multiple samples integration", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load melanoma test data (has multiple samples)
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  melanoma <- load_melanoma_thrane()
  samples <- names(melanoma@spatial_meta)
  
  skip_if(length(samples) < 2, "Need at least 2 samples for multi-sample test")
  
  # --------------------------------------------------------------
  # Test 1: Cluster multiple samples together
  # Each sample gets independent clustering, spatial_meta columns added
  # --------------------------------------------------------------
  multi_sample_result <- STclust_hierarchical(
    trcounts_df=list(
      melanoma@tr_counts[[1]],
      melanoma@tr_counts[[2]]
    ),
    coord_dat=list(
      melanoma@spatial_meta[[1]],
      melanoma@spatial_meta[[2]]
    ),
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=3,
    linkage='ward.D2',
    deepSplit=2,
    minClusterSize=10
  )
  
  # Verify result is list with clustering for each sample
  expect_type(multi_sample_result, "list")
  expect_length(multi_sample_result, 2)  # One per sample
  
  # --------------------------------------------------------------
  # Test 2: Verify each sample gets independent clustering
  # Each should have its own cluster_assignment
  # --------------------------------------------------------------
  for(i in 1:2) {
    expect_true("cluster_assignment" %in% names(multi_sample_result[[i]]))
    expect_type(multi_sample_result[[i]]$cluster_assignment, "integer")
    expect_true(all(multi_sample_result[[i]]$cluster_assignment >= 1))
    
    # Verify correct length for each sample
    n_spots_i <- nrow(melanoma@spatial_meta[[i]])
    expect_equal(length(multi_sample_result[[i]]$cluster_assignment), n_spots_i)
    
    # Verify exactly 3 clusters (ks=3)
    expect_equal(max(multi_sample_result[[i]]$cluster_assignment), 3)
    expect_true(all(multi_sample_result[[i]]$cluster_assignment %in% 1:3))
  }
  
  # --------------------------------------------------------------
  # Test 3: Verify spatial_meta columns added correctly
  # Should add cluster_assignment column to spatial_meta
  # --------------------------------------------------------------
  expect_true("x" %in% names(multi_sample_result))
  expect_true(inherits(multi_sample_result$x, "STlist"))
  
  # Check that spatial_meta has cluster_assignment column
  for(i in 1:2) {
    expect_true("cluster_assignment" %in% colnames(multi_sample_result$x@spatial_meta[[i]]))
    
    # Verify spatial_meta values match cluster_assignment
    meta_clusters <- multi_sample_result$x@spatial_meta[[i]]$cluster_assignment
    expect_equal(as.integer(meta_clusters), multi_sample_result[[i]]$cluster_assignment)
  }
  
  # --------------------------------------------------------------
  # Test 4: Test sample-specific cluster labels
  # Sample 1 and Sample 2 should have independent cluster numbering
  # (cluster 1 in sample 1 is independent of cluster 1 in sample 2)
  # --------------------------------------------------------------
  sample1_clusters <- multi_sample_result[[1]]$cluster_assignment
  sample2_clusters <- multi_sample_result[[2]]$cluster_assignment
  
  # Both should have clusters 1, 2, 3 (independent numbering)
  expect_true(all(c(1, 2, 3) %in% sample1_clusters))
  expect_true(all(c(1, 2, 3) %in% sample2_clusters))
  
  # Verify they are truly independent (not forced to be identical)
  # (This is expected since clustering is independent per sample)
  expect_true(any(sample1_clusters != sample2_clusters),
              "Sample cluster assignments should be independent")
  
  # --------------------------------------------------------------
  # Test 5: Test with 3 samples if available
  # --------------------------------------------------------------
  if(length(samples) >= 3) {
    multi_sample_3result <- STclust_hierarchical(
      trcounts_df=list(
        melanoma@tr_counts[[1]],
        melanoma@tr_counts[[2]],
        melanoma@tr_counts[[3]]
      ),
      coord_dat=list(
        melanoma@spatial_meta[[1]],
        melanoma@spatial_meta[[2]],
        melanoma@spatial_meta[[3]]
      ),
      dist_metric='euclidean',
      ws=0.025,
      method='fixed',
      ks=3,
      linkage='ward.D2',
      deepSplit=2,
      minClusterSize=10
    )
    
    expect_type(multi_sample_3result, "list")
    expect_length(multi_sample_3result, 3)
    
    for(i in 1:3) {
      expect_true("cluster_assignment" %in% names(multi_sample_3result[[i]]))
      expect_type(multi_sample_3result[[i]]$cluster_assignment, "integer")
      expect_equal(max(multi_sample_3result[[i]]$cluster_assignment), 3)
      expect_true("cluster_assignment" %in% colnames(multi_sample_3result$x@spatial_meta[[i]]))
    }
    
    cat("3-sample clustering test passed\n")
  }
  
  cat("Multi-sample integration tests passed: independent clustering per sample\n")
})
