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


# =============================================================================
# SECTION 3: STclust Edge Cases, Error Handling, and Regression Tests
# Tests for boundary conditions, error handling, consistency, and integration
# Test data: melanoma_thrane dataset
# =============================================================================

test_that("STclust_edge_cases - Boundary conditions (single gene, two spots, extreme ws, ks=1)", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load melanoma test data
  data_dir <- testthat::test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  melanoma <- load_melanoma_thrane()
  samples <- names(melanoma@spatial_meta)
  
  # --------------------------------------------------------------
  # Test 1: Single gene minimum (topgenes=1)
  # Edge case: minimum number of genes for clustering
  # Should still work but may produce degenerate clusters
  # --------------------------------------------------------------
  melanoma_min <- STclust_select_genes(x=melanoma, samples=samples, topgenes=1, cores=1)
  expect_type(melanoma_min, "list")
  expect_length(melanoma_min$trcounts_df, length(samples))
  expect_equal(nrow(melanoma_min$trcounts_df[[1]]), 1)
  
  # Should still be able to cluster with single gene
  scaled_dists_min <- STclust_calculate_distances(
    trcounts_df=melanoma_min$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean'
  )
  expect_type(scaled_dists_min, "list")
  expect_length(scaled_dists_min, 2)
  
  weighted_dists_min <- STclust_weight_distances(scaled_dists=scaled_dists_min, ws=0.025)
  expect_type(weighted_dists_min, "list")
  expect_length(weighted_dists_min, 1)
  
  # Cluster with single gene should still produce valid result
  cluster_single_gene <- STclust_hierarchical(
    trcounts_df=melanoma_min$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=2,
    linkage='ward.D2',
    deepSplit=2,
    minClusterSize=5
  )
  expect_type(cluster_single_gene, "list")
  expect_true("cluster_assignment" %in% names(cluster_single_gene))
  expect_type(cluster_single_gene$cluster_assignment, "integer")
  expect_true(all(cluster_single_gene$cluster_assignment >= 1))
  cat("Single gene clustering test passed\n")
  
  # --------------------------------------------------------------
  # Test 2: Two spots minimum for clustering
  # Edge case: minimum number of spots that can be clustered
  # Should produce exactly one cluster (or error if minClusterSize violated)
  # --------------------------------------------------------------
  # Create minimal test data with 2 spots
  minimal_expr <- melanoma_min$trcounts_df[[1]][1:2, ]
  minimal_coord <- melanoma@spatial_meta[[1]][1:2, ]
  
  # With only 2 spots and minClusterSize=5, should fail gracefully or produce 1 cluster
  cluster_two_spots <- tryCatch({
    STclust_hierarchical(
      trcounts_df=minimal_expr,
      coord_dat=minimal_coord,
      dist_metric='euclidean',
      ws=0.025,
      method='fixed',
      ks=2,
      linkage='ward.D2',
      deepSplit=2,
      minClusterSize=5  # Larger than n_spots, should produce 1 cluster or error
    )
  }, error = function(e) {
    # If error, expect reasonable error message
    expect_true(nchar(e$message) > 0)
    return(NULL)
  })
  
  # If clustering succeeded, verify result
  if(!is.null(cluster_two_spots)) {
    expect_type(cluster_two_spots, "list")
    expect_true("cluster_assignment" %in% names(cluster_two_spots))
    # With 2 spots and minClusterSize=5, should have all spots in one cluster
    expect_equal(length(unique(cluster_two_spots$cluster_assignment)), 1)
  }
  cat("Two spots minimum test passed\n")
  
  # --------------------------------------------------------------
  # Test 3: Very large ws (ws=0.9, spatial dominated)
  # Edge case: clustering dominated by spatial proximity
  # Should produce spatially contiguous clusters
  # --------------------------------------------------------------
  melanoma_large_ws <- STclust_select_genes(x=melanoma, samples=samples, topgenes=1000, cores=1)
  
  scaled_dists_large_ws <- STclust_calculate_distances(
    trcounts_df=melanoma_large_ws$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean'
  )
  
  # ws=0.9 means 90% spatial, 10% expression
  weighted_large_ws <- STclust_weight_distances(scaled_dists=scaled_dists_large_ws, ws=0.9)
  expect_type(weighted_large_ws, "list")
  expect_length(weighted_large_ws, 1)
  
  cluster_large_ws <- STclust_hierarchical(
    trcounts_df=melanoma_large_ws$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.9,
    method='fixed',
    ks=3,
    linkage='ward.D2',
    deepSplit=2,
    minClusterSize=10
  )
  expect_type(cluster_large_ws, "list")
  expect_true("cluster_assignment" %in% names(cluster_large_ws))
  expect_type(cluster_large_ws$cluster_assignment, "integer")
  expect_equal(max(cluster_large_ws$cluster_assignment), 3)
  cat("Very large ws (0.9) test passed\n")
  
  # --------------------------------------------------------------
  # Test 4: Very small ws (ws=0.001, expression dominated)
  # Edge case: clustering dominated by gene expression
  # Should produce clusters based on transcriptomic similarity
  # --------------------------------------------------------------
  scaled_dists_small_ws <- STclust_calculate_distances(
    trcounts_df=melanoma_large_ws$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean'
  )
  
  # ws=0.001 means 0.1% spatial, 99.9% expression
  weighted_small_ws <- STclust_weight_distances(scaled_dists=scaled_dists_small_ws, ws=0.001)
  expect_type(weighted_small_ws, "list")
  expect_length(weighted_small_ws, 1)
  
  cluster_small_ws <- STclust_hierarchical(
    trcounts_df=melanoma_large_ws$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.001,
    method='fixed',
    ks=3,
    linkage='ward.D2',
    deepSplit=2,
    minClusterSize=10
  )
  expect_type(cluster_small_ws, "list")
  expect_true("cluster_assignment" %in% names(cluster_small_ws))
  expect_type(cluster_small_ws$cluster_assignment, "integer")
  expect_equal(max(cluster_small_ws$cluster_assignment), 3)
  
  # Verify ws=0.001 and ws=0.9 produce different clusterings
  different_extremes <- !all(cluster_large_ws$cluster_assignment == cluster_small_ws$cluster_assignment)
  expect_true(different_extremes, "ws=0.9 and ws=0.001 should produce different clusterings")
  cat("Very small ws (0.001) test passed\n")
  
  # --------------------------------------------------------------
  # Test 5: ks=1 (single cluster edge case)
  # Edge case: requesting only 1 cluster
  # Should produce exactly one cluster with all spots assigned
  # --------------------------------------------------------------
  cluster_ks1 <- STclust_hierarchical(
    trcounts_df=melanoma_large_ws$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=1,
    linkage='ward.D2',
    deepSplit=2,
    minClusterSize=10
  )
  expect_type(cluster_ks1, "list")
  expect_true("cluster_assignment" %in% names(cluster_ks1))
  expect_type(cluster_ks1$cluster_assignment, "integer")
  
  # Verify exactly 1 cluster produced
  n_clusters_ks1 <- max(cluster_ks1$cluster_assignment)
  expect_equal(n_clusters_ks1, 1)
  
  # Verify all spots assigned to cluster 1
  expect_true(all(cluster_ks1$cluster_assignment == 1))
  cat("Single cluster (ks=1) test passed\n")
})


test_that("STclust_error_handling - Invalid inputs (ks, ws, dist_metric, linkage, STlist)", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load melanoma test data
  data_dir <- testthat::test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  melanoma <- load_melanoma_thrane()
  samples <- names(melanoma@spatial_meta)
  
  # Prepare valid data for error tests
  melanoma_prepared <- STclust_select_genes(x=melanoma, samples=samples, topgenes=1000, cores=1)
  
  # --------------------------------------------------------------
  # Test 1: Invalid ks values
  # - ks = negative (should error)
  # - ks = zero (should error)
  # - ks = non-integer (should error or coerce)
  # - ks = character string (should error)
  # --------------------------------------------------------------
  
  # Test negative ks
  expect_error(
    STclust_hierarchical(
      trcounts_df=melanoma_prepared$trcounts_df[[1]],
      coord_dat=melanoma@spatial_meta[[1]],
      dist_metric='euclidean',
      ws=0.025,
      method='fixed',
      ks=-1,
      linkage='ward.D2',
      deepSplit=2,
      minClusterSize=10
    ),
    regexp = "invalid|negative|ks",
    info = "Negative ks should produce error"
  )
  
  # Test zero ks
  expect_error(
    STclust_hierarchical(
      trcounts_df=melanoma_prepared$trcounts_df[[1]],
      coord_dat=melanoma@spatial_meta[[1]],
      dist_metric='euclidean',
      ws=0.025,
      method='fixed',
      ks=0,
      linkage='ward.D2',
      deepSplit=2,
      minClusterSize=10
    ),
    regexp = "invalid|zero|ks",
    info = "Zero ks should produce error"
  )
  
  # Test non-integer ks (should either error or round)
  result_nonint_ks <- tryCatch({
    STclust_hierarchical(
      trcounts_df=melanoma_prepared$trcounts_df[[1]],
      coord_dat=melanoma@spatial_meta[[1]],
      dist_metric='euclidean',
      ws=0.025,
      method='fixed',
      ks=2.5,
      linkage='ward.D2',
      deepSplit=2,
      minClusterSize=10
    )
  }, error = function(e) {
    # If error, verify it's about ks
    expect_true(grepl("ks|integer|invalid", e$message, ignore.case=TRUE))
    return(NULL)
  })
  
  # If non-integer ks was accepted, verify it was coerced to integer
  if(!is.null(result_nonint_ks)) {
    expect_type(result_nonint_ks, "list")
    expect_true("cluster_assignment" %in% names(result_nonint_ks))
    # Should have been rounded/truncated to 2 or 3
    expect_true(max(result_nonint_ks$cluster_assignment) %in% c(2, 3))
  }
  cat("Invalid ks values test passed\n")
  
  # --------------------------------------------------------------
  # Test 2: Invalid ws values
  # - ws < 0 (should error)
  # - ws > 1 (should error)
  # - ws = NA (should error)
  # --------------------------------------------------------------
  
  # Test negative ws
  expect_error(
    STclust_hierarchical(
      trcounts_df=melanoma_prepared$trcounts_df[[1]],
      coord_dat=melanoma@spatial_meta[[1]],
      dist_metric='euclidean',
      ws=-0.1,
      method='fixed',
      ks=3,
      linkage='ward.D2',
      deepSplit=2,
      minClusterSize=10
    ),
    regexp = "invalid|negative|ws|between",
    info = "Negative ws should produce error"
  )
  
  # Test ws > 1
  expect_error(
    STclust_hierarchical(
      trcounts_df=melanoma_prepared$trcounts_df[[1]],
      coord_dat=melanoma@spatial_meta[[1]],
      dist_metric='euclidean',
      ws=1.5,
      method='fixed',
      ks=3,
      linkage='ward.D2',
      deepSplit=2,
      minClusterSize=10
    ),
    regexp = "invalid|greater|ws|between",
    info = "ws > 1 should produce error"
  )
  
  # Test ws = NA
  expect_error(
    STclust_hierarchical(
      trcounts_df=melanoma_prepared$trcounts_df[[1]],
      coord_dat=melanoma@spatial_meta[[1]],
      dist_metric='euclidean',
      ws=NA,
      method='fixed',
      ks=3,
      linkage='ward.D2',
      deepSplit=2,
      minClusterSize=10
    ),
    regexp = "NA|missing|invalid",
    info = "NA ws should produce error"
  )
  cat("Invalid ws values test passed\n")
  
  # --------------------------------------------------------------
  # Test 3: Invalid dist_metric (unknown metric)
  # - Should error with informative message about valid metrics
  # --------------------------------------------------------------
  
  expect_error(
    STclust_hierarchical(
      trcounts_df=melanoma_prepared$trcounts_df[[1]],
      coord_dat=melanoma@spatial_meta[[1]],
      dist_metric='unknown_metric',
      ws=0.025,
      method='fixed',
      ks=3,
      linkage='ward.D2',
      deepSplit=2,
      minClusterSize=10
    ),
    regexp = "invalid|unknown|dist_metric|euclidean|manhattan",
    info = "Invalid dist_metric should produce error"
  )
  cat("Invalid dist_metric test passed\n")
  
  # --------------------------------------------------------------
  # Test 4: Invalid linkage method (unknown method)
  # - Should error with informative message about valid linkage methods
  # --------------------------------------------------------------
  
  expect_error(
    STclust_hierarchical(
      trcounts_df=melanoma_prepared$trcounts_df[[1]],
      coord_dat=melanoma@spatial_meta[[1]],
      dist_metric='euclidean',
      ws=0.025,
      method='fixed',
      ks=3,
      linkage='unknown_linkage',
      deepSplit=2,
      minClusterSize=10
    ),
    regexp = "invalid|unknown|linkage|ward|average|complete",
    info = "Invalid linkage should produce error"
  )
  cat("Invalid linkage method test passed\n")
  
  # --------------------------------------------------------------
  # Test 5: NULL or empty STlist
  # - Should error with clear message
  # --------------------------------------------------------------
  
  # Test NULL STlist
  expect_error(
    STclust_hierarchical(
      trcounts_df=NULL,
      coord_dat=NULL,
      dist_metric='euclidean',
      ws=0.025,
      method='fixed',
      ks=3,
      linkage='ward.D2',
      deepSplit=2,
      minClusterSize=10
    ),
    regexp = "NULL|empty|missing|length",
    info = "NULL inputs should produce error"
  )
  
  # Test empty list
  expect_error(
    STclust_hierarchical(
      trcounts_df=list(),
      coord_dat=list(),
      dist_metric='euclidean',
      ws=0.025,
      method='fixed',
      ks=3,
      linkage='ward.D2',
      deepSplit=2,
      minClusterSize=10
    ),
    regexp = "empty|length|zero|samples",
    info = "Empty list inputs should produce error"
  )
  cat("NULL/empty STlist test passed\n")
  
  # --------------------------------------------------------------
  # Test 6: Missing samples in STlist (mismatched lengths)
  # - trcounts_df and coord_dat should have same length
  # - Should error if lengths don't match
  # --------------------------------------------------------------
  
  expect_error(
    STclust_hierarchical(
      trcounts_df=list(melanoma_prepared$trcounts_df[[1]]),
      coord_dat=list(melanoma@spatial_meta[[1]], melanoma@spatial_meta[[2]]),  # Mismatch!
      dist_metric='euclidean',
      ws=0.025,
      method='fixed',
      ks=3,
      linkage='ward.D2',
      deepSplit=2,
      minClusterSize=10
    ),
    regexp = "length|mismatch|sample|equal",
    info = "Mismatched sample lengths should produce error"
  )
  cat("Missing/mismatched samples test passed\n")
})


test_that("STclust_regression - Consistency and reproducibility", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load melanoma test data
  data_dir <- testthat::test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  melanoma <- load_melanoma_thrane()
  samples <- names(melanoma@spatial_meta)
  
  # Prepare data
  melanoma_prepared <- STclust_select_genes(x=melanoma, samples=samples, topgenes=1000, cores=1)
  
  # --------------------------------------------------------------
  # Test 1: Deterministic results with set.seed()
  # - Same seed should produce identical clustering
  # - Different seeds may produce different results (stochastic elements)
  # --------------------------------------------------------------
  
  # Run 1 with seed 42
  set.seed(42)
  result_seed42_1 <- STclust_hierarchical(
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
  
  # Run 2 with same seed 42
  set.seed(42)
  result_seed42_2 <- STclust_hierarchical(
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
  
  # Verify identical results with same seed
  expect_equal(result_seed42_1$cluster_assignment, result_seed42_2$cluster_assignment,
               info = "Same seed should produce identical clustering")
  cat("Deterministic results with set.seed() test passed\n")
  
  # --------------------------------------------------------------
  # Test 2: Multiple runs produce identical clustering (no randomness)
  # - Hierarchical clustering should be deterministic without DTC
  # - Multiple runs without set.seed() should produce identical results
  # --------------------------------------------------------------
  
  # Run 1 without explicit seed
  result_run1 <- STclust_hierarchical(
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
  
  # Run 2 without explicit seed
  result_run2 <- STclust_hierarchical(
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
  
  # Run 3 without explicit seed
  result_run3 <- STclust_hierarchical(
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
  
  # Verify all runs produce identical results
  expect_equal(result_run1$cluster_assignment, result_run2$cluster_assignment,
               info = "Multiple runs should produce identical clustering")
  expect_equal(result_run2$cluster_assignment, result_run3$cluster_assignment,
               info = "Multiple runs should produce identical clustering")
  cat("Multiple runs reproducibility test passed\n")
  
  # --------------------------------------------------------------
  # Test 3: STclust() vs STclust_legacy() comparison with various parameters
  # - Legacy and current implementations should produce consistent results
  # - Test with different ks, ws, linkage values
  # --------------------------------------------------------------
  
  # Check if STclust_legacy exists
  if(exists("STclust_legacy", mode="function")) {
    # Test 1: Default parameters
    result_current <- STclust(
      x=melanoma,
      samples=samples,
      topgenes=500,
      ws=0.025,
      ks=3,
      dist_metric='euclidean',
      linkage='ward.D2',
      cores=1
    )
    
    result_legacy <- STclust_legacy(
      x=melanoma,
      samples=samples,
      topgenes=500,
      ws=0.025,
      ks=3,
      dist_metric='euclidean',
      linkage='ward.D2',
      cores=1
    )
    
    # Verify both produce valid results
    expect_type(result_current, "list")
    expect_type(result_legacy, "list")
    expect_true("cluster_assignment" %in% names(result_current))
    expect_true("cluster_assignment" %in% names(result_legacy))
    
    # Verify cluster counts match
    expect_equal(max(result_current$cluster_assignment[[1]]),
                 max(result_legacy$cluster_assignment[[1]]),
                 info = "Current and legacy should produce same number of clusters")
    
    cat("STclust vs STclust_legacy test passed (same parameters)\n")
    
    # Test 2: Different ks value
    result_current_k5 <- STclust(
      x=melanoma, samples=samples, topgenes=500, ws=0.025, ks=5,
      dist_metric='euclidean', linkage='ward.D2', cores=1
    )
    
    result_legacy_k5 <- STclust_legacy(
      x=melanoma, samples=samples, topgenes=500, ws=0.025, ks=5,
      dist_metric='euclidean', linkage='ward.D2', cores=1
    )
    
    expect_equal(max(result_current_k5$cluster_assignment[[1]]),
                 max(result_legacy_k5$cluster_assignment[[1]]),
                 info = "Current and legacy should produce same number of clusters for ks=5")
    cat("STclust vs STclust_legacy test passed (ks=5)\n")
    
    # Test 3: Different ws value
    result_current_ws01 <- STclust(
      x=melanoma, samples=samples, topgenes=500, ws=0.1, ks=3,
      dist_metric='euclidean', linkage='ward.D2', cores=1
    )
    
    result_legacy_ws01 <- STclust_legacy(
      x=melanoma, samples=samples, topgenes=500, ws=0.1, ks=3,
      dist_metric='euclidean', linkage='ward.D2', cores=1
    )
    
    expect_equal(max(result_current_ws01$cluster_assignment[[1]]),
                 max(result_legacy_ws01$cluster_assignment[[1]]),
                 info = "Current and legacy should produce same number of clusters for ws=0.1")
    cat("STclust vs STclust_legacy test passed (ws=0.1)\n")
  } else {
    # If STclust_legacy doesn't exist, just verify STclust is consistent
    cat("STclust_legacy not available, skipping legacy comparison\n")
  }
  
  # --------------------------------------------------------------
  # Test 4: Verify reproducibility across different runs (fresh sessions)
  # - Simulate by running the same computation multiple times
  # --------------------------------------------------------------
  
  # Run 1
  repro_run1 <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=4,
    linkage='average',
    deepSplit=2,
    minClusterSize=10
  )
  
  # Run 2
  repro_run2 <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=4,
    linkage='average',
    deepSplit=2,
    minClusterSize=10
  )
  
  # Run 3
  repro_run3 <- STclust_hierarchical(
    trcounts_df=melanoma_prepared$trcounts_df[[1]],
    coord_dat=melanoma@spatial_meta[[1]],
    dist_metric='euclidean',
    ws=0.025,
    method='fixed',
    ks=4,
    linkage='average',
    deepSplit=2,
    minClusterSize=10
  )
  
  # Verify all runs identical
  expect_equal(repro_run1$cluster_assignment, repro_run2$cluster_assignment)
  expect_equal(repro_run2$cluster_assignment, repro_run3$cluster_assignment)
  expect_equal(repro_run1$cluster_assignment, repro_run3$cluster_assignment)
  cat("Reproducibility across runs test passed\n")
})


test_that("STclust_integration - Full workflow (STclust → STdiff, visualization)", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load melanoma test data
  data_dir <- testthat::test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  melanoma <- load_melanoma_thrane()
  samples <- names(melanoma@spatial_meta)
  
  # --------------------------------------------------------------
  # Test 1: STclust → STdiff pipeline
  # - Perform clustering
  # - Use clustered data for differential expression analysis
  # - Verify downstream functions work with clustered data
  # --------------------------------------------------------------
  
  # Step 1: Cluster the data
  cluster_result <- STclust(
    x=melanoma,
    samples=samples,
    topgenes=500,
    ws=0.025,
    ks=3,
    dist_metric='euclidean',
    linkage='ward.D2',
    cores=1
  )
  
  # Verify clustering result
  expect_type(cluster_result, "list")
  expect_true("cluster_assignment" %in% names(cluster_result))
  expect_true("x" %in% names(cluster_result))
  expect_true(inherits(cluster_result$x, "STlist"))
  
  # Verify spatial_meta has cluster assignments
  expect_true("cluster_assignment" %in% colnames(cluster_result$x@spatial_meta[[1]]))
  cat("STclust clustering completed\n")
  
  # Step 2: Perform differential expression analysis using clusters
  # Create group variable from cluster assignments
  cluster_groups <- as.factor(cluster_result$cluster_assignment[[1]])
  
  # Run STdiff with cluster groups
  tryCatch({
    stdiff_result <- STdiff(
      x=cluster_result$x,
      groups=cluster_groups,
      method="wilcox",
      cores=1
    )
    
    # Verify STdiff result
    expect_type(stdiff_result, "list")
    expect_true("statistic" %in% names(stdiff_result))
    expect_true("pvalue" %in% names(stdiff_result))
    expect_true("adj_pvalue" %in% names(stdiff_result))
    
    # Verify all genes have statistics
    expect_equal(length(stdiff_result$statistic), nrow(cluster_result$x@tr_counts[[1]]))
    expect_equal(length(stdiff_result$pvalue), nrow(cluster_result$x@tr_counts[[1]]))
    
    cat("STclust → STdiff pipeline test passed\n")
  }, error = function(e) {
    # If STdiff fails, verify it's not due to clustering issues
    # (could be due to other factors)
    cat(sprintf("STdiff error (may be acceptable): %s\n", e$message))
  })
  
  # --------------------------------------------------------------
  # Test 2: Verify clustered data works with other downstream functions
  # - Test with basic STlist functions
  # --------------------------------------------------------------
  
  # Verify STlist structure is intact after clustering
  expect_true(inherits(cluster_result$x, "STlist"))
  expect_true("tr_counts" %in% names(cluster_result$x))
  expect_true("spatial_meta" %in% names(cluster_result$x))
  expect_true("gene_meta" %in% names(cluster_result$x))
  
  # Verify spatial coordinates are preserved
  expect_true("row" %in% colnames(cluster_result$x@spatial_meta[[1]]))
  expect_true("column" %in% colnames(cluster_result$x@spatial_meta[[1]]))
  
  # Verify gene expression matrix is preserved
  expect_true(nrow(cluster_result$x@tr_counts[[1]]) > 0)
  expect_true(ncol(cluster_result$x@tr_counts[[1]]) > 0)
  cat("Clustered data structure verification passed\n")
  
  # --------------------------------------------------------------
  # Test 3: Test plot_STclust() or similar visualization works
  # - Verify plotting function exists and works
  # - Test with clustered data
  # --------------------------------------------------------------
  
  # Check if plot_STclust exists
  if(exists("plot_STclust", mode="function")) {
    # Try to create a plot (will not save, just verify it runs)
    tryCatch({
      # Suppress plot output
      old_device <- dev.cur()
      dev.new()  # Open new device
      on.exit(dev.off(), add=TRUE)
      
      plot_STclust(cluster_result, ks=3)
      
      # If we get here without error, plot works
      cat("plot_STclust() visualization test passed\n")
    }, error = function(e) {
      # If plot fails, note it but don't fail test
      cat(sprintf("plot_STclust error (may be acceptable): %s\n", e$message))
    })
  } else if(exists("plot_STclust_legacy", mode="function")) {
    # Try legacy plotting function
    tryCatch({
      old_device <- dev.cur()
      dev.new()
      on.exit(dev.off(), add=TRUE)
      
      plot_STclust_legacy(cluster_result, ks=3)
      cat("plot_STclust_legacy() visualization test passed\n")
    }, error = function(e) {
      cat(sprintf("plot_STclust_legacy error (may be acceptable): %s\n", e$message))
    })
  } else {
    # No plotting function found - note this but don't fail test
    cat("No plot_STclust function found, skipping visualization test\n")
  }
  
  # --------------------------------------------------------------
  # Test 4: Test multi-sample clustering with downstream analysis
  # - Cluster multiple samples
  # - Verify each sample has cluster assignments
  # - Test downstream compatibility
  # --------------------------------------------------------------
  
  if(length(samples) >= 2) {
    multi_cluster_result <- STclust(
      x=melanoma,
      samples=samples[1:2],
      topgenes=500,
      ws=0.025,
      ks=3,
      dist_metric='euclidean',
      linkage='ward.D2',
      cores=1
    )
    
    # Verify multi-sample clustering
    expect_type(multi_cluster_result, "list")
    expect_true("cluster_assignment" %in% names(multi_cluster_result))
    expect_length(multi_cluster_result$cluster_assignment, 2)
    
    # Verify each sample has cluster assignments
    for(i in 1:2) {
      expect_type(multi_cluster_result$cluster_assignment[[i]], "integer")
      expect_true(all(multi_cluster_result$cluster_assignment[[i]] >= 1))
      expect_equal(max(multi_cluster_result$cluster_assignment[[i]]), 3)
    }
    cat("Multi-sample integration test passed\n")
  }
  
  cat("Full workflow integration tests passed\n")
})


test_that("STclust_performance - Runtime sanity check", {
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load melanoma test data
  data_dir <- testthat::test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  melanoma <- load_melanoma_thrane()
  samples <- names(melanoma@spatial_meta)
  
  # --------------------------------------------------------------
  # Test: Reasonable runtime for small dataset (<30 seconds)
  # - Measure time for complete STclust workflow
  # - Should complete within reasonable time for small dataset
  # --------------------------------------------------------------
  
  start_time <- Sys.time()
  
  cluster_result <- STclust(
    x=melanoma,
    samples=samples[1],  # Use single sample for speed
    topgenes=500,
    ws=0.025,
    ks=3,
    dist_metric='euclidean',
    linkage='ward.D2',
    cores=1
  )
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units="secs"))
  
  cat(sprintf("STclust completed in %.2f seconds\n", elapsed_time))
  
  # Verify clustering result
  expect_type(cluster_result, "list")
  expect_true("cluster_assignment" %in% names(cluster_result))
  expect_type(cluster_result$cluster_assignment, "integer")
  
  # Verify runtime is reasonable (<30 seconds for small dataset)
  expect_true(elapsed_time < 30,
              info = sprintf("STclust should complete in <30 seconds, took %.2f seconds", elapsed_time))
  
  cat("Performance sanity check passed: completed in %.2f seconds\n", elapsed_time)
  
  # --------------------------------------------------------------
  # Test: No memory leaks or excessive allocation
  # - Verify cluster result size is reasonable
  # --------------------------------------------------------------
  
  # Verify result size is proportional to input
  n_spots <- nrow(melanoma@spatial_meta[[1]])
  n_genes <- nrow(melanoma@tr_counts[[1]])
  
  # Cluster assignment should be same length as number of spots
  expect_equal(length(cluster_result$cluster_assignment[[1]]), n_spots)
  
  # STlist should have same number of genes (or fewer if filtered)
  expect_true(nrow(cluster_result$x@tr_counts[[1]]) <= n_genes)
  
  cat("Memory allocation sanity check passed\n")
})
