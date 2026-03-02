##
# STclust Tests
# Tests for the refactored STclust modular implementation

test_that("STclust basic clustering works", {
  # Unload installed package to avoid conflicts
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  # Load local development version
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load data
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  
  count_files <- list.files(data_dir, full.names=TRUE, pattern='counts')
  coord_files <- list.files(data_dir, full.names=TRUE, pattern='mapping')
  clin_file <- list.files(data_dir, full.names=TRUE, pattern='clinical')
  
  melanoma <- STlist(rnacounts = count_files[1], spotcoords = coord_files[1], samples = clin_file)
  melanoma <- transform_data(melanoma)
  
  # Test basic clustering with ks=2:3, ws=0.025
  res <- STclust(melanoma, ks=2:3, ws=0.025, verbose=0)
  
  expect_true(inherits(res, "STlist"))
  
  # Check that clustering columns were added
  expect_true(any(grepl("stclust_spw0.025_k", colnames(res@spatial_meta[[1]]))))
})


test_that("STclust DTC clustering works", {
  # Unload installed package to avoid conflicts
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  # Load local development version
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load data
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  
  count_files <- list.files(data_dir, full.names=TRUE, pattern='counts')
  coord_files <- list.files(data_dir, full.names=TRUE, pattern='mapping')
  clin_file <- list.files(data_dir, full.names=TRUE, pattern='clinical')
  
  melanoma <- STlist(rnacounts = count_files[1], spotcoords = coord_files[1], samples = clin_file)
  melanoma <- transform_data(melanoma)
  
  # Test DTC clustering (adaptive k)
  res <- STclust(melanoma, ks='dtc', ws=0.025, deepSplit=2, verbose=0)
  
  expect_true(inherits(res, "STlist"))
  
  # Check that DTC clustering column was added
  expect_true(any(grepl("stclust_spw0.025_dspl", colnames(res@spatial_meta[[1]]))))
})


test_that("STclust fixed k clustering works", {
  # Unload installed package to avoid conflicts
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  # Load local development version
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load data
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  
  count_files <- list.files(data_dir, full.names=TRUE, pattern='counts')
  coord_files <- list.files(data_dir, full.names=TRUE, pattern='mapping')
  clin_file <- list.files(data_dir, full.names=TRUE, pattern='clinical')
  
  melanoma <- STlist(rnacounts = count_files[1], spotcoords = coord_files[1], samples = clin_file)
  melanoma <- transform_data(melanoma)
  
  # Test fixed k clustering
  res <- STclust(melanoma, ks=2, ws=0.025, verbose=0)
  
  expect_true(inherits(res, "STlist"))
  
  # Check that fixed k clustering column was added
  expect_true(any(grepl("stclust_spw0.025_k2", colnames(res@spatial_meta[[1]]))))
})


test_that("STclust() produces identical results to STclust_legacy()", {
  # Unload installed package to avoid conflicts
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  # Load local development version
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load data
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  
  count_files <- list.files(data_dir, full.names=TRUE, pattern='counts')
  coord_files <- list.files(data_dir, full.names=TRUE, pattern='mapping')
  clin_file <- list.files(data_dir, full.names=TRUE, pattern='clinical')
  
  melanoma <- STlist(rnacounts = count_files[1], spotcoords = coord_files[1], samples = clin_file)
  melanoma <- transform_data(melanoma)
  
  # Run both versions with same parameters
  res_new <- STclust(melanoma, ks=2:3, ws=0.025, verbose=0)
  res_legacy <- STclust_legacy(melanoma, ks=2:3, ws=0.025, verbose=0)
  
  # Compare results
  
  # Check that both have same samples
  expect_identical(names(res_new@spatial_meta), names(res_legacy@spatial_meta))
  
  # Check that both have same clustering columns
  new_cols = grep("stclust", colnames(res_new@spatial_meta[[1]]), value=TRUE)
  legacy_cols = grep("stclust", colnames(res_legacy@spatial_meta[[1]]), value=TRUE)
  expect_identical(sort(new_cols), sort(legacy_cols))
  
  # Compare cluster assignments for each column
  for(col_name in new_cols){
    expect_identical(
      res_new@spatial_meta[[1]][[col_name]],
      res_legacy@spatial_meta[[1]][[col_name]]
    )
  }
})


test_that("STclust helper functions work correctly", {
  # Unload installed package to avoid conflicts
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  # Load local development version
  devtools::load_all('../../.', export_all=TRUE)
  
  # Load data
  data_dir <- test_path("data/melanoma_thrane")
  skip_if_not(dir.exists(data_dir), "Melanoma data not found")
  
  count_files <- list.files(data_dir, full.names=TRUE, pattern='counts')
  coord_files <- list.files(data_dir, full.names=TRUE, pattern='mapping')
  clin_file <- list.files(data_dir, full.names=TRUE, pattern='clinical')
  
  melanoma <- STlist(rnacounts = count_files[1], spotcoords = coord_files[1], samples = clin_file)
  melanoma <- transform_data(melanoma)
  
  # Test calculate_dist_matrices
  # Need to filter to ensure enough spots for clustering
  # Use all genes to ensure enough spots
  expr_dat = melanoma@tr_counts[[1]]
  coord_dat = melanoma@spatial_meta[[1]]
  dists = calculate_dist_matrices(expr_dat=expr_dat, coord_dat=coord_dat, dist_metric='euclidean')
  
  expect_type(dists, "list")
  expect_true(all(c("scale_exp", "scale_coord") %in% names(dists)))
  expect_true(all(dim(dists$scale_exp) == c(nrow(coord_dat), nrow(coord_dat))))
  expect_true(all(dim(dists$scale_coord) == c(nrow(coord_dat), nrow(coord_dat))))
  # Check scaled to [0,1]
  expect_true(all(dists$scale_exp >= 0 & dists$scale_exp <= 1))
  expect_true(all(dists$scale_coord >= 0 & dists$scale_coord <= 1))
  
  # Test calculate_weighted_dist
  ws = c(0.025, 0.05)
  weighted = calculate_weighted_dist(scaled_dists=dists, ws=ws)
  
  expect_type(weighted, "list")
  expect_length(weighted, length(ws))
  
  # Test get_hier_clusters_dtc
  # Need to filter to ensure enough spots
  if(nrow(coord_dat) < 5){
    skip("Insufficient spots for clustering test")
  }
  # Pass the full weighted list (not a single matrix)
  dtc_result = get_hier_clusters_dtc(weighted_dists=weighted, ws=ws, deepSplit=2, linkage='ward.D2')
  
  expect_type(dtc_result, "list")
  expect_length(dtc_result, length(ws))
  expect_true(inherits(dtc_result[[1]], "data.frame"))
  expect_true(any(grepl("stclust_spw", colnames(dtc_result[[1]]))))
  
  # Test get_hier_clusters_ks
  ks = 2:3
  ks_result = get_hier_clusters_ks(weighted_dists=weighted, ws=ws, ks=ks, linkage='ward.D2')
  
  expect_type(ks_result, "list")
  expect_length(ks_result, length(ws))
  expect_true(inherits(ks_result[[1]], "data.frame"))
  expect_true(all(c("stclust_spw0.025_k2", "stclust_spw0.025_k3") %in% colnames(ks_result[[1]])))
})


test_that("STclust error handling works", {
  # Unload installed package to avoid conflicts
  if("package:spatialGE" %in% search()) {
    detach("package:spatialGE", unload=TRUE, force=TRUE)
  }
  # Load local development version
  devtools::load_all('../../.', export_all=TRUE)
  
  # Test invalid input
  expect_error(STclust(), "input must be a STlist")
  
  # Test invalid weights
  expect_error(STclust(x=STlist(), ws=-0.1), "spatial weight between 0 and 1")
  expect_error(STclust(x=STlist(), ws=1.5), "spatial weight between 0 and 1")
  
  # Test invalid ks
  expect_error(STclust(x=STlist(), ks=1), "Refusing to generate")
})