##
# STclust Helper Functions
# Utility functions for spatial clustering

##
#' @title calculate_dist_matrices: Calculate and scale distance matrices
#' @description Calculate Euclidean distance matrices for gene expression and spatial coordinates,
#' then scale both to [0,1] range
#' @param expr_dat a sparse matrix with gene expression (genes in rows, spots/cells in columns)
#' @param coord_dat a data frame with columns: 'libname', 'xpos', 'ypos'
#' @param dist_metric character string indicating distance metric (e.g., 'euclidean')
#' @return a list with two matrices: `scale_exp` (scaled expression distances) and 
#' `scale_coord` (scaled coordinate distances)
#' @keywords internal
#' @importFrom Matrix t
#' @importFrom dplyr %>%
#' @importFrom tibble column_to_rownames
#' @importFrom wordspace dist.matrix
#' @importFrom stats dist as.matrix
calculate_dist_matrices = function(expr_dat=NULL, coord_dat=NULL, dist_metric=NULL){
  a = Matrix::t(expr_dat)
  b = coord_dat[, c('libname', 'xpos', 'ypos')] %>% tibble::column_to_rownames(var='libname')
  b = as.matrix(b[match(rownames(b), rownames(a)), ])

  # Get distance matrices
  da = wordspace::dist.matrix(a, method=dist_metric)
  db = dist(b, upper=T, diag=T, method=dist_metric)
  dam = as.matrix(da)
  dbm = as.matrix(db)

  rm(a, b, da, db) # Clean env

  # Scale matrices to [0,1]
  dam = dam/max(dam)
  dbm = dbm/max(dbm)

  return(list(scale_exp=dam, scale_coord=dbm))
}


##
#' @title calculate_weighted_dist: Calculate weighted distance matrices
#' @description Combine scaled expression and spatial distance matrices using user-specified weights
#' @param scaled_dists a list with two matrices: scaled expression and spatial distance matrices
#' @param ws vector of spatial weights (0-1) to apply to spatial distances
#' @return a list of weighted distance matrices, one for each ws value
#' @keywords internal
#' @importFrom utils combn
calculate_weighted_dist = function(scaled_dists=NULL, ws=NULL){
  weight_mtx_ls = lapply(1:length(ws), function(w){
    weight_d = ws[w]
    weight_g = 1-weight_d

    # Create vector of weights for Reduce
    weight_ls = c(weight_g, weight_d)
    dmxs = list(scaled_dists[[1]], scaled_dists[[2]])

    # Apply weight element-wise
    weight_mtx = Reduce('+', Map('*', dmxs, weight_ls))

    return(weight_mtx)
  })

  return(weight_mtx_ls)
}


##
#' @title get_hier_clusters_dtc: Hierarchical clustering with DynamicTreeCut
#' @description Perform hierarchical clustering followed by DynamicTreeCut for adaptive cluster detection
#' @param weighted_dists a list of distance matrices (NOT dist objects) for each spatial weight
#' @param ws a vector with spatial weights
#' @param deepSplit a logical or integer (1-4) for cutreeDynamic deepSplit parameter
#' @param linkage a string with the linkage method for hclust (e.g., 'ward.D2')
#' @return a list of data frames with spot/cell cluster assignments for each weight
#' @keywords internal
#' @importFrom dplyr %>%
#' @importFrom tibble tibble add_column
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom stats as.dist hclust
get_hier_clusters_dtc = function(weighted_dists=NULL, ws=NULL, deepSplit=NULL, linkage=NULL){
  grp_df_ls = lapply(1:length(ws), function(w){
    # Construct column name to be put in `spatial_meta` based on weight and deepSplit
    if(is.logical(deepSplit)){
      dspl = 'False'
      if(deepSplit){
        dspl = 'True'
      }
    } else{
      dspl = deepSplit
    }
    col_name = paste0('stclust_spw', ws[w], '_dspl', dspl)

    # Run hierarchical clustering
    hierclusters = hclust(as.dist(weighted_dists[[w]]), method=linkage)

    # Use DynamicTreeClusters
    grp_df = dynamicTreeCut::cutreeDynamic(hierclusters, method='hybrid', distM=weighted_dists[[w]], deepSplit=deepSplit, verbose=F)

    # Create data frame with cluster assignments
    grp_df = tibble::tibble(libname=colnames(weighted_dists[[w]]), !!col_name:=as.factor(as.vector(grp_df)))
    # Convert zeroes to NAs
    grp_df[[col_name]][grp_df[[col_name]] == 0] = NA

    return(grp_df)
  })

  return(grp_df_ls)
}


##
#' @title get_hier_clusters_ks: Hierarchical clustering with fixed k values
#' @description Perform hierarchical clustering followed by cutree for fixed number of clusters
#' @param weighted_dists a list of distance matrices (NOT dist objects) for each spatial weight
#' @param ws a vector with spatial weights
#' @param ks a vector with k values for cluster detection
#' @param linkage a string with the linkage method for hclust (e.g., 'ward.D2')
#' @return a list of data frames with spot/cell cluster assignments for each weight
#' @keywords internal
#' @importFrom dplyr %>%
#' @importFrom tibble tibble add_column
#' @importFrom stats as.dist hclust cutree
get_hier_clusters_ks = function(weighted_dists=NULL, ws=NULL, ks=NULL, linkage=NULL){

  # To prevent NOTES in R CMD check
  . = NULL

  grp_df_ls = lapply(1:length(ws), function(w){
    grp_df = tibble::tibble(libname=colnames(weighted_dists[[w]]))
    for(k in ks){
      # Construct column name to be put in `spatial_meta` based on weight and k
      col_name = paste0('stclust_spw', ws[w], '_k', k)

      # Run hierarchical clustering
      hierclusters = hclust(as.dist(weighted_dists[[w]]), method=linkage)

      # Cut the dendrogram
      grp_df_tmp = cutree(hierclusters, k=k)
      # Create data frame with cluster assignments
      grp_df_tmp = tibble::tibble(libname=colnames(weighted_dists[[w]]), !!col_name:=as.factor(as.vector(grp_df_tmp)))

      grp_df = grp_df %>% dplyr::left_join(., grp_df_tmp, by='libname')

      rm(grp_df_tmp, hierclusters, col_name) # Clean env
    }
    return(grp_df)
  })

  return(grp_df_ls)
}