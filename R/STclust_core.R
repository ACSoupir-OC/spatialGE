##
# STclust Core Modular Functions
# Main refactored clustering functions

##
#' @title STclust_select_genes: Select variable genes for clustering
#' @description Selects top variable genes based on Seurat's VST method and filters expression data
#' @param x an STlist object with normalized expression data
#' @param samples vector of sample names to process
#' @param topgenes number of top variable genes to select per sample
#' @param cores number of cores for parallelization
#' @return an STlist with updated `@gene_meta` (VST variance) and `@tr_counts` filtered to top genes
#' @keywords internal
STclust_select_genes = function(x, samples, topgenes, cores){

  # Identify variable genes (Seurat's VST)
  x = calculate_vst(x=x, samples=samples, cores=cores)

  # Subset variable genes
  trcounts_df = parallel::mclapply(samples, function(i){
    topgenenames_tmp = x@gene_meta[[i]] %>%
      dplyr::arrange(dplyr::desc(vst.variance.standardized)) %>%
      dplyr::slice_head(n=topgenes) %>%
      dplyr::select(gene) %>%
      unlist() %>%
      as.vector()

    # Get transformed counts
    trcounts_df_tmp = x@tr_counts[[i]][rownames(x@tr_counts[[i]]) %in% topgenenames_tmp, ]

    return(trcounts_df_tmp)
  }, mc.cores=cores)
  names(trcounts_df) = samples

  return(list(x=x, trcounts_df=trcounts_df))
}


##
#' @title STclust_calculate_distances: Calculate scaled distance matrices
#' @description Calculate and scale gene expression and spatial distance matrices for a single sample
#' @param trcounts_df a filtered expression matrix for one sample
#' @param coord_dat a data frame with spot coordinates for one sample
#' @param dist_metric distance metric to use
#' @return a list with scaled distance matrices (expression + spatial)
#' @keywords internal
STclust_calculate_distances = function(trcounts_df, coord_dat, dist_metric){
  scaled_dists = calculate_dist_matrices(expr_dat=trcounts_df, coord_dat=coord_dat, dist_metric=dist_metric)
  return(scaled_dists)
}


##
#' @title STclust_weight_distances: Create weighted distance matrices
#' @description Combine scaled expression and spatial distance matrices using user-specified weights
#' @param scaled_dists list of scaled distance matrices (expression + spatial)
#' @param ws vector of spatial weights to apply
#' @return list of weighted distance matrices
#' @keywords internal
STclust_weight_distances = function(scaled_dists, ws){
  weighted_dists = calculate_weighted_dist(scaled_dists=scaled_dists, ws=ws)
  return(weighted_dists)
}


##
#' @title STclust_hierarchical: Core hierarchical clustering logic
#' @description Performs hierarchical clustering on weighted distance matrices using either DTC or fixed k
#' @param weighted_dists list of weighted distance matrices for one sample (one per ws value)
#' @param ws vector of spatial weights
#' @param ks clustering method: 'dtc' for DynamicTreeCut or numeric vector for fixed k values
#' @param linkage linkage method for hclust (e.g., 'ward.D2')
#' @param deepSplit deepSplit parameter for cutreeDynamic (if ks='dtc')
#' @param verbose verbosity level
#' @return list of data frames with cluster assignments for each ws value
#' @keywords internal
STclust_hierarchical = function(weighted_dists, ws, ks, linkage, deepSplit, verbose){

  # Identify clustering method... ONLY HIERARCHICAL CLUSTERING IMPLEMENTED SO FAR
  clmethod = 'hclust'

  if(clmethod == 'hclust'){

    # Apply dtc or split to k
    if(as.character(ks[1]) == 'dtc'){
      # Hierarchical clustering using DynamicTreeClusters
      hierclusters_ls = get_hier_clusters_dtc(weighted_dists=weighted_dists, ws=ws, deepSplit=deepSplit, linkage=linkage)
    } else if(is.numeric(ks)){
      # Hierarchical clustering using range of Ks
      hierclusters_ls = get_hier_clusters_ks(weighted_dists=weighted_dists, ws=ws, ks=ks, linkage=linkage)
    } else{
      stop('Enter a valid number of k values to evaluate or \'dtc\' to apply cutreeDynamic.')
    }

  } else{
    stop('Currently, only spatially-informed hierarchical clustering is supported.')
  }

  return(hierclusters_ls)
}


##
#' @title STclust: Spatial clustering (modular implementation)
#' @description
#' `r lifecycle::badge("stable")`
#' 
#' Perform unsupervised spatially-informed clustering on the spots/cells of an ST sample
#' @details
#' The function calculates Euclidean distances between cells or spots based on:
#' 1. Expression of top variable genes (Seurat's FindVariableFeatures)
#' 2. Spatial coordinates (x,y)
#' 
#' These distances are weighted and combined, then hierarchical clustering is performed.
#' Supports two clustering modes:
#' - DTC (DynamicTreeCut): Adaptive cluster detection
#' - Fixed k: User-specified number of clusters
#' 
#' \strong{Note:} This function was refactored in version 2.0.0 to use a modular architecture
#' with 5 helper functions (\code{STclust_select_genes}, \code{STclust_calculate_distances}, etc.).
#' The original monolithic implementation is available as \code{STclust_legacy()} for reproducibility
#' with version 1.x results.
#'
#' @param x an STlist with normalized expression data
#' @param samples a vector with strings or integers indicating samples to cluster
#' @param ws a double (0-1) indicating weight for spatial distances (default: 0.025)
#' @param dist_metric distance metric to use (default: 'euclidean')
#' @param linkage linkage method for hclust (default: 'ward.D2')
#' @param ks clustering method: 'dtc' for DynamicTreeCut, or numeric vector for fixed k values
#' @param topgenes number of top variable genes to use (default: 2000)
#' @param deepSplit deepSplit parameter for cutreeDynamic (if ks='dtc')
#' @param cores number of cores for parallelization
#' @param verbose verbosity level (0, 1, or 2)
#' @return an STlist with cluster assignments added to @spatial_meta
#'
#' @export
#'
#' @examples
#' \donttest{
#' # TNBC dataset (Bassiouni et al.)
#' tnbc_tmp = tempdir()
#' unlink(tnbc_tmp, recursive=TRUE)
#' dir.create(tnbc_tmp)
#' lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/tnbc_bassiouni.zip?download='
#' tryCatch({
#'   download.file(lk, destfile=file.path(tnbc_tmp, 'tnbc.zip'), mode='wb')
#'   unzip(file.path(tnbc_tmp, 'tnbc.zip'), exdir=tnbc_tmp)
#'   tnbc = STlist(rnacounts=list.files(file.path(tnbc_tmp, 'tnbc_bassiouni'), pattern='counts', full.names=TRUE)[1],
#'                 spotcoords=list.files(file.path(tnbc_tmp, 'tnbc_bassiouni'), pattern='mapping', full.names=TRUE)[1],
#'                 samples='TNBC')
#'   tnbc = transform_data(tnbc)
#'   tnbc = STclust(tnbc, ks=3, ws=0.025)
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#' })
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom methods as is new
#' @importFrom stats as.dist complete.cases cutree dist hclust prcomp sd na.omit
#' @importFrom parallel mclapply detectCores
STclust = function(x=NULL, samples=NULL, ws=0.025, dist_metric='euclidean', linkage='ward.D2', ks='dtc', topgenes=2000, deepSplit=FALSE, cores=NULL, verbose=TRUE){

  # To prevent NOTES in R CMD check
  . = NULL

  # Record time
  zero_t = Sys.time()
  verbose = 1L
  if(verbose){
    cat(paste0('STclust started...\n'))
  }

  # Force ws and topgenes as numeric in case entered as character
  ws = as.double(ws)
  topgenes = as.integer(topgenes)

  # Do not allow weights higher than 1
  if(any(ws < 0) | any(ws > 1)){
    stop('Please select a spatial weight between 0 and 1.')
  }

  # Check to ensure number of ks is acceptable
  if(is.numeric(ks)){
    ks = as.integer(ks)
    if(length(ks) == 1 & ks[1] < 2){
      raise_err(err_code='error0016')
    } else if(any(ks < 2)){
      warning('Refusing to generate < 2 clusters. Skipping any k < 2.')
      ks = ks[ks >= 2]
    }
  }

  # Test if an STList has been input
  if(is.null(x) | !is(x, 'STlist')){
    stop("The input must be a STlist.")
  }

  # Check data has been normalized
  if(length(x@tr_counts) < 1){
    raise_err(err_code='error0007')
  }

  # Define samples using names (convert indexes to names if necessary)
  if(is.null(samples)){
    samples = names(x@spatial_meta)
  } else{
    if(is.numeric(samples)){
      samples = as.vector(na.omit(names(x@spatial_meta)[samples]))
    } else{
      samples = samples[samples %in% names(x@spatial_meta)]
    }
    # Verify that sample names exist
    if(length(samples) == 0 | !any(samples %in% names(x@spatial_meta))){
      stop('None of the requested samples are present in the STlist.')
    }
  }

  # Define number of cores for parallelization
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(samples))
  } else{
    cores = ceiling(cores)
  }

  # Step 1: Select variable genes
  if(verbose){
    cat('Step 1: Selecting top variable genes...\n')
  }
  select_res = STclust_select_genes(x=x, samples=samples, topgenes=topgenes, cores=cores)
  x = select_res[['x']]
  trcounts_df = select_res[['trcounts_df']]

  # Step 2 & 3: Calculate and weight distances for each sample
  if(verbose){
    cat('Step 2: Calculating distance matrices...\n')
  }

  # Step 4: Hierarchical clustering
  if(verbose){
    cat('Step 3: Performing hierarchical clustering...\n')
  }

  res_ls = parallel::mclapply(samples, function(i){
    # Calculate scaled expression and spatial distance matrices
    scaled_dists = STclust_calculate_distances(trcounts_df=trcounts_df[[i]], coord_dat=x@spatial_meta[[i]], dist_metric=dist_metric)

    # Calculate weighted distance matrices
    weighted_dists = STclust_weight_distances(scaled_dists=scaled_dists, ws=ws)

    rm(scaled_dists) # Clean env

    # Perform hierarchical clustering
    hierclusters_ls = STclust_hierarchical(weighted_dists=weighted_dists, ws=ws, ks=ks, linkage=linkage, deepSplit=deepSplit, verbose=verbose)

    if(verbose > 1L){
      system(sprintf('echo "%s"', paste0("\tClustering completed for ", i, "...")))
    }

    return(hierclusters_ls)
  }, mc.cores=cores)
  names(res_ls) = samples

  rm(trcounts_df) # Clean env

  # Add results to STlist
  if(verbose){
    cat('Step 4: Updating STlist with results...\n')
  }
  # Use lapply to merge results, then assign back to x@spatial_meta
  for(i in samples){
    # Pre-merge all results into a single data frame (5-10% speedup)
    merged_result = res_ls[[i]][[1]]
    if(length(res_ls[[i]]) > 1){
      for(j in 2:length(res_ls[[i]])){
        # Left join to add new columns (more stable than full_join)
        merged_result = dplyr::left_join(merged_result, res_ls[[i]][[j]], by='libname')
      }
    }
    
    # Remove any existing clustering columns first
    if(any(colnames(x@spatial_meta[[i]])[-c(1:5)] %in% colnames(merged_result))){
      col_names = intersect(colnames(x@spatial_meta[[i]])[-1], colnames(merged_result))
      x@spatial_meta[[i]] = dplyr::select(x@spatial_meta[[i]], -!!col_names)
    }
    
    # Single left_join with pre-merged result
    x@spatial_meta[[i]] = dplyr::left_join(x@spatial_meta[[i]], merged_result, by='libname')
  }

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('STclust completed in ', round(end_t, 2), ' min.\n'))
  }

  return(x)
}