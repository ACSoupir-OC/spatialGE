##
# STenrich Core Implementation
# Core workflow logic for spatial enrichment analysis
#
# This file contains the main STenrich workflow that orchestrates
# the enrichment analysis by calling helper functions.
#

##
#' @title STenrich_core: Core workflow for spatial enrichment analysis
#' @description Main workflow function that orchestrates the enrichment analysis
#' @details This function performs the core workflow:
#' 1. Prepare data (extract tissue spots, coordinates, match gene sets)
#' 2. Calculate gene set expression (average or GSVA)
#' 3. Perform permutation testing
#' 4. Format results with adjusted p-values
#'
#' @param x an STlist with transformed gene expression
#' @param samples a vector with sample names or indexes to run analysis
#' @param gene_sets a named list of gene sets to test
#' @param score_type Controls how gene set expression is calculated. Options: 'avg' or 'gsva'
#' @param reps the number of random samples to be extracted. Default is 1000 replicates
#' @param annot name of the annotation within `x@spatial_meta` containing spot/cell categories
#' @param domain the domain to restrict the analysis
#' @param num_sds number of standard deviations to set minimum gene set expression threshold (default: 1)
#' @param min_units Minimum number of spots with high expression (default: 20)
#' @param min_genes the minimum number of genes of a gene set present in data (default: 5)
#' @param pval_adj_method the method for multiple comparison adjustment (default: 'BH')
#' @param seed the seed number for random sampling (default: 12345)
#' @param cores the number of cores for parallelization. If NULL, auto-detected
#' @param verbose logical, whether to print text to console
#' @return list of data frames with p-values and adjusted p-values per sample
#' @keywords internal
STenrich_core = function(x, samples, gene_sets, score_type='avg', reps=1000,
                         annot=NULL, domain=NULL, num_sds=1, min_units=20, min_genes=5,
                         pval_adj_method='BH', seed=12345, cores=NULL, verbose=TRUE){
  
  # Prepare data
  if(verbose){
    cat('\tPreparing data...\n')
  }
  data_res = STenrich_prepare_data(x, samples, gene_sets, annot, domain, min_units)
  
  # Update samples after filtering
  samples = data_res$samples
  
  # Define number of cores for parallelization
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(samples))
    user_cores = FALSE
  } else{
    cores = ceiling(cores)
    user_cores = TRUE
  }
  
  # Calculate gene set expression
  if(verbose){
    cat('\tCalculating gene set expression...\n')
  }
  
  combo = data_res$combo
  pw_genes = data_res$pw_genes
  delayed_x = lapply(samples, function(i){
    mtx_tmp = x@tr_counts[[i]]
    if(length(data_res$tissue_spots) > 0){
      mtx_tmp = mtx_tmp[, data_res$tissue_spots[[i]]]
    }
    DelayedArray::DelayedArray(mtx_tmp)
  })
  names(delayed_x) = samples
  
  if(score_type == 'avg'){
    if(verbose){
      cat("\tCalculating average gene set expression...\n")
    }
    if(!user_cores){
      cores = count_cores(nrow(combo))
    }
    result_df = STenrich_calculate_gs_mean_exp(delayed_x, combo, pw_genes, min_genes, cores)
  } else if(score_type == 'gsva'){
    if(verbose){
      cat("\tCalculating GSVA score...\n")
    }
    if(!user_cores){
      cores = count_cores(length(gene_sets))
    }
    result_df = STenrich_calculate_gs_gsva_score(delayed_x, pw_genes, gene_sets, min_genes, cores, verbose)
  }
  
  # Perform permutation testing
  if(verbose){
    cat('\tPerforming permutation testing...\n')
  }
  
  if(!user_cores){
    cores = count_cores(nrow(combo))
  }
  pval_res = STenrich_permutation_test(result_df, data_res$coords_df, combo, pw_genes, samples, 
                                       gene_sets, num_sds, min_units, reps, seed, cores, verbose)
  
  # Format results
  pval_res = STenrich_format_results(pval_res, samples, pval_adj_method)
  
  return(pval_res)
}
