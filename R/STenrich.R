##
# STenrich: Spatial enrichment of gene expression sets
#
# This file contains the public STenrich function (modular implementation)
# and imports from STenrich_helpers.R for modular helper functions
#
# @importFrom magrittr %>%
# @importFrom stats p.adjust
# @importFrom tibble add_column column_to_rownames
# @importFrom dplyr mutate select
#

# Source helper functions
helper_path <- system.file('R', 'STenrich_helpers.R', package='spatialGE', mustWork=FALSE)
if(helper_path == '' || !file.exists(helper_path)){
  source('R/STenrich_helpers.R')
} else {
  source(helper_path)
}

##
#' @title STenrich: Test for spatial enrichment of gene expression sets
#' @description Test for spatial enrichment of gene expression sets in spatial transcriptomics data
#' @details The function performs a randomization test to assess if the sum of
#' distances between cells/spots with high expression of a gene set is lower than
#' the sum of distances among randomly selected cells/spots. The cells/spots are
#' considered as having high gene set expression if the average expression of genes in a
#' set is higher than the average expression plus `num_sds` times the standard deviation.
#' Control over the size of regions with high expression is provided by setting the
#' minimum number of cells/spots (`min_units`). This method is a modification of
#' the method devised by Hunter et al. 2021 (zebrafish melanoma study).
#'
#' @param x an STlist with transformed gene expression
#' @param samples a vector with sample names or indexes to run analysis
#' @param gene_sets a named list of gene sets to test. The names of the list should
#' identify the gene sets to be tested
#' @param score_type Controls how gene set expression is calculated. The options are
#' the average expression among genes in a set ('avg'), or a GSEA score ('gsva'). The
#' default is 'avg'
#' @param reps the number of random samples to be extracted. Default is 1000 replicates
#' @param annot name of the annotation within `x@spatial_meta` containing the spot/cell
#' categories. Needs to be used in conjunction with `domain`
#' @param domain the domain to restrict the analysis. Must exist within the spot/cell
#' categories included in the selected annotation (i.e., `annot`)
#' @param num_sds the number of standard deviations to set the minimum gene set
#' expression threshold. Default is one (1) standard deviation
#' @param min_units Minimum number of spots with high expression of a pathway for
#' that gene set to be considered in the analysis. Defaults to 20 spots or cells
#' @param min_genes the minimum number of genes of a gene set present in the data set
#' for that gene set to be included. Default is 5 genes
#' @param pval_adj_method the method for multiple comparison adjustment of p-values.
#' Options are the same as that of `p.adjust`. Default is 'BH'
#' @param seed the seed number for the selection of random samples. Default is 12345
#' @param cores the number of cores used during parallelization. If NULL (default),
#' the number of cores is defined automatically
#' @param verbose logical, whether to print text to console
#' @return a list of data frames with the results of the test
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Using included melanoma example (Thrane et al.)
#' # Download example data set from spatialGE_Data
#' thrane_tmp = tempdir()
#' unlink(thrane_tmp, recursive=TRUE)
#' dir.create(thrane_tmp)
#' lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
#' tryCatch({
#'   download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
#'   unzip(file=paste0(thrane_tmp, '/melanoma_thrane.zip'), exdir=thrane_tmp)
#'   count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='counts')
#'   coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='mapping')
#'   clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='clinical')
#'   st_obj <- STlist(rnacounts=count_files, spotcoords=coord_files, samples=clin_file)
#'   st_obj <- transform_data(st_obj)
#'
#'   # Define gene sets (using first 10 genes from the data as example)
#'   genes <- rownames(st_obj@tr_counts[[1]])
#'   gene_sets <- list(
#'     GS1 = genes[1:5],
#'     GS2 = genes[6:10]
#'   )
#'
#'   # Run STenrich
#'   res <- STenrich(st_obj, gene_sets = gene_sets, samples = 1, reps = 100)
#' }
#' }
#'
#' @seealso STenrich_legacy for reproducibility with original implementation
#'
#' @importFrom magrittr %>%
#' @importFrom stats p.adjust
#' @importFrom tibble add_column column_to_rownames
#' @importFrom dplyr mutate select
#' @importFrom parallel mclapply count
STenrich = function(x, samples=NULL, gene_sets=NULL, score_type='avg', reps=1000,
                    annot=NULL, domain=NULL, num_sds=1, min_units=20, min_genes=5,
                    pval_adj_method='BH', seed=12345, cores=NULL, verbose=TRUE){

  # Record time
  zero_t = Sys.time()
  
  if(verbose){
    cat("Running STenrich...\n")
  }
  
  # Validate inputs
  validated = STenrich_validate_input(x, samples, gene_sets, score_type, annot, domain, num_sds, min_units, min_genes, pval_adj_method)
  x = validated$x
  samples = validated$samples
  gene_sets = validated$gene_sets
  score_type = validated$score_type
  annot = validated$annot
  domain = validated$domain
  num_sds = validated$num_sds
  min_units = validated$min_units
  min_genes = validated$min_genes
  pval_adj_method = validated$pval_adj_method
  
  # Prepare data
  if(verbose){
    cat('\tPreparing data...\n')
  }
  data_res = STenrich_prepare_data(x, samples, gene_sets, annot, domain, min_units)
  
  # Calculate gene set expression
  if(verbose){
    cat('\tCalculating gene set expression...\n')
  }
  
  # Define number of cores for parallelization
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(samples))
    user_cores = F
  } else{
    cores = ceiling(cores)
    user_cores = T
  }
  
  # Calculate expression
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
  
  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('STenrich completed in ', round(end_t, 2), ' min.\n'))
  }
  
  return(pval_res)
}


##
#' @title STenrich_calculate_gs_mean_exp: Calculate average gene set expression
#' @description Calculate mean expression across genes in each gene set
#' @param x an STlist with transformed gene expression (as DelayedArray)
#' @param combo data frame with combinations of samples and gene sets
#' @param pw_genes list of available genes per gene set
#' @param min_genes the minimum number of genes of a gene set present in data (default: 5)
#' @param cores the number of cores for parallelization
#' @return list of data frames with average gene set expression
#' @keywords internal
STenrich_calculate_gs_mean_exp = function(x, combo, pw_genes, min_genes, cores){
  
  result_df = parallel::mclapply(1:nrow(combo), function(i){
    sample_tmp = as.vector(combo[[1]][i])
    geneset_tmp = as.vector(combo[[2]][i])
    pw_genes_tmp = pw_genes[[i]]
    
    if(length(pw_genes_tmp) >= min_genes){
      expr_subset = x[[sample_tmp]][pw_genes_tmp, ]
      pw_avg_exp = DelayedMatrixStats::colMeans2(expr_subset)
    } else{
      pw_avg_exp = rep(NA, ncol(x[[sample_tmp]]))
      names(pw_avg_exp) = colnames(x[[sample_tmp]])
    }
    
    pw_avg_exp = data.frame(as.list(pw_avg_exp), check.names=F)
    rownames(pw_avg_exp) = geneset_tmp
    return(pw_avg_exp)
  }, mc.cores=cores)
  names(result_df) = paste0(combo[[1]], '&&', combo[[2]])
  
  result_df = lapply(1:length(unique(combo[[1]])), function(i){
    sample_tmp = unique(combo[[1]])[i]
    idx = which(combo[[1]] == sample_tmp)
    df_tmp = do.call(rbind, result_df[idx])
    rownames(df_tmp) = gsub(paste0('^', sample_tmp, '&&'), '', rownames(df_tmp))
    return(df_tmp)
  })
  names(result_df) = unique(combo[[1]])
  
  return(result_df)
}


##
#' @title STenrich_calculate_gs_gsva_score: Calculate GSVA scores
#' @description Calculate GSVA enrichment scores for gene sets
#' @param x an STlist with transformed gene expression (as DelayedArray)
#' @param pw_genes list of available genes per gene set
#' @param gene_sets a named list of gene sets to test
#' @param min_genes the minimum number of genes of a gene set present in data (default: 5)
#' @param cores the number of cores for parallelization
#' @param verbose verbosity level
#' @return list of data frames with GSVA scores
#' @keywords internal
STenrich_calculate_gs_gsva_score = function(x, pw_genes, gene_sets, min_genes, cores, verbose=TRUE){
  
  samples = as.vector(unique(combo[[1]]))
  result_df = lapply(1:length(samples), function(i){
    sample_tmp = as.vector(unique(samples))[i]
    
    if(!is.null(pw_genes)){
      pw_genes_tmp = pw_genes[combo[[1]] == sample_tmp]
      names(pw_genes_tmp) = as.vector(combo[[2]][combo[[1]] == sample_tmp])
    }
    gene_sets_tmp = gene_sets[ names(pw_genes_tmp)[unlist(lapply(pw_genes_tmp, length)) >= min_genes] ]
    
    if(length(gene_sets_tmp) > 0){
      gsvapar = GSVA::gsvaParam(as.array(x[[sample_tmp]]), geneSets=gene_sets_tmp)
      pw_gsva_exp = GSVA::gsva(gsvapar, BPPARAM=BiocParallel::MulticoreParam(workers=cores), verbose=verbose)
      pw_gsva_exp = as.data.frame(pw_gsva_exp)
    } else{
      pw_gsva_exp = data.frame(matrix(nrow=length(pw_genes_tmp), ncol=ncol(x[[sample_tmp]])))
      rownames(pw_gsva_exp) = names(pw_genes_tmp)
      colnames(pw_gsva_exp) = colnames(x[[sample_tmp]])
    }
    
    if(length(gene_sets_tmp) != length(gene_sets)){
      pw_gsva_exp_list = lapply(names(gene_sets), function(j){
        if(j %in% rownames(pw_gsva_exp)){
          return(pw_gsva_exp[j, , drop=F])
        } else{
          df_tmp = as.data.frame(as.list(setNames(rep(NA, ncol(pw_gsva_exp)), colnames(pw_gsva_exp))), check.names=F)
          rownames(df_tmp) = j
          return(df_tmp)
        }
      })
      pw_gsva_exp = do.call(rbind, pw_gsva_exp_list)
    }
    
    return(pw_gsva_exp)
  })
  names(result_df) = unique(samples)
  
  return(result_df)
}


##
#' @title STenrich_permutation_test: Perform permutation testing
#' @description Calculate observed distances and perform random permutations
#' @param result_df list of data frames with expression scores
#' @param coords_df list of coordinate matrices for each sample
#' @param combo data frame with combinations of samples and gene sets
#' @param pw_genes list of available genes per gene set
#' @param samples a vector with sample names to run analysis
#' @param gene_sets a named list of gene sets to test
#' @param num_sds number of standard deviations to set minimum gene set expression threshold (default: 1)
#' @param min_units Minimum number of spots with high expression (default: 20)
#' @param reps the number of random samples to be extracted (default: 1000)
#' @param seed the seed number for random sampling (default: 12345)
#' @param cores the number of cores for parallelization
#' @param verbose verbosity level
#' @return list of data frames with p-values per sample
#' @keywords internal
STenrich_permutation_test = function(result_df, coords_df, combo, pw_genes, samples, gene_sets, num_sds, min_units, reps, seed, cores, verbose){
  
  pval_res = parallel::mclapply(1:nrow(combo), function(i){
    sample_tmp = as.vector(combo[i, 1])
    gs_tmp = as.vector(combo[i, 2])
    expr_vals = unlist(result_df[[sample_tmp]][gs_tmp, ])
    
    res_df = data.frame(sample_name=sample_tmp, gene_set=gs_tmp,
                        size_test=length(pw_genes[[i]]),
                        size_gene_set=length(gene_sets[[gs_tmp]]),
                        p_value=NA)
    
    if(!all(is.na(expr_vals))){
      if(verbose >= 2L){
        system(sprintf('echo "%s"', paste0("\tTesting sample ", sample_tmp, ", ", gs_tmp, "...")))
      }
      
      exp_thresh = mean(expr_vals, na.rm=T) + (num_sds*sd(expr_vals, na.rm=T))
      high_spots_bc = names(which(expr_vals >= exp_thresh))
      
      if(length(high_spots_bc) >= min_units){
        coords_df_tmp = coords_df[[sample_tmp]][high_spots_bc, ]
        sum_high_distances = computeSubsampleSums(coords_df_tmp, n_subsample=nrow(coords_df_tmp), n_samples=1)
        rm(coords_df_tmp)
        
        set.seed(seed)
        sum_rand_distances = computeSubsampleSums(coords_df[[sample_tmp]], n_subsample=length(high_spots_bc), n_samples=reps)
        
        count_test = sum(sum_rand_distances < sum_high_distances)
        p_val = count_test / reps
        
        res_df[['p_value']] = p_val
        
        rm(high_spots_bc, sum_high_distances, sum_rand_distances, count_test, p_val)
      }
    }
    
    return(res_df)
  }, mc.cores=cores)
  
  pval_res = do.call(rbind, pval_res)
  
  return(pval_res)
}


##
#' @title STenrich_format_results: Format final results
#' @description Compile results and adjust p-values
#' @param pval_res data frame with p-values per sample and gene set
#' @param samples a vector with sample names to run analysis
#' @param pval_adj_method the method for multiple comparison adjustment (default: 'BH')
#' @return list of data frames with p-values and adjusted p-values per sample
#' @keywords internal
STenrich_format_results = function(pval_res, samples, pval_adj_method){
  
  pval_res = pval_res %>%
    tibble::add_column(prop_size_test=.[['size_test']]/.[['size_gene_set']], .before='p_value') %>%
    dplyr::mutate(prop_size_test=round(prop_size_test, 3))
  
  sample_names_tmp = unique(pval_res[['sample_name']])
  pval_res = lapply(sample_names_tmp, function(i){
    df_tmp = pval_res[pval_res[['sample_name']] == i, ]
    df_tmp[['adj_p_value']] = p.adjust(df_tmp[['p_value']], method=pval_adj_method)
    df_tmp = df_tmp[order(df_tmp[['adj_p_value']]), ]
    return(df_tmp)
  })
  names(pval_res) = sample_names_tmp
  
  return(pval_res)
}