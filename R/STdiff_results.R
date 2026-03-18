##
# STdiff Result Compilation Module
# Functions for compiling, formatting, and enriching differential expression results

##
#' @title STdiff_add_metadata: Add cluster and sample metadata to results
#' @description Enriches differential expression results with metadata from meta_dict
#' including original annotation names, sample information, and cluster mappings
#' @param results a list of data frames containing DE results per sample
#' @param meta_dict a data frame with columns (orig_annot, coded_annot) mapping
#' coded annotations to original annotation names
#' @param pairwise logical indicating whether tests are pairwise or reference-based
#' @return enriched list of data frames with complete metadata columns
#' @export
#' @keywords internal
STdiff_add_metadata = function(results=NULL, meta_dict=NULL, pairwise=FALSE){
  
  # To prevent NOTES in R CMD check
  . = NULL
  
  enriched_results = list()
  
  for(sample_name in names(results)){
    res_tmp = results[[sample_name]]
    
    # Add original annotation names for cluster_1
    res_tmp = res_tmp %>%
      dplyr::mutate(cluster_1_orig = meta_dict[['orig_annot']][
        match(cluster_1, meta_dict[['coded_annot']])
      ])
    
    # Add original annotation names for cluster_2 if pairwise
    if(pairwise){
      res_tmp = res_tmp %>%
        dplyr::mutate(cluster_2_orig = meta_dict[['orig_annot']][
          match(cluster_2, meta_dict[['coded_annot']])
        ])
    }
    
    # Add sample metadata
    res_tmp = res_tmp %>%
      dplyr::mutate(sample_name = sample_name)
    
    enriched_results[[sample_name]] = res_tmp
  }
  
  return(enriched_results)
}


##
#' @title STdiff_format_output: Format results for user consumption
#' @description Formats DE results into a user-friendly structure with consistent
#' column naming, proper ordering, and optional filtering
#' @param results a list of data frames containing DE results per sample
#' @param select_cols optional vector of columns to include in output
#' @param filter_sig logical whether to filter for significant results only
#' @param adj_pval_thr p-value threshold for filtering (if filter_sig=TRUE)
#' @param order_by column name to order results by
#' @param ascending logical indicating sort order
#' @return formatted list of data frames ready for user consumption
#' @export
STdiff_format_output = function(results=NULL, select_cols=NULL, 
                                 filter_sig=FALSE, adj_pval_thr=0.05,
                                 order_by='exp_adj_p_val', ascending=TRUE){
  
  # To prevent NOTES in R CMD check
  . = NULL
  
  formatted_results = list()
  
  for(sample_name in names(results)){
    res_tmp = results[[sample_name]]
    
    # Select specified columns or use default set
    if(is.null(select_cols)){
      select_cols = c('sample', 'gene', 'cluster_1', 'cluster_2',
                      'avg_log2fc', 'mm_p_val', 'ttest_p_val', 'wilcox_p_val',
                      'adj_p_val', 'exp_p_val', 'exp_adj_p_val', 'comments')
    }
    
    # Keep only requested columns that exist
    available_cols = intersect(select_cols, colnames(res_tmp))
    res_tmp = res_tmp[, available_cols, drop=FALSE]
    
    # Filter for significant results if requested
    if(filter_sig && 'exp_adj_p_val' %in% colnames(res_tmp)){
      res_tmp = res_tmp %>%
        dplyr::filter(!is.na(exp_adj_p_val) & exp_adj_p_val < adj_pval_thr)
    }
    
    # Order results if specified column exists
    if(order_by %in% colnames(res_tmp)){
      res_tmp = res_tmp %>%
        dplyr::arrange(.data[[order_by]], .data[[order_by]], .by_group=FALSE)
      if(!ascending){
        res_tmp = res_tmp %>%
          dplyr::mutate(!!order_by := -exp_adj_p_val)
      }
    }
    
    formatted_results[[sample_name]] = res_tmp
  }
  
  return(formatted_results)
}


##
#' @title STdiff_compile_results: Compile non-spatial and spatial results
#' @description Main entry point for compiling differential expression results.
#' Merges non-spatial model results (linear models, t-tests, or Wilcoxon tests)
#' with spatial mixed model results containing Matern covariance structure.
#' Handles convergence issues, time-outs, and applies appropriate p-value
#' adjustments. Returns a list of data frames with complete DE statistics
#' per sample.
#'
#' @param prep list from STdiff_select_genes containing:
#'   - non_spatial_results: list of data frames with non-spatial model results
#'   - meta_dict: data frame mapping coded annotations to original names
#'   - pairwise: logical indicating test type
#'   - pval_adj: p-value adjustment method
#' @param spatial list from STdiff_fit_spatial containing:
#'   - sp_models: list of spatial model results (or NULL if not run)
#'   - meta_dict: same meta_dict from prep
#'   - combo_df: combinations data frame
#' @return list of data frames with final differential expression results per sample.
#' Each data frame contains columns: sample, gene, cluster_1, cluster_2 (if pairwise),
#' avg_log2fc, mm_p_val/ttest_p_val/wilcox_p_val, adj_p_val, exp_p_val, exp_adj_p_val,
#' comments (indicating spatial model status: 'no_spatial_test', 'time_out',
#' 'no_convergence', or NA for successful spatial models)
#' @export
#' @examples
#' # Example usage (internal function, called by STdiff())
#' # result_de = STdiff_compile_results(prep=prep, spatial=spatial_res)
#' @keywords internal
STdiff_compile_results = function(prep=NULL, spatial=NULL){
  
  # To prevent NOTES in R CMD check
  . = NULL
  
  result_de = prep[['non_spatial_results']]
  meta_dict = prep[['meta_dict']]
  pairwise = prep[['pairwise']]
  
  # If spatial models were run, merge results
  if(!is.null(spatial) && length(spatial) > 0 && 
     !is.null(spatial[['sp_models']]) && length(spatial[['sp_models']]) > 0){
    
    sp_models = spatial[['sp_models']]
    
    # Compile spatial results
    sp_de = tibble::tibble()
    for(i in names(sp_models)){
      for(cl in 1:length(sp_models[[i]])){
        for(test_name in names(sp_models[[i]][[cl]])){
          
          # Extract results
          if(any(class(sp_models[[i]][[cl]][[test_name]][['spmod']]) == 'HLfit')){
            exp_summ = spaMM::summary.HLfit(sp_models[[i]][[cl]][[test_name]][['spmod']],
                                            details=c(p_value='Wald'), verbose=F)
            exp_summ = exp_summ[['beta_table']]
            exp_warn = NA_character_
          } else if(sp_models[[i]][[cl]][[test_name]][['spmod']] == 'time_out'){
            exp_summ = data.frame(c(NA,NA), c(NA,NA), c(NA,NA), c(NA, -9999))
            exp_warn = 'time_out'
          } else if(sp_models[[i]][[cl]][[test_name]][['spmod']] == 'no_conv'){
            exp_summ = data.frame(c(NA,NA), c(NA,NA), c(NA,NA), c(NA, -9999))
            exp_warn = 'no_convergence'
          } else{
            exp_summ = data.frame(c(NA,NA), c(NA,NA), c(NA,NA), c(NA, -9999))
            exp_warn = 'unknown_error'
          }
          
          meta1_exp_tmp = stringr::str_replace(test_name, 
                                                paste0(i, '_', sp_models[[i]][[cl]][[test_name]][['gene']], '_'), '') %>%
            stringr::str_replace(., '_other$', '')
          if(pairwise){
            meta2_exp_tmp = stringr::str_extract(meta1_exp_tmp, 'c[0-9]+$')
            meta1_exp_tmp = stringr::str_extract(meta1_exp_tmp, '^c[0-9]+')
          }
          
          tibble_tmp = tibble::tibble(sample=i,
                                      gene=sp_models[[i]][[cl]][[test_name]][['gene']],
                                      cluster_1=as.vector(unlist(meta_dict[['orig_annot']][
                                        meta_dict[['coded_annot']] == meta1_exp_tmp
                                      ])),
                                      exp_p_val=exp_summ[[2, 4]],
                                      comments_spatial=exp_warn)
          if(pairwise){
            tibble_tmp = tibble_tmp %>%
              tibble::add_column(cluster_2=as.vector(unlist(meta_dict[['orig_annot']][
                meta_dict[['coded_annot']] == meta2_exp_tmp
              ])), .after='cluster_1')
          }
          sp_de = dplyr::bind_rows(sp_de, tibble_tmp)
        }
      }
    }
    
    # Merge spatial and non-spatial results
    for(sample_name in names(result_de)){
      if(pairwise){
        tibble_tmp = sp_de %>%
          dplyr::select(-c('exp_p_val')) %>%
          dplyr::filter(sample == sample_name) %>%
          dplyr::left_join(., sp_de %>%
                             dplyr::select(-c('comments_spatial')) %>%
                             dplyr::filter(sample == sample_name) %>%
                             dplyr::filter(!is.na(exp_p_val) & exp_p_val != -9999) %>%
                             tibble::add_column(exp_adj_p_val=p.adjust(.[['exp_p_val']],
                                        method=prep[['pval_adj']])),
                           by=c('sample', 'gene', 'cluster_1', 'cluster_2')) %>%
          dplyr::relocate(comments_spatial, .after='exp_adj_p_val')
        result_de[[sample_name]] = dplyr::full_join(result_de[[sample_name]], tibble_tmp,
                                                     by=c('sample', 'gene', 'cluster_1', 'cluster_2'))
      } else{
        tibble_tmp = sp_de %>%
          dplyr::select(-c('exp_p_val')) %>%
          dplyr::filter(sample == sample_name) %>%
          dplyr::left_join(., sp_de %>%
                             dplyr::select(-c('comments_spatial')) %>%
                             dplyr::filter(sample == sample_name) %>%
                             dplyr::filter(!is.na(exp_p_val) & exp_p_val != -9999) %>%
                             tibble::add_column(exp_adj_p_val=p.adjust(.[['exp_p_val']],
                                        method=prep[['pval_adj']])),
                           by=c('sample', 'gene', 'cluster_1')) %>%
          dplyr::relocate(comments_spatial, .after='exp_adj_p_val')
        result_de[[sample_name]] = dplyr::full_join(result_de[[sample_name]], tibble_tmp,
                                                     by=c('sample', 'gene', 'cluster_1'))
      }
      
      result_de[[sample_name]] = result_de[[sample_name]] %>%
        dplyr::mutate(comments=dplyr::case_when(is.na(comments) ~ comments_spatial,
                                                 TRUE ~ comments)) %>%
        dplyr::select(-c('comments_spatial')) %>%
        dplyr::relocate(comments, .after='exp_adj_p_val')
      
      result_de[[sample_name]] = result_de[[sample_name]] %>%
        dplyr::mutate(spatialmod=dplyr::case_when(is.na(exp_p_val) ~ '2',
                                                    TRUE ~ '1'))
      
      if(pairwise){
        result_de[[sample_name]] = result_de[[sample_name]] %>%
          dplyr::arrange(spatialmod, cluster_1, cluster_2, exp_adj_p_val)
      } else{
        result_de[[sample_name]] = result_de[[sample_name]] %>%
          dplyr::arrange(spatialmod, cluster_1, exp_adj_p_val)
      }
      result_de[[sample_name]] = result_de[[sample_name]] %>%
        dplyr::select(-c('spatialmod'))
    }
  }
  
  return(result_de)
}
