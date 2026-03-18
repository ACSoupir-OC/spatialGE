##
# STdiff Spatial Module
# 
# This file contains spatial differential expression testing functions
# extracted from STdiff_helpers.R for better modularity and maintainability.
# The spatial models use spaMM with Matern covariance structure to test
# for spatial autocorrelation in gene expression between clusters.
#
# @importFrom magrittr %>%
# @importFrom rlang :=
# @importFrom spaMM fitme summary.HLfit
# @importFrom stats as.formula
# @importFrom tibble tibble
# @importFrom dplyr filter mutate select rename arrange full_join
# @importFrom stringr str_extract str_replace str_split
# @importFrom parallel mclapply

##
#' @title STdiff_run_spatial: Run spatial differential expression analysis
#' @description Main entry point for spatial mixed model fitting. This function
#' orchestrates the fitting of spatial models with Matern covariance structure
#' for genes selected from non-spatial differential expression testing.
#'
#' The spatial models test for spatial autocorrelation in gene expression
#' between clusters using the formula: exprval ~ meta + Matern(1|xpos + ypos)
#' with fixed nu parameter (nu=0.5) for Matern covariance.
#'
#' @param x an STlist object containing the spatial transcriptomics data
#' @param prep list from STdiff_select_genes containing combo_df, meta_dict,
#'   non_spatial_results, pval_thr, test_type, and pairwise
#' @param sp_topgenes proportion (0-1) of top DE genes to fit spatial models on
#'   (default 0.2 = 20%)
#' @param cores number of CPU cores for parallelization (NULL = auto-detect)
#' @param verbose verbosity level (0 = silent, 1 = progress, 2 = detailed)
#' @return list containing:
#'   - sp_models: nested list of fitted spatial models (HLfit objects) or status strings
#'   - meta_dict: annotation dictionary mapping coded to original cluster names
#'   - combo_df: combinations dataframe for reference
#' @export
#' @keywords spatial differential expression spaMM Matern covariance
#' @examples
#' \dontrun{
#' # Run spatial DE analysis on top 20% of DE genes
#' spatial_res = STdiff_run_spatial(x=stlist_obj, prep=prep_res,
#'                                  sp_topgenes=0.2, cores=4, verbose=1)
#' }
STdiff_run_spatial = function(x=NULL, prep=NULL, sp_topgenes=0.2, cores=NULL, verbose=1L){

  # To prevent NOTES in R CMD check
  . = NULL

  # If sp_topgenes is 0 or negative, return empty result
  if(sp_topgenes <= 0){
    if(verbose >= 1L){
      cat('\tSkipping spatial tests (sp_topgenes <= 0).\n')
    }
    return(list(sp_models=list(), meta_dict=prep[['meta_dict']], combo_df=prep[['combo_df']]))
  }

  if(verbose >= 1L){
    cat('\tRunning spatial tests...\n')
  }

  # Extract components from prep
  non_sp_de_tmp = prep[['non_spatial_results']]
  meta_dict = prep[['meta_dict']]
  combo_df = prep[['combo_df']]
  sp_topgenes = as.double(sp_topgenes)
  
  # Set number of cores
  if(is.null(cores)){
    cores = count_cores(length(non_sp_de_tmp))
  } else{
    cores = ceiling(cores)
  }

  # Initialize spatial models list
  sp_models = list()
  
  # Process each sample
  for(sample_name in names(non_sp_de_tmp)){
    # Subset genes based on adjusted p-values and select top genes for spatial testing
    models_keep = non_sp_de_tmp[[sample_name]] %>%
      dplyr::filter(adj_p_val < prep[['pval_thr']]) %>%
      dplyr::arrange(adj_p_val)
    
    # Group by clusters and select top proportion
    if(prep[['pairwise']]){
      models_keep = models_keep %>%
        dplyr::group_by(cluster_1, cluster_2)
    } else{
      models_keep = models_keep %>%
        dplyr::group_by(cluster_1)
    }
    models_keep = models_keep %>%
      dplyr::slice_head(prop=sp_topgenes) %>%
      dplyr::ungroup()

    # Mark genes for spatial testing by joining with comments
    if(prep[['pairwise']]){
      non_sp_de_tmp[[sample_name]] = non_sp_de_tmp[[sample_name]] %>%
        dplyr::full_join(., models_keep %>%
                           tibble::add_column(comments='spatial_test') %>%
                           dplyr::select(c("sample", "gene", "cluster_1", "cluster_2", "comments")), by=c("sample", "gene", "cluster_1", "cluster_2"))
    } else{
      non_sp_de_tmp[[sample_name]] = non_sp_de_tmp[[sample_name]] %>%
        dplyr::full_join(., models_keep %>%
                           tibble::add_column(comments='spatial_test') %>%
                           dplyr::select(c("sample", "gene", "cluster_1", "comments")), by=c("sample", "gene", "cluster_1"))
    }
    non_sp_de_tmp[[sample_name]] = non_sp_de_tmp[[sample_name]] %>%
      dplyr::mutate(comments=dplyr::case_when(comments == 'spatial_test' ~ NA_character_, TRUE ~ 'no_spatial_test'))

    # Recode annotation names from original to coded format
    for(spotrow in 1:nrow(models_keep)){
      models_keep[['cluster_1']][spotrow] = meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == models_keep[['cluster_1']][spotrow] ]
      if(prep[['pairwise']]){
        models_keep[['cluster_2']][spotrow] = meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == models_keep[['cluster_2']][spotrow] ]
      }
    }

    if(nrow(models_keep) > 0){
      # Create gene x cluster combinations with unique identifiers
      if(prep[['pairwise']]){
        models_keep = models_keep %>%
          dplyr::select(c('sample', 'gene', 'cluster_1', 'cluster_2')) %>%
          tibble::add_column(genexcluster=paste0(.[['sample']], '_', .[['gene']], '_', .[['cluster_1']], '_', .[['cluster_2']]))
      } else{
        models_keep = models_keep %>%
          dplyr::select(c('sample', 'gene', 'cluster_1')) %>%
          tibble::add_column(genexcluster=paste0(.[['sample']], '_', .[['gene']], '_', .[['cluster_1']], '_other'))
      }
      clusters_tmp = unique(models_keep[['cluster_1']])

      # Create expression data frame with coordinates
      expr_df = expandSparse(x@tr_counts[[sample_name]])
      expr_df = t(expr_df) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var='libname') %>%
        dplyr::right_join(x@spatial_meta[[sample_name]] %>%
                            tibble::add_column(group=1, .after='libname') %>%
                            dplyr::select(c('libname', 'group', 'ypos', 'xpos')), . , by='libname') %>%
        dplyr::select(c('libname', 'group', 'ypos', 'xpos'))

      # Run spatial models in parallel for each cluster
      sp_models[[sample_name]] = parallel::mclapply(1:length(clusters_tmp), function(i){
        models_keep_tmp = models_keep %>%
          dplyr::filter(cluster_1 == clusters_tmp[i]) %>%
          dplyr::select('genexcluster') %>%
          unlist() %>% as.vector()
        non_sp_models_tmp = non_sp_de_tmp[[sample_name]][ models_keep_tmp ]

        # Print progress message
        if(verbose == 1L){
          if(prep[['pairwise']]){
            stdout_print = models_keep_tmp %>%
              stringr::str_extract(., paste0('(?<=', clusters_tmp[i], '_)c[0-9]+$')) %>%
              stringr::str_split(., '_') %>% unlist() %>% unique()
            stdout_print_tmp = c()
            for(m2 in stdout_print){
              stdout_print_tmp = append(stdout_print_tmp, meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == m2 ])
            }
            stdout_print = paste0(stdout_print_tmp, collapse=', ')
            system(sprintf('echo "%s"', paste0("\t\tSample: ", sample_name, ", ",
                                               meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == clusters_tmp[i] ],
                                               " (", length(models_keep_tmp), " tests)...")))
          } else {
            system(sprintf('echo "%s"', paste0("\t\tSample: ", sample_name, ", ",
                                               meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == clusters_tmp[i] ],
                                               " (", length(models_keep_tmp), " tests)...")))
          }
        }

        # Run spatial models
        sp_mods = STdiff_fit_spatial_models(expr_dat=expr_df, non_sp_mods=non_sp_models_tmp,
                                            annot_dict=meta_dict, verb=verbose,
                                            prep_pairwise=prep[['pairwise']])
        return(sp_mods)
      }, mc.cores=cores)
    } else{
      if(verbose >= 1L){
        warning(paste0('No genes to test for ', sample_name, '. Try increasing sp_topgenes.'))
      }
    }
  }

  return(list(sp_models=sp_models, meta_dict=meta_dict, combo_df=combo_df))
}


##
#' @title STdiff_fit_spatial_models: Fit spaMM spatial models with Matern covariance
#' @description Fits spatial mixed models using spaMM::fitme with Matern covariance
#' structure to test for spatial autocorrelation in gene expression. The model
#' formula is: exprval ~ meta + Matern(1|xpos + ypos) with fixed nu=0.5.
#'
#' This function handles both pairwise and reference-based comparisons, applying
#' appropriate data filtering and recoding for each model fit.
#'
#' @param expr_dat a data frame with expression data and spatial coordinates
#'   (columns: meta, group, ypos, xpos, and gene expression columns)
#' @param non_sp_mods a list of non-spatial models containing metadata for
#'   each gene-cluster combination to test
#' @param annot_dict a data frame mapping coded annotations to original annotation names
#' @param verb verbosity level (0, 1, or 2)
#' @param prep_pairwise whether tests are pairwise (TRUE) or reference-based (FALSE)
#' @return list of spatial model results with the following structure:
#'   - For successful fits: 'spmod' contains HLfit object
#'   - For timeouts: 'spmod' = 'time_out'
#'   - For convergence failures: 'spmod' = 'no_conv'
#'   - Each element contains: spmod, sample, gene
#' @export
#' @keywords spatial models spaMM Matern covariance mixed models
#' @importFrom spaMM fitme summary.HLfit
#' @importFrom stats as.formula
STdiff_fit_spatial_models = function(expr_dat=NULL, non_sp_mods=NULL, annot_dict=NULL,
                                     verb=NULL, prep_pairwise=NULL){

  # To prevent NOTES in R CMD check
  . = NULL

  res_ls = list()
  
  # Iterate over each gene-cluster combination
  for(i in names(non_sp_mods)){
    # Extract metadata from non-spatial models
    sample_tmp = non_sp_mods[[i]][['samplename']]
    gene_tmp = non_sp_mods[[i]][['gene']]
    meta1_tmp = non_sp_mods[[i]][['meta1']]
    meta2_tmp = non_sp_mods[[i]][['meta2']]

    # Print progress information for verbose level 2
    if(verb == 2L){
      if(grepl('_c[0-9]+_other$', i)){
        stdout_print = i %>%
          stringr::str_extract(., 'c[0-9]+_other$') %>%
          stringr::str_replace(., '_other$', '')
        stdout_print = annot_dict[['orig_annot']][ annot_dict[['coded_annot']] == stdout_print ]
      } else{
        stdout_print = i %>%
          stringr::str_extract(., paste0('c[0-9]+_c[0-9]+$')) %>%
          stringr::str_split(., '_') %>% unlist() %>% unique()
        stdout_print = c(annot_dict[['orig_annot']][ annot_dict[['coded_annot']] == stdout_print[1] ],
                         annot_dict[['orig_annot']][ annot_dict[['coded_annot']] == stdout_print[2] ])
        stdout_print = paste0(stdout_print[1], ' vs. ', stdout_print[2])
      }
      system(sprintf('echo "%s"', paste0('\t\tSample: ', sample_tmp, '; ', gene_tmp, ', ', stdout_print)))
    }

    # Subset expression data for this gene
    expr_subset = expr_dat[, c("meta" , "group", "ypos", "xpos", gene_tmp)] %>%
      dplyr::rename(exprval := !!gene_tmp)
    
    # Apply filtering based on pairwise or reference-based comparison
    if(prep_pairwise){
      expr_subset = expr_subset[ expr_subset[['meta']] %in% c(meta1_tmp, meta2_tmp), ]
    } else{
      expr_subset[['meta']][ expr_subset[['meta']] != meta1_tmp ] = 'other'
    }

    # Fit spatial mixed model with Matern covariance
    exp_out = tryCatch({
      spaMM::fitme(formula=stats::as.formula(paste0("exprval~meta+Matern(1|xpos+ypos)")),
                   data=expr_subset,
                   fixed=list(nu=0.5), method="REML",
                   control.HLfit=list(algebra="decorr"))
    }, error=function(err){
      # Return error object if fitting fails
      return(err)
    })

    # Store results
    res_ls[[i]] = list()
    
    # Check for timeout exception
    if(any(class(exp_out) == 'TimeoutException')){
      res_ls[[i]][['spmod']] = 'time_out'
    } else if(any(class(exp_out) == 'simpleError')){
      # Mark as non-convergent
      res_ls[[i]][['spmod']] = 'no_conv'
    } else{
      # Store successful model fit
      res_ls[[i]][['spmod']] = exp_out
    }
    
    res_ls[[i]][['sample']] = sample_tmp
    res_ls[[i]][['gene']] = gene_tmp
  }
  
  return(res_ls)
}


##
#' @title STdiff_check_convergence: Check model convergence and extract summary statistics
#' @description Checks the convergence status of fitted spatial models and extracts
#' summary statistics including p-values, coefficient estimates, and warnings.
#'
#' This function handles multiple convergence states:
#' - Successful fits (HLfit objects)
#' - Timeouts ('time_out')
#' - Non-convergence ('no_convergence')
#' - Unknown errors ('unknown_error')
#'
#' @param sp_models list of spatial models from STdiff_run_spatial()
#' @param meta_dict annotation dictionary mapping coded to original cluster names
#' @param pairwise whether tests are pairwise (TRUE) or reference-based (FALSE)
#' @param pval_adj p-value adjustment method for spatial p-values
#' @return list containing:
#'   - sp_de: data frame with spatial model results (p-values, coefficients, convergence status)
#'   - convergence_summary: summary of convergence status across all models
#' @export
#' @keywords model convergence spaMM summary statistics
#' @examples
#' \dontrun{
#' # After running STdiff_run_spatial, check convergence
#' conv_result = STdiff_check_convergence(sp_models=spatial_res$sp_models,
#'                                        meta_dict=prep$meta_dict,
#'                                        pairwise=FALSE)
#' print(conv_result$convergence_summary)
#' }
STdiff_check_convergence = function(sp_models=NULL, meta_dict=NULL, pairwise=NULL,
                                    pval_adj='fdr'){

  # To prevent NOTES in R CMD check
  . = NULL

  # Initialize results data frame
  sp_de = tibble::tibble()
  
  # Iterate through samples and clusters
  for(i in names(sp_models)){
    for(cl in 1:length(sp_models[[i]])){
      for(test_name in names(sp_models[[i]][[cl]])){
        # Extract model or status
        model_obj = sp_models[[i]][[cl]][[test_name]][['spmod']]
        
        # Determine convergence status and extract statistics
        if(any(class(model_obj) == 'HLfit')){
          # Successful fit - extract summary statistics
          exp_summ = spaMM::summary.HLfit(model_obj, details=c(p_value='Wald'), verbose=F)
          exp_summ = exp_summ[['beta_table']]
          exp_p_val = exp_summ[[2, 4]]
          exp_estimate = exp_summ[[2, 2]]
          exp_std_error = exp_summ[[2, 3]]
          exp_warn = NA_character_
          conv_status = 'converged'
        } else if(model_obj == 'time_out'){
          # Timeout - model did not complete
          exp_p_val = NA_real_
          exp_estimate = NA_real_
          exp_std_error = NA_real_
          exp_warn = 'time_out'
          conv_status = 'timeout'
        } else if(model_obj == 'no_conv'){
          # Non-convergence - model failed to converge
          exp_p_val = NA_real_
          exp_estimate = NA_real_
          exp_std_error = NA_real_
          exp_warn = 'no_convergence'
          conv_status = 'non_convergent'
        } else{
          # Unknown error
          exp_p_val = NA_real_
          exp_estimate = NA_real_
          exp_std_error = NA_real_
          exp_warn = 'unknown_error'
          conv_status = 'error'
        }

        # Extract cluster names from test name
        meta1_exp_tmp = stringr::str_replace(test_name, paste0(i, '_', sp_models[[i]][[cl]][[test_name]][['gene']], '_'), '')
        meta1_exp_tmp = stringr::str_replace(meta1_exp_tmp, '_other$', '')
        
        if(pairwise){
          meta2_exp_tmp = stringr::str_extract(meta1_exp_tmp, 'c[0-9]+$')
          meta1_exp_tmp = stringr::str_extract(meta1_exp_tmp, '^c[0-9]+')
        }

        # Create result tibble
        tibble_tmp = tibble::tibble(
          sample = i,
          gene = sp_models[[i]][[cl]][[test_name]][['gene']],
          cluster_1 = as.vector(unlist(meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == meta1_exp_tmp ])),
          exp_p_val = exp_p_val,
          exp_estimate = exp_estimate,
          exp_std_error = exp_std_error,
          comments_spatial = exp_warn,
          convergence_status = conv_status
        )
        
        if(pairwise){
          tibble_tmp = tibble_tmp %>%
            tibble::add_column(
              cluster_2 = as.vector(unlist(meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == meta2_exp_tmp ])),
              .after = 'cluster_1'
            )
        }
        
        sp_de = dplyr::bind_rows(sp_de, tibble_tmp)
      }
    }
  }

  # Calculate adjusted p-values for converged models
  if(!is.null(sp_de) && nrow(sp_de) > 0){
    # Filter to converged models with valid p-values
    converged_models = sp_de %>%
      dplyr::filter(convergence_status == 'converged', !is.na(exp_p_val))
    
    if(nrow(converged_models) > 0){
      # Calculate adjusted p-values per sample
      for(sample_name in unique(sp_de[['sample']])){
        sample_converged = converged_models %>% dplyr::filter(sample == sample_name)
        adjusted_pvals = p.adjust(sample_converged[['exp_p_val']], method=pval_adj)
        sp_de[[sample_name, 'exp_adj_p_val']] = NA_real_  # Initialize
        sp_de[sp_de[['sample']] == sample_name & sp_de[['convergence_status']] == 'converged', 'exp_adj_p_val'] = adjusted_pvals
      }
    }
  }

  # Generate convergence summary
  if(!is.null(sp_de) && nrow(sp_de) > 0){
    conv_summary = sp_de %>%
      dplyr::group_by(convergence_status) %>%
      dplyr::summarize(
        n_models = dplyr::n(),
        n_genes = dplyr::n_distinct(gene),
        .groups = 'drop'
      ) %>%
      dplyr::arrange(dplyr::desc(n_models))
  } else{
    conv_summary = tibble::tibble(
      convergence_status = 'no_models',
      n_models = 0,
      n_genes = 0
    )
  }

  return(list(sp_de=sp_de, convergence_summary=conv_summary))
}


##
#' @title STdiff_extract_spatial_pvalues: Extract spatial p-values from model results
#' @description Extracts spatial p-values and associated statistics from fitted
#' spatial models, handling convergence failures and generating adjusted p-values.
#'
#' This is a convenience function that wraps STdiff_check_convergence and provides
#' a simplified output focused on p-values for downstream analysis.
#'
#' @param sp_models list of spatial models from STdiff_run_spatial()
#' @param meta_dict annotation dictionary from STdiff_select_genes()
#' @param combo_df combinations dataframe from STdiff_select_genes()
#' @param pairwise whether tests are pairwise (TRUE) or reference-based (FALSE)
#' @param pval_adj p-value adjustment method (default: 'fdr')
#' @return data frame with columns:
#'   - sample: sample name
#'   - gene: gene name
#'   - cluster_1: first cluster (original annotation name)
#'   - cluster_2: second cluster (only if pairwise=TRUE)
#'   - exp_p_val: raw spatial p-value
#'   - exp_adj_p_val: adjusted spatial p-value
#'   - exp_estimate: coefficient estimate for spatial effect
#'   - convergence_status: model convergence status
#'   - comments_spatial: warning/error message if any
#' @export
#' @keywords p-values spatial statistics extraction
#' @examples
#' \dontrun{
#' # Extract spatial p-values for DE gene ranking
#' spatial_pvals = STdiff_extract_spatial_pvalues(
#'   sp_models=spatial_res$sp_models,
#'   meta_dict=prep$meta_dict,
#'   combo_df=prep$combo_df,
#'   pairwise=FALSE,
#'   pval_adj='fdr'
#' )
#' # View top spatial DE genes
#' head(spatial_pvals[order(spatial_pvals$exp_adj_p_val), ], 10)
#' }
STdiff_extract_spatial_pvalues = function(sp_models=NULL, meta_dict=NULL,
                                          combo_df=NULL, pairwise=NULL,
                                          pval_adj='fdr'){

  # Use STdiff_check_convergence to extract all statistics
  conv_result = STdiff_check_convergence(
    sp_models=sp_models,
    meta_dict=meta_dict,
    pairwise=pairwise,
    pval_adj=pval_adj
  )
  
  sp_de = conv_result[['sp_de']]
  
  # Select and reorder columns for clean output
  if(pairwise){
    result_df = sp_de %>%
      dplyr::select(
        sample, gene, cluster_1, cluster_2,
        exp_p_val, exp_adj_p_val, exp_estimate,
        convergence_status, comments_spatial
      )
  } else{
    result_df = sp_de %>%
      dplyr::select(
        sample, gene, cluster_1,
        exp_p_val, exp_adj_p_val, exp_estimate,
        convergence_status, comments_spatial
      )
  }
  
  return(result_df)
}
