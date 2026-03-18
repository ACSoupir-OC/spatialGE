##
# STdiff Non-Spatial Module
# 
# This file contains non-spatial differential expression testing functions
# extracted from STdiff.R for better modularity and maintainability.
#
# @importFrom magrittr %>%
# @importFrom rlang :=
# @importFrom stats p.adjust t.test wilcox.test
# @importFrom tibble tibble bind_rows
# @importFrom dplyr filter mutate select arrange bind_rows case_when
# @importFrom dplyr rename left_join full_join right_join slice_head group_by ungroup
# @importFrom dplyr relocate relocate column_to_rownames
# @importFrom spaMM fitme summary.HLfit

##
#' @title STdiff_run_nonspatial: Run non-spatial differential expression tests
#' @description Main entry point for non-spatial DE testing. Runs linear models,
#' t-tests, or Wilcoxon tests for selected genes between groups of spots/cells
#' @param x an STlist object
#' @param samples vector of sample names to test
#' @param annot column name in spatial_meta for cluster annotations
#' @param w spatial weight parameter (used if annot is NULL)
#' @param k number of clusters (used if annot is NULL)
#' @param deepSplit deepSplit parameter for dynamicTreeCut clusters
#' @param topgenes number of top variable genes to select
#' @param pval_thr p-value threshold for selecting DE genes
#' @param pval_adj p-value adjustment method
#' @param test_type type of test: 'mm', 't_test', or 'wilcoxon'
#' @param clusters optional vector of specific clusters to test
#' @param pairwise whether to perform pairwise tests (TRUE) or reference-based (FALSE)
#' @param verbose verbosity level (0, 1, or 2)
#' @param cores number of cores for parallelization
#' @return list containing: combo_df, meta_dict, non_spatial_results, pval_thr, test_type, pairwise
#' @export
#' @keywords non-spatial differential expression
STdiff_run_nonspatial = function(x=NULL, samples=NULL, annot=NULL, w=NULL, k=NULL, deepSplit=NULL,
                                 topgenes=5000, pval_thr=0.05, pval_adj='fdr', test_type='mm',
                                 clusters=NULL, pairwise=FALSE, verbose=1L, cores=NULL){

  # To prevent NOTES in R CMD check
  . = NULL

  # Convert user inputs to expected types
  topgenes = as.integer(ceiling(topgenes))
  pval_thr = as.double(pval_thr)
  verbose = as.integer(verbose)
  if(!is.integer(verbose) | !(verbose %in% c(0L, 1L, 2L))){
    verbose = 1L
  }

  # Check pairwise requirement
  if(pairwise & !is.null(clusters)){
    if(length(clusters) < 2){
      stop('If pairwise tests requested, at least two clusters are required.')
    }
  }

  # Validate topgenes and pval_thr
  if((round(topgenes, 0) <= 0) | (pval_thr < 0 | pval_thr > 1)){
    stop('topgenes or pval_thr contain invalid values.')
  }

  # Define samples
  if(is.null(samples)){
    samples = names(x@spatial_meta)
  } else{
    if(is.numeric(samples)){
      samples = as.vector(na.omit(names(x@spatial_meta)[samples]))
    } else{
      samples = samples[samples %in% names(x@spatial_meta)]
    }
    if(length(samples) == 0 | !any(samples %in% names(x@spatial_meta))){
      stop('None of the requested samples are present in the STlist.')
    }
  }

  # Define annotation column
  if(is.null(annot)){
    if(!is.null(w) & !is.null(k)){
      if(length(w) != 1 | length(k) != 1){
        stop('Please specify a single value for w and a single value for k.')
      }
      annot = paste0('stclust_spw', as.character(w))
      if(k == 'dtc'){
        if(is.null(deepSplit)){
          stop('If k="dtc", then specify deepSplit.')
        } else if(is.logical(deepSplit)){
          annot = paste0(annot, '_dspl', ifelse(deepSplit, 'True', 'False'))
        } else{
          annot = paste0(annot, '_dspl', deepSplit)
        }
      } else{
        if(!is.numeric(k)){
          stop('Specify a valid k value.')
        }
        annot = paste0(annot, '_k', as.character(k))
      }
    } else{
      stop('If no specific annotation is specified, please specify both w and k (STclust parameters).')
    }
  } else{
    if(length(annot) == 0){
      stop('If w and k are not specified, one annotation column from @spatial_meta should be specified.')
    } else if(length(annot) > 1){
      stop('Only one annotation column from @spatial_meta can be tested at a time.')
    }
  }

  # Validate annotation exists in samples
  samples_tmp = samples
  for(i in samples){
    if(!(annot %in% colnames(x@spatial_meta[[i]]))){
      samples_tmp = grep(i, samples_tmp, value=T, invert=T)
      if(verbose){
        cat(paste0('Skipping ', i, '. Annotation not available for this sample.\n'))
      }
    }
    if(length(samples_tmp) == 0){
      stop('No samples left to test. Are the requested annotations/clusters present in at least one sample?')
    }
  }
  samples = samples_tmp

  # Extract annotations for each sample
  metas = tibble::tibble()
  for(sample_name in samples){
    meta_tmp = as.character(unique(x@spatial_meta[[sample_name]][[annot]]))
    metas = dplyr::bind_rows(metas, tibble::tibble(samplename=sample_name, orig_annot=meta_tmp))
  }

  # Find variable genes for each sample
  genes = tibble::tibble()
  for(sample_name in samples){
    genes_tmp = x@gene_meta[[sample_name]] %>% dplyr::arrange(dplyr::desc(vst.variance.standardized))
    if(!is.null(topgenes)){
      genes_tmp = genes_tmp %>% dplyr::slice_head(n=topgenes)
    }
    genes_tmp = as.vector(genes_tmp[['gene']])
    genes = dplyr::bind_rows(genes, tibble::tibble(samplename=sample_name, gene=genes_tmp))
  }

  # Merge sample names, genes, and annotations
  gene_and_meta = dplyr::full_join(genes, metas, by='samplename', multiple='all')

  # Create annotation dictionary
  meta_dict = tibble::tibble(orig_annot=unique(gene_and_meta[['orig_annot']]),
                             coded_annot=paste0('c', seq(1, length(unique(gene_and_meta[['orig_annot']])))))

  gene_and_meta = gene_and_meta %>% dplyr::left_join(., meta_dict, by='orig_annot') %>% dplyr::rename(meta=coded_annot)

  # Create combinations to test
  combo_df = prepare_stdiff_combo(to_expand=gene_and_meta, user_clusters=clusters, pairwise=pairwise, verbose=verbose)

  # Check test type
  if(test_type == 'mm'){
    test_print = 'non-spatial mixed models'
  } else if(test_type == 't_test'){
    test_print = 'T-tests'
  } else if(test_type == 'wilcoxon'){
    test_print = "Wilcoxon's tests"
  } else{
    stop('Test type is one of "mm", "t_test", or "wilcoxon".')
  }

  if(verbose){
    cat(paste0('\tRunning ', test_print, '...\n'))
  }

  # Subset combo_df to specified clusters if provided
  if(!is.null(clusters)){
    coded_metas = meta_dict[['coded_annot']][ meta_dict[['orig_annot']] %in% clusters ]
    combo_df = combo_df[ combo_df[['meta1']] %in% coded_metas, ]
    if(pairwise){
      combo_df = combo_df[ combo_df[['meta2']] %in% coded_metas, ]
    }
    if(nrow(combo_df) < 1){
      stop('All cluster x gene comparisons were removed. Is at least one of the specified clusters present in one of the samples?')
    }
  }

  # Set number of cores
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(unique(combo_df[['samplename']])))
  } else{
    cores = ceiling(cores)
  }

  # Run non-spatial tests
  non_sp_models = parallel::mclapply(1:length(unique(combo_df[['samplename']])), function(i){
    sample_name = unique(combo_df[['samplename']])[i]
    combo_clust = combo_df[ combo_df[['samplename']] == sample_name, ]

    # Create expression data frame
    expr_df = expandSparse(x@tr_counts[[sample_name]])
    expr_df = expr_df[ rownames(expr_df) %in% unique(combo_clust[['gene']]), ]
    expr_df = t(expr_df) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var='libname') %>%
      dplyr::right_join(x@spatial_meta[[sample_name]] %>%
                          tibble::add_column(group=1, .after='libname') %>%
                          dplyr::select(c('libname', 'group', 'ypos', 'xpos'), orig_annot:=!!annot),. , by='libname') %>%
      dplyr::left_join(., meta_dict, by='orig_annot') %>%
      dplyr::rename(meta=coded_annot) %>%
      tibble::column_to_rownames(var='libname') %>%
      dplyr::relocate(meta, .before=1) %>%
      dplyr::select(-c('orig_annot'))

    # Run tests
    res = list()
    for(cluster_tmp in unique(combo_clust[['meta1']])){
      combo_clust_cl = combo_clust[ combo_clust[['meta1']] == cluster_tmp, ]
      
      if(verbose >= 1L && !pairwise){
        system(sprintf('echo "%s"', paste0("\t\tSample: ", sample_name, ", ",
                                           meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == cluster_tmp ],
                                           " (", nrow(combo_clust_cl), " tests)")))
      }
      
      if(test_type == 'mm'){
        res_tmp = STdiff_fit_mm(expr_data=expr_df, combo=combo_clust_cl, pairwise=pairwise)
      } else if(test_type == 't_test'){
        res_tmp = STdiff_fit_ttest(expr_data=expr_df, combo=combo_clust_cl, pairwise=pairwise)
      } else if(test_type == 'wilcoxon'){
        res_tmp = STdiff_fit_wilcoxon(expr_data=expr_df, combo=combo_clust_cl, pairwise=pairwise)
      }
      res = append(res, res_tmp)
    }
    return(res)
  }, mc.cores=cores)

  names(non_sp_models) = unique(combo_df[['samplename']])

  # Compile non-spatial results
  result_de = tibble::tibble()
  for(i in names(non_sp_models)){
    for(mod in names(non_sp_models[[i]])){
      sample_tmp = non_sp_models[[i]][[mod]][['samplename']]
      gene_tmp = non_sp_models[[i]][[mod]][['gene']]
      avglfc_tmp = non_sp_models[[i]][[mod]][['avglfc']]
      cluster_1_tmp = as.vector(unlist(meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == non_sp_models[[i]][[mod]][['meta1']] ]))

      df_tmp = tibble::tibble(sample=sample_tmp,
                              gene=gene_tmp,
                              avg_log2fc=avglfc_tmp,
                              cluster_1=cluster_1_tmp)

      if(pairwise){
        cluster_2_tmp = as.vector(unlist(meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == non_sp_models[[i]][[mod]][['meta2']] ]))
        df_tmp[['cluster_2']] = cluster_2_tmp
      }

      if(test_type == 'mm'){
        df_tmp[['mm_p_val']] = non_sp_models[[i]][[mod]][['nonspmod']]
      } else if(test_type == 't_test'){
        df_tmp[['ttest_p_val']] = non_sp_models[[i]][[mod]][['mean_test']][['p.value']]
      } else if(test_type == 'wilcoxon'){
        df_tmp[['wilcox_p_val']] = non_sp_models[[i]][[mod]][['mean_test']][['p.value']]
      }

      result_de = dplyr::bind_rows(result_de, df_tmp)
    }
  }

  # Adjust p-values
  non_sp_de_tmp = list()
  for(i in unique(result_de[['sample']])){
    non_sp_de_tmp[[i]] = result_de %>% dplyr::filter(sample == i)
    if(test_type == 'mm'){
      non_sp_de_tmp[[i]][['adj_p_val']] = STdiff_adjust_pvalues(pvals=non_sp_de_tmp[[i]][['mm_p_val']], method=pval_adj)
    } else if(test_type == 't_test'){
      non_sp_de_tmp[[i]][['adj_p_val']] = STdiff_adjust_pvalues(pvals=non_sp_de_tmp[[i]][['ttest_p_val']], method=pval_adj)
    } else if(test_type == 'wilcoxon'){
      non_sp_de_tmp[[i]][['adj_p_val']] = STdiff_adjust_pvalues(pvals=non_sp_de_tmp[[i]][['wilcox_p_val']], method=pval_adj)
    }
    
    if(pairwise){
      non_sp_de_tmp[[i]] = non_sp_de_tmp[[i]] %>% dplyr::arrange(cluster_1, cluster_2, adj_p_val, dplyr::desc(avg_log2fc))
    } else{
      non_sp_de_tmp[[i]] = non_sp_de_tmp[[i]] %>% dplyr::arrange(cluster_1, adj_p_val, dplyr::desc(avg_log2fc))
    }
  }

  return(list(combo_df=combo_df,
              meta_dict=meta_dict,
              non_spatial_results=non_sp_de_tmp,
              pval_thr=pval_thr,
              test_type=test_type,
              pairwise=pairwise))
}


##
#' @title STdiff_fit_mm: Fit non-spatial mixed models
#' @description Fits non-spatial linear mixed models using spaMM for testing
#' differential expression between groups of spots/cells
#' @param expr_data a data frame (rows=spots/cells) with normalized expression data
#' @param combo a data frame with columns (samplename, meta1, meta2, gene) containing
#' all combinations of gene x cluster to test
#' @param pairwise whether to perform pairwise tests (TRUE) or reference-based (FALSE)
#' @return a list containing non-spatial model results with p-values, log-fold changes,
#' and metadata
#' @export
#' @keywords non-spatial mixed models
STdiff_fit_mm = function(expr_data=NULL, combo=NULL, pairwise=NULL){
  models_ls = list()
  for(i in 1:nrow(combo)){
    sample_tmp = combo[i, 1]
    meta1_tmp = combo[i, 2]
    meta2_tmp = combo[i, 3]
    gene_tmp = combo[i, 4]

    # Filter data if pairwise or recode annotations if not pairwise
    if(pairwise){
      expr_tmp = expr_data %>%
        dplyr::filter(meta %in% c(meta1_tmp, meta2_tmp))
    } else{
      expr_tmp = expr_data %>%
        dplyr::mutate(meta=factor(dplyr::case_when(meta != meta1_tmp ~ 'other', TRUE ~ meta),
                                  levels=c('other', meta1_tmp)))
    }

    # Get relevant gene column
    expr_tmp = expr_tmp %>%
      dplyr::select(c('group', 'ypos', 'xpos', 'meta'), exprval:=!!gene_tmp)

    # Calculate log-fold change
    avgexpr_cl1 = mean(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select(exprval) %>% unlist())
    avgexpr_cl2 = mean(expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select(exprval) %>% unlist())
    avglogfold_tmp = (avgexpr_cl1 - avgexpr_cl2)

    # Create non-spatial model
    error_message = c()
    warning_message = c()
    res_mod = withCallingHandlers(
      tryCatch(spaMM::fitme(formula=stats::as.formula('exprval~meta'), data=expr_tmp, method="REML"),
               error=function(e){error_message <<- conditionMessage(e)}),
      warning=function(w){warning_message <<- conditionMessage(w)
      invokeRestart("muffleWarning")}
    )

    # Calculate p-value
    res_mod = spaMM::summary.HLfit(res_mod, details=c(p_value="Wald"), verbose=F)[['beta_table']][2, 4]

    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]] = list()
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['nonspmod']] = res_mod
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['samplename']] = sample_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['gene']] = gene_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['avglfc']] = avglogfold_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta1']] = meta1_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta2']] = meta2_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['warning']] = warning_message
  }
  return(models_ls)
}


##
#' @title STdiff_fit_ttest: Fit non-spatial t-tests
#' @description Performs Welch's t-tests for testing differential expression
#' between groups of spots/cells
#' @param expr_data a data frame with expression data and metadata
#' @param combo a data frame with columns (samplename, meta1, meta2, gene)
#' @param pairwise whether to perform pairwise tests (TRUE) or reference-based (FALSE)
#' @return a list containing test results with p-values, log-fold changes, and metadata
#' @export
#' @keywords non-spatial t-tests
STdiff_fit_ttest = function(expr_data=NULL, combo=NULL, pairwise=NULL){
  test_ls = list()
  for(i in 1:nrow(combo)){
    sample_tmp = combo[i, 1]
    meta1_tmp = combo[i, 2]
    meta2_tmp = combo[i, 3]
    gene_tmp = combo[i, 4]

    # Filter data if pairwise or recode annotations if not pairwise
    if(pairwise){
      expr_tmp = expr_data %>%
        dplyr::filter(meta %in% c(meta1_tmp, meta2_tmp))
    } else{
      expr_tmp = expr_data %>%
        dplyr::mutate(meta=factor(dplyr::case_when(meta != meta1_tmp ~ 'other', TRUE ~ meta),
                                  levels=c('other', meta1_tmp)))
    }

    # Get relevant gene column
    expr_tmp = expr_tmp %>%
      dplyr::select(c('group', 'ypos', 'xpos', 'meta'), exprval:=!!gene_tmp)

    # Calculate log-fold change
    avgexpr_cl1 = mean(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select(exprval) %>% unlist())
    avgexpr_cl2 = mean(expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select(exprval) %>% unlist())
    avglogfold_tmp = (avgexpr_cl1 - avgexpr_cl2)

    # Perform t-test for differences in mean
    res_test = stats::t.test(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select('exprval') %>% unlist(),
                             expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select('exprval') %>% unlist())

    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]] = list()
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['mean_test']] = res_test
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['samplename']] = sample_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['gene']] = gene_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['avglfc']] = avglogfold_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta1']] = meta1_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta2']] = meta2_tmp
  }
  return(test_ls)
}


##
#' @title STdiff_fit_wilcoxon: Fit non-spatial Wilcoxon tests
#' @description Performs Wilcoxon rank-sum tests for testing differential expression
#' between groups of spots/cells (non-parametric alternative to t-test)
#' @param expr_data a data frame with expression data and metadata
#' @param combo a data frame with columns (samplename, meta1, meta2, gene)
#' @param pairwise whether to perform pairwise tests (TRUE) or reference-based (FALSE)
#' @return a list containing test results with p-values, log-fold changes, and metadata
#' @export
#' @keywords non-spatial Wilcoxon tests
STdiff_fit_wilcoxon = function(expr_data=NULL, combo=NULL, pairwise=NULL){
  test_ls = list()
  for(i in 1:nrow(combo)){
    sample_tmp = combo[i, 1]
    meta1_tmp = combo[i, 2]
    meta2_tmp = combo[i, 3]
    gene_tmp = combo[i, 4]

    # Filter data if pairwise or recode annotations if not pairwise
    if(pairwise){
      expr_tmp = expr_data %>%
        dplyr::filter(meta %in% c(meta1_tmp, meta2_tmp))
    } else{
      expr_tmp = expr_data %>%
        dplyr::mutate(meta=factor(dplyr::case_when(meta != meta1_tmp ~ 'other', TRUE ~ meta),
                                  levels=c('other', meta1_tmp)))
    }

    # Get relevant gene column
    expr_tmp = expr_tmp %>%
      dplyr::select(c('group', 'ypos', 'xpos', 'meta'), exprval:=!!gene_tmp)

    # Calculate log-fold change
    avgexpr_cl1 = mean(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select(exprval) %>% unlist())
    avgexpr_cl2 = mean(expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select(exprval) %>% unlist())
    avglogfold_tmp = (avgexpr_cl1 - avgexpr_cl2)

    # Perform Wilcoxon test for differences in mean
    res_test = stats::wilcox.test(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select('exprval') %>% unlist(),
                                  expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select('exprval') %>% unlist())

    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]] = list()
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['mean_test']] = res_test
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['samplename']] = sample_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['gene']] = gene_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['avglfc']] = avglogfold_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta1']] = meta1_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta2']] = meta2_tmp
  }
  return(test_ls)
}


##
#' @title STdiff_adjust_pvalues: Adjust p-values using multiple testing correction
#' @description Applies multiple testing correction to p-values using various methods
#' @param pvals numeric vector of p-values to adjust
#' @param method adjustment method (default: 'fdr' for Benjamini-Hochberg)
#' @return numeric vector of adjusted p-values
#' @export
#' @keywords p-value adjustment multiple testing
STdiff_adjust_pvalues = function(pvals=NULL, method='fdr'){
  adjusted_pvals = stats::p.adjust(pvals, method=method)
  return(adjusted_pvals)
}
