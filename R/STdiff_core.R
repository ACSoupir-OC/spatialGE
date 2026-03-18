##
# STdiff Core Modular Functions
# Refactored differential expression analysis functions

##
#' @title STdiff_select_genes: Gene filtering and non-spatial differential testing
#' @description Selects variable genes and performs non-spatial differential expression
#' testing (linear models, t-tests, or Wilcoxon tests)
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
#' @param pairwise whether to perform pairwise tests
#' @param verbose verbosity level
#' @param cores number of cores for parallelization
#' @return list containing: combo_df, meta_dict, non_spatial_results, pval_thr, test_type
#' @keywords internal
STdiff_select_genes = function(x=NULL, samples=NULL, annot=NULL, w=NULL, k=NULL, deepSplit=NULL,
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
        res_tmp = non_spatial_de(expr_data=expr_df, combo=combo_clust_cl, pairwise=pairwise)
      } else if(test_type == 't_test' | test_type == 'wilcoxon'){
        res_tmp = stdiff_mean_test(expr_data=expr_df, combo=combo_clust_cl, pairwise=pairwise, test_type=test_type)
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
      non_sp_de_tmp[[i]][['adj_p_val']] = p.adjust(non_sp_de_tmp[[i]][['mm_p_val']], method=pval_adj)
    } else if(test_type == 't_test'){
      non_sp_de_tmp[[i]][['adj_p_val']] = p.adjust(non_sp_de_tmp[[i]][['ttest_p_val']], method=pval_adj)
    } else if(test_type == 'wilcoxon'){
      non_sp_de_tmp[[i]][['adj_p_val']] = p.adjust(non_sp_de_tmp[[i]][['wilcox_p_val']], method=pval_adj)
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
#' @title STdiff_fit_spatial: Fit spatial mixed models for selected genes
#' @description Fits spatial mixed models with Matern covariance structure for genes
#' selected from non-spatial testing
#' @param x an STlist object
#' @param prep list from STdiff_select_genes containing combo_df, meta_dict, non_spatial_results
#' @param sp_topgenes proportion of top DE genes to fit spatial models on
#' @param cores number of cores for parallelization
#' @param verbose verbosity level
#' @return list containing spatial model results
#' @keywords internal
STdiff_fit_spatial = function(x=NULL, prep=NULL, sp_topgenes=0.2, cores=NULL, verbose=1L){

  # If sp_topgenes is 0, return empty
  if(sp_topgenes <= 0){
    return(list())
  }

  if(verbose){
    cat('\tRunning spatial tests...\n')
  }

  non_sp_de_tmp = prep[['non_spatial_results']]
  meta_dict = prep[['meta_dict']]
  combo_df = prep[['combo_df']]
  sp_topgenes = as.double(sp_topgenes)
  
  if(is.null(cores)){
    cores = count_cores(length(non_sp_de_tmp))
  } else{
    cores = ceiling(cores)
  }

  sp_models = list()
  for(sample_name in names(non_sp_de_tmp)){
    # Subset genes based on adjusted p-values
    models_keep = non_sp_de_tmp[[sample_name]] %>%
      dplyr::filter(adj_p_val < prep[['pval_thr']]) %>%
      dplyr::arrange(adj_p_val)
    
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

    # Mark genes for spatial testing
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

    # Recode annotations
    for(spotrow in 1:nrow(models_keep)){
      models_keep[['cluster_1']][spotrow] = meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == models_keep[['cluster_1']][spotrow] ]
      if(prep[['pairwise']]){
        models_keep[['cluster_2']][spotrow] = meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == models_keep[['cluster_2']][spotrow] ]
      }
    }

    if(nrow(models_keep) > 0){
      # Create gene x cluster combinations
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

      # Create expression data frame
      expr_df = expandSparse(x@tr_counts[[sample_name]])
      expr_df = t(expr_df) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var='libname') %>%
        dplyr::right_join(x@spatial_meta[[sample_name]] %>%
                            tibble::add_column(group=1, .after='libname') %>%
                            dplyr::select(c('libname', 'group', 'ypos', 'xpos')), . , by='libname') %>%
        dplyr::select(c('libname', 'group', 'ypos', 'xpos'))

      # Run spatial models in parallel
      sp_models[[sample_name]] = parallel::mclapply(1:length(clusters_tmp), function(i){
        models_keep_tmp = models_keep %>%
          dplyr::filter(cluster_1 == clusters_tmp[i]) %>%
          dplyr::select('genexcluster') %>%
          unlist() %>% as.vector()
        non_sp_models_tmp = non_sp_de_tmp[[sample_name]][ models_keep_tmp ]

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
        sp_mods = spatial_de(expr_dat=expr_df, non_sp_mods=non_sp_models_tmp, annot_dict=meta_dict, verb=verbose)
        return(sp_mods)
      }, mc.cores=cores)
    } else{
      warning(paste0('No genes to test for ', sample_name, '. Try increasing sp_topgenes.'))
    }
  }

  return(list(sp_models=sp_models, meta_dict=meta_dict, combo_df=combo_df))
}


##
#' @title STdiff: Differential expression analysis (modular wrapper)
#' @description Public wrapper function that calls the modular STdiff functions
#' @details Maintains the original API while using the new modular implementation
#' @inheritParams STdiff
#' @export
STdiff = function(x=NULL, samples=NULL, annot=NULL, w=NULL, k=NULL, deepSplit=NULL,
                  topgenes=5000, pval_thr=0.05, pval_adj='fdr', test_type='mm', sp_topgenes=0.2,
                  clusters=NULL, pairwise=FALSE, verbose=1L, cores=NULL){

  # Record time
  zero_t = Sys.time()

  # Step 1: Select genes and run non-spatial tests
  if(verbose){
    cat('Step 1: Selecting genes and running non-spatial tests...\n')
  }
  prep = STdiff_select_genes(x=x, samples=samples, annot=annot, w=w, k=k, deepSplit=deepSplit,
                             topgenes=topgenes, pval_thr=pval_thr, pval_adj=pval_adj, test_type=test_type,
                             clusters=clusters, pairwise=pairwise, verbose=verbose, cores=cores)

  # Step 2: Fit spatial models if requested
  spatial_res = NULL
  if(test_type == 'mm' && sp_topgenes > 0){
    if(verbose){
      cat('Step 2: Fitting spatial mixed models...\n')
    }
    spatial_res = STdiff_fit_spatial(x=x, prep=prep, sp_topgenes=sp_topgenes, cores=cores, verbose=verbose)
  }

  # Step 3: Compile results
  if(verbose){
    cat('Step 3: Compiling results...\n')
  }
  result_de = STdiff_compile_results(prep=prep, spatial=spatial_res)

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('STdiff completed in ', round(end_t, 2), ' min.\n'))
  }

  return(result_de)
}