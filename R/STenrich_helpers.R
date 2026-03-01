##
# STenrich Helper Functions
# Modular helper functions for spatial enrichment analysis
#
# All functions are marked internal and not exported
#

##
#' @title STenrich_validate_input: Validate input parameters
#' @description Check all input parameters and return validated values
#' @param x an STlist with transformed gene expression
#' @param samples a vector with sample names or indexes to run analysis
#' @param gene_sets a named list of gene sets to test
#' @param score_type Controls how gene set expression is calculated. Options: 'avg' or 'gsva'
#' @param annot name of the annotation within `x@spatial_meta` containing spot/cell categories
#' @param domain the domain to restrict the analysis
#' @param num_sds number of standard deviations to set minimum gene set expression threshold (default: 1)
#' @param min_units Minimum number of spots with high expression (default: 20)
#' @param min_genes the minimum number of genes of a gene set present in data (default: 5)
#' @param pval_adj_method the method for multiple comparison adjustment (default: 'BH')
#' @return list with validated parameters
#' @keywords internal
STenrich_validate_input = function(x, samples, gene_sets, score_type, annot, domain, num_sds, min_units, min_genes, pval_adj_method){
  
  # Validate STlist
  if(is.null(x) | !is(x, 'STlist')){
    stop("The input must be a STlist.")
  }
  
  # Validate samples
  if(is.null(samples)){
    samples = names(x@spatial_meta)
  } else{
    if(is.numeric(samples)){
      samples = as.vector(na.omit(names(x@spatial_meta)[samples]))
    } else{
      samples = samples[samples %in% names(x@spatial_meta)]
    }
    if(length(samples) == 0 | !any(samples %in% names(x@spatial_meta))){
      raise_err(err_code="error0020")
    }
  }
  
  # Validate score_type
  score_type = tolower(score_type)
  if(!(score_type %in% c('avg', 'gsva'))){
    warning('Only `avg` or `gsva` are allowed in `score_type`. Using `avg`.')
    score_type = 'avg'
  }
  
  # Type conversions
  num_sds = as.double(num_sds)
  min_units = as.integer(min_units)
  min_genes = as.integer(min_genes)
  
  return(list(x=x, samples=samples, gene_sets=gene_sets, score_type=score_type,
              annot=annot, domain=domain, num_sds=num_sds, min_units=min_units, 
              min_genes=min_genes, pval_adj_method=pval_adj_method))
}


##
#' @title STenrich_prepare_data: Prepare data for enrichment analysis
#' @description Extract tissue spots, coordinates, and gene sets from STlist
#' @param x an STlist with transformed gene expression
#' @param samples a vector with sample names to run analysis
#' @param gene_sets a named list of gene sets to test
#' @param annot name of the annotation within `x@spatial_meta` containing spot/cell categories
#' @param domain the domain to restrict the analysis
#' @param min_units Minimum number of spots with high expression (default: 20)
#' @return list with prepared data including tissue_spots, coords_df, combo, pw_genes
#' @keywords internal
STenrich_prepare_data = function(x, samples, gene_sets, annot, domain, min_units){
  
  # Domain filtering
  tissue_spots = NULL
  sample_rm = c()
  if(!is.null(annot) & !is.null(domain)){
    for(i in samples){
      if(any(colnames(x@spatial_meta[[i]]) == annot)){
        if(any(x@spatial_meta[[i]][[annot]] %in% domain)){
          tissue_spots_tmp = x@spatial_meta[[i]][[1]][x@spatial_meta[[i]][[annot]] %in% domain]
          if(length(tissue_spots_tmp) >= min_units){
            tissue_spots[[i]] = tissue_spots_tmp
          } else{
            message(paste0('\tSample ', i, ' has less than two spots/cells assigned to domain. Skipping.\n'))
            sample_rm = append(sample_rm, i)
          }
        } else{
          message(paste0('\tSample ', i, ' does not contain the specified domains. Skipping.\n'))
          sample_rm = append(sample_rm, i)
        }
      } else{
        message(paste0('\tSample ', i, ' does not contain the specified annotation. Skipping.\n'))
        sample_rm = append(sample_rm, i)
      }
    }
  }
  samples = samples[!(samples %in% sample_rm)]
  
  if(length(samples) < 1){
    raise_err(err_code="error0021")
  }
  
  # Extract coordinates
  coords_df = lapply(samples, function(i){
    df_tmp = x@spatial_meta[[i]][, c('libname', 'xpos', 'ypos')]
    if(length(tissue_spots) > 0){
      df_tmp = df_tmp[df_tmp[['libname']] %in% tissue_spots[[i]], ]
    }
    df_tmp = df_tmp %>% tibble::column_to_rownames(var='libname')
    as.matrix(df_tmp)
  })
  names(coords_df) = samples
  
  # Match gene sets to available genes
  combo = expand.grid(sample_name=samples, gene_set=names(gene_sets))
  pw_genes = lapply(1:nrow(combo), function(i){
    sample_tmp = as.vector(combo[i, 1])
    gs_tmp = as.vector(combo[i, 2])
    pw_genes_tmp = unique(gene_sets[[gs_tmp]])
    intersect(pw_genes_tmp, rownames(x@tr_counts[[sample_tmp]]))
  })
  
  return(list(x=x, samples=samples, gene_sets=gene_sets, tissue_spots=tissue_spots,
              coords_df=coords_df, combo=combo, pw_genes=pw_genes))
}


# Note: STenrich_calculate_expression removed - calling functions directly in main STenrich


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
  
  # Loop through combinations of samples and gene sets
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
  
  # Create dataframes per sample
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
STenrich_permutation_test = function(result_df, coords_df, pw_genes, samples, gene_sets, num_sds, min_units, reps, seed, cores, verbose){
  
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
      
      # Calculate expression threshold
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
  
  # Calculate gene set coverage
  pval_res = pval_res %>%
    tibble::add_column(prop_size_test=.[['size_test']]/.[['size_gene_set']], .before='p_value') %>%
    dplyr::mutate(prop_size_test=round(prop_size_test, 3))
  
  # Split results among samples and adjust p-values
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