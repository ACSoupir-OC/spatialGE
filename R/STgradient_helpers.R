##
# STgradient_helpers.R - Helper functions for STgradient spatial gradient analysis
#

##
#' @title STgradient_validate_input: Validate STgradient inputs
#' @description Validate inputs for STgradient function
#' @param x STlist object
#' @param samples samples to process
#' @param topgenes number of top genes
#' @param annot annotation column
#' @param ref reference cluster
#' @param exclude optional exclude cluster
#' @param min_nb minimum neighbors
#' @param nb_dist_thr neighborhood distance threshold
#' @param cores number of cores
#' @param verbose verbosity
#' @return validated parameters list
#' @keywords internal
STgradient_validate_input = function(x, samples, topgenes, annot, ref, exclude, min_nb, nb_dist_thr, cores, verbose){
  # Define samples
  if(is.null(samples)){
    samples = names(x@spatial_meta)
  } else{
    if(is.numeric(samples)){
      samples = as.vector(na.omit(names(x@spatial_meta)[samples]))
    } else{
      samples = samples[samples %in% names(x@spatial_meta)]
    }
    if(length(samples) == 0 || !any(samples %in% names(x@spatial_meta))){
      stop('None of the requested samples are present in the STlist.')
    }
  }
  
  # Remove samples without annotation
  sample_rm = c()
  for(i in samples){
    if(!(annot %in% colnames(x@spatial_meta[[i]]))){
      sample_rm = append(sample_rm, i)
      message(paste0('The annotation specified in `annot` is not present in ', i, '\n'))
    }
  }
  samples = samples[!(samples %in% sample_rm)]
  
  # Remove samples without reference cluster
  sample_rm = c()
  for(i in samples){
    if(sum(x@spatial_meta[[i]][[annot]] == ref) < 1){
      sample_rm = append(sample_rm, i)
    }
  }
  samples = samples[!(samples %in% sample_rm)]
  
  # Define cores
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = parallel::detectCores()
  }
  cores = ceiling(cores)
  
  # Define neighborhood tolerance
  if(is.null(nb_dist_thr) || !is.numeric(nb_dist_thr) || length(nb_dist_thr) != 2){
    nb_dist_thr = c(0.75, 1.25)
    if(x@misc[['platform']] != 'visium'){
      nb_dist_thr = c(0.25, 3)
    }
  }
  
  return(list(samples = samples, annot = annot, ref = ref, exclude = exclude, 
              min_nb = min_nb, nb_dist_thr = nb_dist_thr, cores = cores, verbose = verbose))
}

##

##
#' @title STgradient_prepare_distances: Prepare distance data
#' @description Calculate and categorize distances for gradient analysis
#' @param x STlist object
#' @param sample_name name of sample to process
#' @param annot annotation column name
#' @param ref reference cluster
#' @param exclude optional cluster to exclude
#' @return list with distance data frames
#' @keywords internal
STgradient_prepare_distances = function(x, sample_name, annot, ref, exclude){
  # Calculate euclidean distances
  coords_tmp = x@spatial_meta[[sample_name]][, c('libname', 'ypos', 'xpos'), drop=FALSE]
  dist_tmp = as.matrix(stats::dist(coords_tmp[, c('ypos', 'xpos')], method='euclidean'))
  # Set rownames and colnames to libname
  rownames(dist_tmp) = coords_tmp$libname
  colnames(dist_tmp) = coords_tmp$libname
  
  # Categorize spots
  ref_tmp = x@spatial_meta[[sample_name]]$libname[x@spatial_meta[[sample_name]][[annot]] == ref]
  ref_tmp = ref_tmp[!is.na(ref_tmp)]  # Remove NA values
  nonref_tmp = x@spatial_meta[[sample_name]]$libname[!x@spatial_meta[[sample_name]][[annot]] %in% c(ref, exclude)]
  nonref_tmp = nonref_tmp[!is.na(nonref_tmp)]  # Remove NA values
  
  return(list(dist_tmp = dist_tmp, ref_tmp = ref_tmp, nonref_tmp = nonref_tmp))
}

##
#' @title STgradient_filter_neighbors: Filter reference spots by neighbor count
#' @description Identify reference spots with sufficient neighbors
#' @param dist_tmp distance matrix
#' @param ref_tmp reference barcodes
#' @param min_nb minimum neighbors required
#' @param nb_dist_thr neighborhood distance threshold
#' @return vector of reference barcodes to keep
#' @keywords internal
STgradient_filter_neighbors = function(dist_tmp, ref_tmp, min_nb, nb_dist_thr){
  # Check if ref_tmp is empty
  if(length(ref_tmp) == 0){
    return(c())
  }
  
  # Get minimum distance among all spots
  min_sample = min(as.data.frame(dist_tmp[lower.tri(dist_tmp)]))
  
  # Get distances among reference spots
  dists_ref_tmp = dist_tmp[ref_tmp, ref_tmp, drop=FALSE]
  
  # Get number of neighbors within minimum distance
  nbs = colSums(dists_ref_tmp >= min_sample * nb_dist_thr[1] & dists_ref_tmp <= min_sample * nb_dist_thr[2])
  
  if(sum(nbs >= min_nb) < 1){
    nbs_keep = c()
  } else{
    nbs_keep = names(nbs)[nbs >= min_nb]
  }
  
  return(nbs_keep)
}

##
#' @title STgradient_summarize_distances: Summarize distances from reference
#' @description Calculate min or average distance from reference for non-reference spots
#' @param dist_tmp distance matrix
#' @param nonref_tmp non-reference barcodes
#' @param ref_tmp reference barcodes
#' @param nbs_keep reference barcodes to keep
#' @param distsumm distance summary metric ('min' or 'avg')
#' @param limit optional distance limit
#' @return data frame with barcode and dist2ref
#' @keywords internal
STgradient_summarize_distances = function(dist_tmp, nonref_tmp, ref_tmp, nbs_keep, distsumm, limit){
  # Select spots in analysis
  dists_nonref_tmp = as.data.frame(dist_tmp[nonref_tmp, ref_tmp, drop=FALSE])
  dists_nonref_tmp = dists_nonref_tmp[, colnames(dists_nonref_tmp) %in% nbs_keep, drop=FALSE]
  
  if(nrow(dists_nonref_tmp) > 1 & ncol(dists_nonref_tmp) > 0){
    if(distsumm == 'avg'){
      dists_summ_tmp = tibble::tibble(barcode=rownames(dists_nonref_tmp),
                                      dist2ref=apply(dists_nonref_tmp, 1, mean))
    } else{
      dists_summ_tmp = tibble::tibble(barcode=rownames(dists_nonref_tmp),
                                      dist2ref=apply(dists_nonref_tmp, 1, min))
    }
  } else{
    dists_summ_tmp = tibble::tibble()
  }
  
  # Apply distance limit if specified
  if(!is.null(limit) & nrow(dists_summ_tmp) > 1){
    dist2reflower = min(dists_summ_tmp$dist2ref, na.rm=TRUE)
    if(dist2reflower > limit){
      dist2refupper = dist2reflower
    } else{
      dist2refupper = limit
    }
    dists_summ_tmp$dist2ref = ifelse(dists_summ_tmp$dist2ref <= dist2refupper, 
                                     dists_summ_tmp$dist2ref, NA)
  }
  
  return(dists_summ_tmp)
}

##
#' @title STgradient_identify_variable_genes: Identify variable genes
#' @description Select top variable genes for analysis
#' @param x STlist object
#' @param sample_name sample name
#' @param dists_summ distance summary data frame
#' @param topgenes number of top genes to select
#' @return vector of gene names
#' @keywords internal
STgradient_identify_variable_genes = function(x, sample_name, dists_summ, topgenes){
  # Get raw counts for non-reference spots
  raw_cts = x@counts[[sample_name]][, dists_summ$barcode[!is.na(dists_summ$dist2ref)], drop=FALSE]
  
  if(ncol(raw_cts) > 1){
    vargenes = calculate_vst(x=raw_cts)
    vargenes = vargenes[order(-vargenes$vst.variance.standardized), ]
    vargenes = vargenes$gene[1:min(topgenes, nrow(vargenes))]
  } else{
    vargenes = character(0)
  }
  
  return(vargenes)
}

##
#' @title STgradient_extract_expression: Extract expression data
#' @description Extract transformed expression for variable genes with distance data
#' @param x STlist object
#' @param sample_name sample name
#' @param vargenes vector of gene names
#' @param dists_summ distance summary data frame
#' @return data frame with expression and distance data
#' @keywords internal
STgradient_extract_expression = function(x, sample_name, vargenes, dists_summ){
  if(length(vargenes) > 0){
    # Extract expression for variable genes
    vargenes_expr = expandSparse(x@tr_counts[[sample_name]])
    vargenes_expr = vargenes_expr[rownames(vargenes_expr) %in% vargenes, , drop=FALSE]
    
    # Transpose and convert to data.frame
    vargenes_expr = t(vargenes_expr)
    vargenes_expr_df = as.data.frame(vargenes_expr)
    vargenes_expr_df$barcode = rownames(vargenes_expr_df)
    
    # Join with distance data (inner join to keep only spots with valid distances)
    joined_df = merge(vargenes_expr_df, dists_summ, by='barcode', all=FALSE)
    
    # Keep only relevant columns
    gene_cols = vargenes
    joined_df = joined_df[, c('barcode', 'dist2ref', gene_cols), drop=FALSE]
    
    # Join with spatial metadata
    meta_df = x@spatial_meta[[sample_name]][, c('libname', 'ypos', 'xpos'), drop=FALSE]
    names(joined_df)[names(joined_df) == 'barcode'] = 'libname'
    vargenes_expr = merge(joined_df, meta_df, by='libname', all.x=TRUE)
    
    # Filter to valid rows
    vargenes_expr = vargenes_expr[!is.na(vargenes_expr$libname), ]
  } else{
    vargenes_expr = data.frame()
  }
  
  return(vargenes_expr)
}

##
#' @title STgradient_detect_outliers: Detect expression outliers
#' @description Identify outlier spots using IQR-based method for each gene
#' @param vargenes_expr data.frame with expression and distance data
#' @param out_rm logical whether to detect outliers
#' @return list with outlier barcodes per gene
#' @keywords internal
STgradient_detect_outliers = function(vargenes_expr, out_rm){
  if(out_rm && nrow(vargenes_expr) > 0){
    outs_dist2ref = list()
    dfdist2ref = vargenes_expr[, !c('ypos', 'xpos', 'dist2ref'), drop=FALSE]
    
    for(gene in colnames(dfdist2ref)){
      quarts = stats::quantile(dfdist2ref[[gene]], probs=c(0.25, 0.75), na.rm=TRUE)
      iqr_dist2ref = stats::IQR(dfdist2ref[[gene]], na.rm=TRUE)
      low_up_limits = c((quarts[1]-1.5*iqr_dist2ref), (quarts[2]+1.5*iqr_dist2ref))
      
      outlier_mask = dfdist2ref[[gene]] < low_up_limits[1] | dfdist2ref[[gene]] > low_up_limits[2]
      outlier_barcodes = vargenes_expr[outlier_mask, 'libname']
      outs_dist2ref[[gene]] = as.character(outlier_barcodes)
    }
  } else{
    outs_dist2ref = list()
  }
  
  return(outs_dist2ref)
}

##
#' @title STgradient_calculate_correlations: Calculate correlations
#' @description For each gene, fit linear models and calculate Spearman correlations
#' @param vargenes_expr data.frame with expression and distance data
#' @param outs_dist2ref list of outlier barcodes per gene
#' @param robust logical whether to use robust regression
#' @param log_dist logical whether to log-transform distances
#' @param sample_name name of sample
#' @return data.frame with correlation results
#' @keywords internal
STgradient_calculate_correlations = function(vargenes_expr, outs_dist2ref, robust, log_dist, sample_name){
  if(nrow(vargenes_expr) == 0){
    return(data.frame())
  }
  
  genes_sample = setdiff(colnames(vargenes_expr), c('libname', 'ypos', 'xpos', 'dist2ref'))
  
  dist_cor = tibble::tibble()
  
  for(gene in genes_sample){
    # Select columns for this gene
    df_gene = vargenes_expr[, c('dist2ref', gene), drop=FALSE]
    # Ensure numeric
    df_gene$dist2ref = as.numeric(df_gene$dist2ref)
    df_gene[[gene]] = as.numeric(df_gene[[gene]])
    
    lm_res = list(estimate=NA, estimate_p=NA)
    cor_res = list(estimate=NA, p.value=NA)
    
    if(!robust && length(outs_dist2ref) > 0 && gene %in% names(outs_dist2ref)){
      # Remove outliers
      if(length(outs_dist2ref[[gene]]) > 0){
        df_gene_outrm = df_gene[!(rownames(df_gene) %in% outs_dist2ref[[gene]]), ]
      } else{
        df_gene_outrm = df_gene
      }
      
      if(nrow(df_gene_outrm) > 1){
        if(log_dist){
          df_gene_outrm$dist2ref = log(df_gene_outrm$dist2ref + 1e-200)
        }
        
        lm_tmp = lm(as.formula(paste0("`", gene, "` ~ dist2ref")), data=df_gene_outrm)
        lm_summ_tmp = summary(lm_tmp)[['coefficients']]
        if(nrow(lm_summ_tmp) > 1){
          lm_res = list(estimate=lm_summ_tmp[2,1], estimate_p=lm_summ_tmp[2,4])
        }
        
        cor_test_result = tryCatch({
          cor.test(df_gene_outrm$dist2ref, df_gene_outrm[[gene]], method='spearman')
        }, warning=function(w){
          if(grepl('standard deviation is zero', w$message)){
            cor.test(df_gene_outrm$dist2ref, df_gene_outrm[[gene]], method='spearman', exact=FALSE)
          } else if(grepl('ties', w$message)){
            # Ties warning is OK, try with exact=FALSE
            cor.test(df_gene_outrm$dist2ref, df_gene_outrm[[gene]], method='spearman', exact=FALSE)
          } else{
            warning(w)
            return(list(estimate = NA_real_, p.value = NA_real_))
          }
        }, error = function(e){
          return(list(estimate = NA_real_, p.value = NA_real_))
        })
        cor_res = cor_test_result
      }
    } else if(robust){
      df_gene_range = df_gene
      if(nrow(df_gene_range) > 1){
        if(log_dist){
          df_gene_range$dist2ref = log(df_gene_range$dist2ref + 1e-200)
        }
        
        lm_tmp = MASS::rlm(as.formula(paste0("`", gene, "` ~ dist2ref")), data=df_gene_range, maxit=100)
        if(lm_tmp$converged && !is.na(lm_tmp$coefficients[2]) && lm_tmp$coefficients[2] != 0){
          lm_res = list(estimate=summary(lm_tmp)[['coefficients']][2,1],
                        estimate_p=NA)  # Skip p-value for robust
          
          cor_test_result = tryCatch({
            cor.test(df_gene_range$dist2ref, df_gene_range[[gene]], method='spearman')
          }, warning=function(w){
            if(grepl('standard deviation is zero', w$message)){
              cor.test(df_gene_range$dist2ref, df_gene_range[[gene]], method='spearman', exact=FALSE)
            } else if(grepl('ties', w$message)){
              cor.test(df_gene_range$dist2ref, df_gene_range[[gene]], method='spearman', exact=FALSE)
            } else{
              warning(w)
              return(list(estimate = NA_real_, p.value = NA_real_))
            }
          }, error = function(e){
            return(list(estimate = NA_real_, p.value = NA_real_))
          })
          cor_res = cor_test_result
        }
      }
    } else{
      df_gene_range = df_gene
      if(nrow(df_gene_range) > 1){
        if(log_dist){
          df_gene_range$dist2ref = log(df_gene_range$dist2ref + 1e-200)
        }
        
        lm_tmp = lm(as.formula(paste0("`", gene, "` ~ dist2ref")), data=df_gene_range)
        lm_summ_tmp = summary(lm_tmp)[['coefficients']]
        if(nrow(lm_summ_tmp) > 1){
          lm_res = list(estimate=lm_summ_tmp[2,1], estimate_p=lm_summ_tmp[2,4])
        }
        
        cor_test_result = tryCatch({
          cor.test(df_gene_range$dist2ref, df_gene_range[[gene]], method='spearman')
        }, warning=function(w){
          if(grepl('standard deviation is zero', w$message)){
            cor.test(df_gene_range$dist2ref, df_gene_range[[gene]], method='spearman', exact=FALSE)
          } else if(grepl('ties', w$message)){
            cor.test(df_gene_range$dist2ref, df_gene_range[[gene]], method='spearman', exact=FALSE)
          } else{
            warning(w)
            return(list(estimate = NA_real_, p.value = NA_real_))
          }
        }, error = function(e){
          return(list(estimate = NA_real_, p.value = NA_real_))
        })
        cor_res = cor_test_result
      }
    }
    
    # Create result row
    tibble_tmp = tibble::tibble(
      sample_name = sample_name,
      gene = gene,
      lm_coef = lm_res$estimate,
      lm_pval = lm_res$estimate_p,
      spearman_r = cor_res$estimate,
      spearman_r_pval = cor_res$p.value,
      pval_comment = NA_character_
    )
    
    if(nrow(tibble_tmp) == 1){
      dist_cor = dplyr::bind_rows(dist_cor, tibble_tmp)
    }
  }
  
  return(dist_cor)
}

##
#' @title STgradient_format_results: Format results
#' @description Format correlation results with proper column names
#' @param dist_cor correlation results data frame
#' @param distsumm distance summary metric
#' @return formatted data frame
#' @keywords internal
STgradient_format_results = function(dist_cor, distsumm){
  if(nrow(dist_cor) > 0 && ncol(dist_cor) >= 7){
    prefix = paste0(distsumm, '_')
    
    # Calculate adjusted p-values
    spearman_r_pval_adj = stats::p.adjust(dist_cor$spearman_r_pval, method='BH')
    
    # Create new columns with prefix
    dist_cor[[paste0(prefix, 'lm_coef')]] = dist_cor$lm_coef
    dist_cor[[paste0(prefix, 'lm_pval')]] = dist_cor$lm_pval
    dist_cor[[paste0(prefix, 'spearman_r')]] = dist_cor$spearman_r
    dist_cor[[paste0(prefix, 'spearman_r_pval')]] = dist_cor$spearman_r_pval
    dist_cor$spearman_r_pval_adj = spearman_r_pval_adj
    
    # Select and reorder columns
    col_order = c('sample_name', 'gene', 
                  paste0(prefix, 'lm_coef'),
                  paste0(prefix, 'lm_pval'),
                  paste0(prefix, 'spearman_r'),
                  paste0(prefix, 'spearman_r_pval'),
                  'spearman_r_pval_adj',
                  'pval_comment')
    dist_cor = dist_cor[, col_order, drop=FALSE]
    
    # Sort by adjusted p-value
    dist_cor = dist_cor[order(dist_cor$spearman_r_pval_adj), ]
  }
  
  return(dist_cor)
}

##
#' @title STgradient_cleanup: Final cleanup
#' @description Remove empty samples and report timing
#' @param results_ls list of results
#' @param zero_t start time
#' @param verbose verbosity flag
#' @return cleaned list of results
#' @keywords internal
STgradient_cleanup = function(results_ls, zero_t, verbose){
  sample_rm = c()
  for(i in names(results_ls)){
    if(nrow(results_ls[[i]]) == 0){
      sample_rm = append(sample_rm, i)
    }
  }
  if(length(sample_rm) > 0){
    results_ls = results_ls[!(names(results_ls) %in% sample_rm)]
  }
  
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('STgradient completed in ', round(end_t, 2), ' min.\n'))
  }
  
  return(results_ls)
}
