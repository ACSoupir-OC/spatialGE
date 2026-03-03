##
# STgradient Core Modular Functions
# Helper functions for spatial gradient analysis

# Load required packages
library(dplyr)
library(tibble)

##
#' @title STgradient_validate_input: Validate input parameters
#' @description Validate sample names, annotation columns, reference clusters, and other parameters
#' @param x an STlist object
#' @param samples vector of sample names to process
#' @param topgenes number of top variable genes
#' @param annot annotation column name
#' @param ref reference cluster
#' @param exclude optional exclude cluster
#' @param min_nb minimum number of neighbors
#' @param nb_dist_thr neighborhood distance threshold
#' @param verbose verbosity level
#' @return list with validated inputs
#' @keywords internal
STgradient_validate_input = function(x, samples, topgenes, annot, ref, exclude, min_nb, nb_dist_thr, cores, verbose){

  # Make sure the reference cluster is character
  ref = as.character(ref)

  # Define samples using names (convert indexes to names if necessary)
  if(is.null(samples)){
    samplenames = names(x@tr_counts)
  } else{
    if(is.numeric(samples)){
      samplenames = as.vector(na.omit(names(x@tr_counts)[samples]))
    } else{
      samplenames = samples[samples %in% names(x@tr_counts)]
    }
    # Verify that sample names exist
    if(length(samplenames) == 0 | !any(samplenames %in% names(x@tr_counts))){
      raise_err(err_code="error0041")
    }
  }

  # Remove samples for which the requested annotation is not present
  sample_rm = c()
  for(i in samplenames){
    if( !(annot %in% colnames(x@spatial_meta[[i]])) ){
      sample_rm = append(sample_rm, i)
      if(verbose){
        message(paste0('The annotation specified in `annot` is not present in ', i, '\n'))
      }
    }
  }
  samplenames = samplenames[ !(samplenames %in% sample_rm) ]

  # Remove samples that do not have the requested reference cluster
  sample_rm = c()
  for(i in samplenames){
    if(sum(x@spatial_meta[[i]][[annot]] == ref) < 1){
      sample_rm = append(sample_rm, i)
    }
  }
  samplenames = samplenames[ !(samplenames %in% sample_rm) ]

  # Define number of cores to use
  if(.Platform$OS.type == 'windows'){
    cores = 1
  } else{
    if(is.null(cores)){
      cores = count_cores(length(samplenames))
    }
  }
  # Ensure cores is defined even if not passed
  if(!exists("cores") || is.null(cores)){
    cores = 1
  }

  # Define neighborhood tolerance
  if(is.null(nb_dist_thr) | !is.numeric(nb_dist_thr) | length(nb_dist_thr) != 2){
    nb_dist_thr = c(0.75, 1.25)
    if(x@misc[['platform']] != 'visium'){
      nb_dist_thr = c(0.25, 3)
    }
  }

  return(list(
    samples = samplenames,
    ref = ref,
    annot = annot,
    exclude = exclude,
    min_nb = min_nb,
    nb_dist_thr = nb_dist_thr,
    cores = cores
  ))
}


##
#' @title STgradient_prepare_distances: Calculate Euclidean distances
#' @description Calculate Euclidean distance matrix and categorize spots by cluster
#' @param x an STlist object
#' @param sample_name name of sample to process
#' @param annot annotation column name
#' @param ref reference cluster
#' @param exclude optional exclude cluster
#' @return list with distance matrix and categorized spots
#' @keywords internal
STgradient_prepare_distances = function(x, sample_name, annot, ref, exclude){
  # Calculate euclidean distances
  # Need to keep annotation column for filtering
  coords_df = x@spatial_meta[[sample_name]][, c('libname', 'ypos', 'xpos', annot), drop=FALSE]
  rownames(coords_df) = coords_df$libname
  coords_mat = coords_df[, c('xpos', 'ypos'), drop=FALSE]
  dist_tmp = as.matrix(stats::dist(coords_mat, method='euclidean'))

  # Save spots in the different categories (ref, nonref, excl)
  ref_tmp = coords_df[['libname']][coords_df[[annot]] == ref]
  nonref_tmp = coords_df[['libname']][!(coords_df[[annot]] %in% c(ref, exclude))]

  return(list(
    dist_tmp = dist_tmp,
    ref_tmp = ref_tmp,
    nonref_tmp = nonref_tmp
  ))
}


##
#' @title STgradient_filter_neighbors: Filter reference spots by neighbor count
#' @description Identify reference spots with sufficient neighbors within neighborhood tolerance
#' @param dist_tmp full distance matrix
#' @param ref_tmp reference spot barcodes
#' @param min_nb minimum number of neighbors required
#' @param nb_dist_thr neighborhood distance threshold (2-element vector)
#' @return vector of barcodes to keep (reference spots with enough neighbors)
#' @keywords internal
STgradient_filter_neighbors = function(dist_tmp, ref_tmp, min_nb, nb_dist_thr){
  # Get minimum distance among all spots within a sample (for Visium would be approximately the same for any sample)
  min_sample = min(as.data.frame(dist_tmp[lower.tri(dist_tmp)]))
  # Get distances among reference spots
  dists_ref_tmp = dist_tmp[ref_tmp, ref_tmp, drop=F]

  # Get number of neighbors within minimum distance
  nbs = colSums(dists_ref_tmp >= min_sample * nb_dist_thr[1] & dists_ref_tmp <= min_sample * nb_dist_thr[2])

  if(sum(nbs >= min_nb) < 1){ # At least 1 cluster of spots to continue with analysis
    nbs_keep = c()
  } else{
    nbs_keep = names(nbs)[nbs >= min_nb] # Save spots to be kept (enough neighbors)
  }

  return(nbs_keep)
}


##
#' @title STgradient_summarize_distances: Summarize distances from reference
#' @description Calculate min or avg distance from non-reference spots to reference spots
#' @param dist_tmp full distance matrix
#' @param nonref_tmp non-reference spot barcodes
#' @param ref_tmp reference spot barcodes
#' @param nbs_keep barcodes of reference spots with enough neighbors
#' @param distsumm distance summary method ('min' or 'avg')
#' @param limit optional distance limit
#' @return data.frame with barcode and distance to reference
#' @keywords internal
STgradient_summarize_distances = function(dist_tmp, nonref_tmp, ref_tmp, nbs_keep, distsumm, limit){
  # Get summarized distances from the reference for each spot in the non-reference
  # Select spots in analysis (non reference in rows, reference in columns)
  dists_nonref_tmp = as.data.frame(dist_tmp[nonref_tmp, ref_tmp, drop=F])
  
  # Remove columns corresponding to spots without enough neighbors
  if(length(nbs_keep) > 0){
    dists_nonref_tmp = dists_nonref_tmp[, nbs_keep, drop=FALSE]
  } else{
    dists_nonref_tmp = data.frame()
  }

  # Check that distances are available for the comparison
  # Number of rows larger than 1, because cannot compute variable genes with a single non-reference spot
  if(nrow(dists_nonref_tmp) > 1 & ncol(dists_nonref_tmp) > 0){
    if(distsumm == 'avg'){ # Summarize non-reference spots using minimum or mean distance
      dists_summ_tmp = data.frame(
        barcode = rownames(dists_nonref_tmp),
        dist2ref = rowMeans(dists_nonref_tmp)
      )
    } else{
      dists_summ_tmp = data.frame(
        barcode = rownames(dists_nonref_tmp),
        dist2ref = do.call(pmin, c(dists_nonref_tmp, na.rm = TRUE))
      )
    }
  } else{
    dists_summ_tmp = data.frame()
  }

  # Remove distances if outside user-specified limit
  if(!is.null(limit) & nrow(dists_summ_tmp) > 1){
    # Get lower and upper distance limits
    # If lower limit is higher than user limit, then set lower limit as upper limit
    if(!all(is.na(dists_summ_tmp$dist2ref))){
      dist2reflower = min(dists_summ_tmp$dist2ref, na.rm=T)
      if(dist2reflower > limit){
        dist2refupper = dist2reflower
      } else{
        dist2refupper = limit
      }
      # Make NA the distances outside range
      dists_summ_tmp$dists2ref[dists_summ_tmp$dist2ref > dist2refupper] = NA
    }
  }

  return(dists_summ_tmp)
}


##
#' @title STgradient_identify_variable_genes: Identify top variable genes
#' @description Identify top genes by variance using VST method
#' @param x an STlist object
#' @param sample_name name of sample to process
#' @param dists_summ data.frame with distance data
#' @param topgenes number of top genes to select
#' @return vector of gene names
#' @keywords internal
STgradient_identify_variable_genes = function(x, sample_name, dists_summ, topgenes){
  # Get expression from variable genes
  # Genes are identified within the range limit
  # Extract expression data (non-transformed counts to be passed to FindVariableFeatures)
  raw_cts = x@counts[[sample_name]]
  # Get spots that have at least 1 distance value
  # However, if only one spot, then sample will be removed from analysis as cannot detect variable genes from single spot
  valid_barcodes = dists_summ[['barcode']][ !is.na(dists_summ[['dist2ref']]) ]
  raw_cts = raw_cts[, valid_barcodes, drop=F]

  # Number of rows larger than 1, because cannot compute variable genes with a single non-reference spot
  # Variable genes in minimum distance range
  if(ncol(raw_cts) > 1){
    vargenes_df = calculate_vst(x=raw_cts)
    # Sort by variance descending (use backticks for column names with dots)
    vargenes_df = vargenes_df[order(-vargenes_df$`vst.variance.standardized`), ]
    vargenes = vargenes_df$gene[1:min(topgenes, nrow(vargenes_df))]
  } else{
    vargenes = c()
  }

  return(vargenes)
}


##
#' @title STgradient_extract_expression: Extract expression and join with distances
#' @description Extract transformed expression data for variable genes and join with distance data
#' @param x an STlist object
#' @param sample_name name of sample to process
#' @param vargenes vector of variable gene names
#' @param dists_summ data.frame with barcode and distance
#' @return data.frame with expression data and distances
#' @keywords internal
STgradient_extract_expression = function(x, sample_name, vargenes, dists_summ){
  if(length(vargenes) > 0){
    # Step 1: Extract expression for variable genes (sparse matrix)
    vargenes_expr = expandSparse(x@tr_counts[[sample_name]])
    vargenes_expr = vargenes_expr[(rownames(vargenes_expr) %in% vargenes), , drop=FALSE]
    
    # Step 2: Transpose and convert to data.frame
    vargenes_expr = t(vargenes_expr)
    vargenes_expr_df = as.data.frame(vargenes_expr)
    vargenes_expr_df$barcode = rownames(vargenes_expr_df)
    
    # Step 3: Prepare distance data
    dists_summ_df = dists_summ
    
    # Step 4: Left join expression with distance data
    joined_df = merge(vargenes_expr_df, dists_summ_df, by='barcode', all.x=TRUE)
    
    # Step 5: Filter to only keep gene columns (and barcode, dist2ref)
    gene_cols = vargenes
    joined_df = joined_df[, c('barcode', 'dist2ref', gene_cols), drop=FALSE]
    
    # Step 6: Prepare spatial metadata
    meta_df = x@spatial_meta[[sample_name]][, c('libname', 'ypos', 'xpos'), drop = FALSE]
    
    # Step 7: Rename barcode in joined data and join with metadata
    names(joined_df)[names(joined_df) == 'barcode'] = 'libname'
    vargenes_expr = merge(joined_df, meta_df, by='libname', all.x=TRUE)
    
    # Remove rows where join failed
    vargenes_expr = vargenes_expr[!is.na(vargenes_expr$libname), ]
    
    # Step 8: Set rownames and keep libname column for later use
    # rownames(vargenes_expr) = vargenes_expr$libname
    # DO NOT remove libname column - calculate_correlations needs it for outlier removal
    
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
  if(out_rm){
    outs_dist2ref = list()
    # Filter to rows with valid dist2ref
    dfdist2ref = vargenes_expr[!is.na(vargenes_expr$dist2ref), ]
    
    # Remove spatial and distance columns for outlier detection
    dfgenes = dfdist2ref[, !c('ypos', 'xpos', 'dist2ref'), drop=FALSE]

    for(gene in colnames(dfgenes)){
      # Calculate gene expression quartiles
      quarts = stats::quantile(dfgenes[[gene]], probs=c(0.25, 0.75))
      # Calculate inter-quartile range
      iqr_dist2ref = stats::IQR(dfgenes[[gene]])
      # Calculate distribution lower and upper limits
      low_up_limits = c((quarts[1]-1.5*iqr_dist2ref),
                        (quarts[2]+1.5*iqr_dist2ref))

      # Save outliers (barcodes)
      outlier_mask = dfgenes[[gene]] < low_up_limits[1] | dfgenes[[gene]] > low_up_limits[2]
      outlier_barcodes = dfdist2ref[outlier_mask, 'libname']
      outs_dist2ref[[gene]] = as.character(outlier_barcodes)
    }
  } else{
    outs_dist2ref = list()
  }

  return(outs_dist2ref)
}


##
#' @title STgradient_calculate_correlations: Calculate linear models and Spearman correlations
#' @description For each gene, fit linear models and calculate Spearman correlations with distance
#' @param vargenes_expr data.frame with expression and distance data
#' @param outs_dist2ref list of outlier barcodes per gene
#' @param robust logical whether to use robust regression
#' @param log_dist logical whether to log-transform distances
#' @param sample_name name of sample
#' @return data.frame with correlation results
#' @keywords internal
STgradient_calculate_correlations = function(vargenes_expr, outs_dist2ref, robust, log_dist, sample_name){
  # Initialize data.frame to store results
  dist_cor = data.frame(
    sample_name = character(),
    gene = character(),
    lm_coef = numeric(),
    lm_pval = numeric(),
    spearman_r = numeric(),
    spearman_r_pval = numeric(),
    pval_comment = character(),
    stringsAsFactors = FALSE
  )

  # CORRELATIONS DISTANCE TO REFERENCE CLUSTER
  genes_sample = setdiff(colnames(vargenes_expr), c('ypos', 'xpos', 'dist2ref'))

  for(gene in genes_sample){
    # Select columns for this gene (keep original column name)
    df_gene = vargenes_expr[, c('libname', 'dist2ref', gene), drop=FALSE]

    lm_res = list(estimate=NA, estimate_p=NA)
    cor_res = list(estimate=NA, p.value=NA)

    if(length(outs_dist2ref) > 0 && gene %in% names(outs_dist2ref) & !robust){ # Regular linear models after removal of outliers
      # Remove outliers (use rownames for barcode matching)
      if(length(outs_dist2ref[[gene]]) > 0){
        df_gene_outrm = df_gene[!(rownames(df_gene) %in% outs_dist2ref[[gene]]), ]
      } else{
        df_gene_outrm = df_gene
      }

      if(nrow(df_gene_outrm) > 1){
        # log-transform distances if selected by user
        if(log_dist){
          df_gene_outrm$dist2ref = log(df_gene_outrm$dist2ref + 1e-200)
        }

        # Run linear model and get summary
        lm_tmp = lm(as.formula(paste(gene, "~ dist2ref")), data=df_gene_outrm)
        lm_summ_tmp = summary(lm_tmp)[['coefficients']]
        if(nrow(lm_summ_tmp) > 1){ # Test a linear model could be run
          lm_res = list(estimate=lm_summ_tmp[2,1],
                        estimate_p=lm_summ_tmp[2,4])
        }
        # Calculate Spearman correlation
        cor_res = tryCatch({cor.test(df_gene_outrm[['dist2ref']], df_gene_outrm[[gene]], method='spearman')}, warning=function(w){return(w)})
        pval_warn = NA_character_
        if(any(class(cor_res) == 'simpleWarning')){ # Let known user if p-value could not be exactly calculated
          if(grepl('standard deviation is zero', cor_res$message)){
            pval_warn = 'zero_st_deviation'
          }
          cor_res = cor.test(df_gene_outrm[['dist2ref']], df_gene_outrm[[gene]], method='spearman', exact=F)
        }

      }

    } else {
      if(robust){ # Robust linear models?
        df_gene_range = df_gene
        if(nrow(df_gene_range) > 1){
          pval_warn = NA_character_

          # log-transform distances if selected by user
          if(log_dist){
            df_gene_range$dist2ref = log(df_gene_range$dist2ref + 1e-200)
          }

          # Run robust linear model and get summary
          lm_tmp = MASS::rlm(as.formula(paste(gene, "~ dist2ref")), data=df_gene_range, maxit=100)
          if(lm_tmp[['converged']] & lm_tmp[['coefficients']][2] != 0){ # Check the model converged and an effect was estimated
            # Run Wald test (MASS::rlm does not provide a p-value)
            lm_test_tmp = sfsmisc::f.robftest(lm_tmp)
            lm_res = list(estimate=summary(lm_tmp)[['coefficients']][2,1],
                          estimate_p=lm_test_tmp[['p.value']])
            # Calculate Spearman correlation
            cor_res = tryCatch({cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman')}, warning=function(w){return(w)})
            if(any(class(cor_res) == 'simpleWarning')){ # Let known user if p-value could not be exactly calculated
              if(grepl('standard deviation is zero', cor_res$message)){
                pval_warn = 'zero_st_deviation'
              }
              cor_res = cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman', exact=F)
            }
          } else{
            pval_warn = 'rob_regr_no_convergence'
          }

        }
      } else{ # Regular linear models without outlier removal
        df_gene_range = df_gene
        if(nrow(df_gene_range) > 1){

          # log-transform distances if selected by user
          if(log_dist){
            df_gene_range$dist2ref = log(df_gene_range$dist2ref + 1e-200)
          }

          lm_tmp = lm(as.formula(paste(gene, "~ dist2ref")), data=df_gene_range)
          lm_summ_tmp = summary(lm_tmp)[['coefficients']]
          if(nrow(lm_summ_tmp) > 1){ # Test a linear model could be run
            lm_res = list(estimate=lm_summ_tmp[2,1],
                          estimate_p=lm_summ_tmp[2,4])
          }
          # Calculate Spearman correlation
          cor_res = tryCatch({cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman')}, warning=function(w){return(w)})
          pval_warn = NA_character_
          if(any(class(cor_res) == 'simpleWarning')){ # Let known user if p-value could not be exactly calculated
            if(grepl('standard deviation is zero', cor_res$message)){
              pval_warn = 'zero_st_deviation'
            }
            cor_res = cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman', exact=F)
          }

        }
      }
    }

    # Create row with results
    tibble_tmp = data.frame(
      sample_name = sample_name,
      gene = gene,
      lm_coef = lm_res[['estimate']],
      lm_pval = lm_res[['estimate_p']],
      spearman_r = as.vector(cor_res[['estimate']]),
      spearman_r_pval = cor_res[['p.value']],
      pval_comment = pval_warn,
      stringsAsFactors = FALSE
    )

    # Add row to result table if there is one row with results
    if(nrow(tibble_tmp) == 1){
      dist_cor = rbind(dist_cor, tibble_tmp)
    }
  }

  return(dist_cor)
}


##
#' @title STgradient_format_results: Format results with p-value adjustment
#' @description Adjust p-values and rename columns with distance summary prefix
#' @param dist_cor data.frame with correlation results
#' @param distsumm distance summary metric ('min' or 'avg')
#' @return formatted data.frame with adjusted p-values
#' @keywords internal
STgradient_format_results = function(dist_cor, distsumm){
  if(nrow(dist_cor) > 0){
    # Adjust p-values for multiple comparison
    dist_cor$spearman_r_pval_adj = p.adjust(dist_cor$spearman_r_pval, method='BH')
    
    # Sort by adjusted p-value
    dist_cor = dist_cor[order(dist_cor$spearman_r_pval_adj), ]
    
    # Add pval_comment column if not present
    if(!'pval_comment' %in% colnames(dist_cor)){
      dist_cor$pval_comment = NA_character_
    }
    
    # Create new column names with distsumm prefix
    old_names = c('sample_name', 'gene', 'lm_coef', 'lm_pval', 'spearman_r', 'spearman_r_pval', 'pval_comment')
    new_names = c('sample_name', 'gene', paste0(distsumm, '_lm_coef'), paste0(distsumm, '_lm_pval'),
                  paste0(distsumm, '_spearman_r'), paste0(distsumm, '_spearman_r_pval'), paste0(distsumm, '_pval_comment'))
    
    # Rename columns
    names(dist_cor)[names(dist_cor) %in% old_names] = new_names
    
    # Rearrange: move spearman_r_pval_adj to end (after spearman_r_pval)
    dist_cor = dist_cor[, c('sample_name', 'gene', 
                            paste0(distsumm, '_lm_coef'), paste0(distsumm, '_lm_pval'),
                            paste0(distsumm, '_spearman_r'), paste0(distsumm, '_spearman_r_pval'),
                            'spearman_r_pval_adj', paste0(distsumm, '_pval_comment')), drop=FALSE]
  }

  return(dist_cor)
}


##
#' @title STgradient_cleanup: Remove empty samples and report timing
#' @description Clean results by removing samples with no output and report completion time
#' @param results_ls list of results data.frames
#' @param zero_t start time
#' @param verbose verbosity level
#' @return cleaned list of results
#' @keywords internal
STgradient_cleanup = function(results_ls, zero_t, verbose){
  names(results_ls) = names(results_ls)

  sample_rm = c()
  for(i in names(results_ls)){
    if(nrow(results_ls[[i]]) == 0){
      sample_rm = append(sample_rm, i)
    }
  }
  if(length(sample_rm) > 0){
    results_ls = results_ls[ !(names(results_ls) %in% sample_rm) ]
  }

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('STgradient completed in ', round(end_t, 2), ' min.\n'))
  }

  return(results_ls)
}