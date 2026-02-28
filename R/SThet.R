##
#' @title SThet: Computes global spatial autocorrelation statistics on gene expression
#' @description Computes the global spatial autocorrelation statistics Moran's I and/or
#' Geary's C for a set of genes
#' @details The function computes global spatial autocorrelation statistics (Moran's I and/or
#' Geary's C) for the requested genes and samples. Then computation uses the
#' package `spdep`. The calculated statistics are stored in the STlist, which can
#' be accessed with the `get_gene_meta` function. For visual comparative analysis,
#' the function `compare_SThet` can be used afterwards.
#'
#' @param x an STlist
#' @param genes a vector of gene names to compute statistics
#' @param samples the samples to compute statistics
#' @param method The spatial statistic(s) to estimate. It can be set to 'moran',
#' 'geary' or both. Default is 'moran'
#' @param k the number of neighbors to estimate weights. By default NULL, meaning that
#' spatial weights will be estimated from Euclidean distances. If an positive integer is
#' entered, then the faster k nearest-neighbors approach is used. Please keep in mind
#' that estimates are not as accurate as when using the default distance-based method.
#' @param overwrite logical indicating if previous statistics should be overwritten.
#' Default to FALSE (do not overwrite)
#' @param cores integer indicating the number of cores to use during parallelization.
#' If NULL, the function uses half of the available cores at a maximum. The parallelization
#' uses `parallel::mclapply` and works only in Unix systems
#' @param verbose logical, whether to print text to console
#' @return an STlist containing spatial statistics
#'
#' @examples
#' \donttest{
#' # Using included melanoma example (Thrane et al.)
#' # Download example data set from spatialGE_Data
#' thrane_tmp = tempdir()
#' unlink(thrane_tmp, recursive=TRUE)
#' dir.create(thrane_tmp)
#' lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
#' tryCatch({ # In case data is not available from network
#'   download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
#'   #' zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
#'   unzip(zipfile=zip_tmp, exdir=thrane_tmp)
#'   # Generate the file paths to be passed to the STlist function
#'   count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='counts')
#'   coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='mapping')
#'   clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                           full.names=TRUE, pattern='clinical')
#'   # Create STlist
#'   library('spatialGE')
#'   melanoma <- STlist(rnacounts=count_files,
#'                      spotcoords=coord_files,
#'                      samples=clin_file)
#'   melanoma <- transform_data(melanoma)
#'   melanoma <- SThet(melanoma, genes=c('MLANA', 'TP53'), method='moran')
#'   get_gene_meta(melanoma, sthet_only=TRUE)
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#'   return(NULL)
#' })
#' }
#'
#' @export
#
SThet = function(x=NULL, genes=NULL, samples=NULL, method='moran', k=NULL, overwrite=TRUE,
                 cores=NULL, verbose=TRUE){
  # Record time
  zero_t = Sys.time()
  if(verbose){
    cat(paste0('SThet started.\n'))
  }

  # Select sample names if NULL or if number entered
  if (is.null(samples)){
    samples = names(x@tr_counts)
  } else{
    if(is.numeric(samples)){
      samples = names(x@tr_counts)[samples]
    }
  }

  # Check that genes have been input
  if(is.null(genes)){
    stop('Please enter one or more genes to calculate statistics.')
  }

  # Generate combination of sample x gene to for.
  combo = tibble::tibble()
  for(i in samples){
    # Check if gene names are in the data set
    subsetgenes = genes[genes %in% rownames(x@tr_counts[[i]])]
    combo = dplyr::bind_rows(combo, expand.grid(i, subsetgenes))

    # Get genes not present.
    notgenes = genes[!(genes %in% rownames(x@tr_counts[[i]]))]

    if(!rlang::is_empty(notgenes)){
      message(paste0(paste(notgenes, collapse=', '), ": Not present in the transformed counts for sample ", i), ".\n")
    }

    rm(subsetgenes, notgenes) # Clean env

    # Add columns in gene meta data if not already present
    if(!('moran_i' %in% colnames(x@gene_meta[[i]]))){
      x@gene_meta[[i]][['moran_i']] = NA
    }
    if(!('geary_c' %in% colnames(x@gene_meta[[i]]))){
      x@gene_meta[[i]][['geary_c']] = NA
    }
  }

  # Check whether or not a list of weights have been created
  if(overwrite | is.null(x@misc[['sthet']][['listws']])){
    if(verbose){
      cat(paste("\tCalculating spatial weights...\n")) ## Mostly added to make sure calculation is happening only when needed.
    }
    if(!is.null(k)){
      k = as.integer(k)
      if(!is.na(k) & k > 0){
        x@misc[['sthet']][['listws']] = create_listw_from_knn(x, ks=k)
      } else{
        stop("If using k nearest-neighbors, please input a positive integer for k.")
      }
    } else{
      x@misc[['sthet']][['listws']] = create_listw_from_dist(x, cores=cores)
    }
  }

  # Perform calculations
  if('moran' %in% method){
    x = gene_moran_i_notest(x=x, combo=combo, overwrite=overwrite, cores=cores)
  }
  if('geary' %in% method){
    x = gene_geary_c_notest(x=x, combo=combo, overwrite=overwrite, cores=cores)
  }

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('SThet completed in ', round(end_t, 2), ' min.\n'))
  }

  return(x)
}


##
# @title calc_spatial_autocorr
# @description Calculate Moran's I or Geary's C spatial autocorrelation statistics
# @details Shared helper that extracts common logic from gene_moran_i_notest and gene_geary_c_notest
#          to reduce code duplication. Dispatches to spdep::moran() or spdep::geary() based on stat_type,
#          extracts the appropriate result field, and stores it in the correct column.
#
# @param x an STlist with normalized gene counts
# @param combo a tibble with combinations of samples and genes to calculate statistics
# @param stat_type character: 'moran' or 'geary'
# @param overwrite logical indicating if previous statistics should be overwritten
# @param cores integer indicating the number of cores to use during parallelization
# @return x an STlist with the calculated spatial autocorrelation statistics stored in gene_meta
#
calc_spatial_autocorr = function(x=NULL, combo=NULL, stat_type='moran', overwrite=T, cores=NULL){
  # Validate input
  if(!stat_type %in% c('moran', 'geary')){
    stop('stat_type must be "moran" or "geary"')
  }
  
  # Define cores available - Windows-compatible parallelization
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(x@spatial_meta))
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      stop('Could not recognize number of cores requested')
    }
  }
  
  # Determine result column name based on stat_type
  result_col = if(stat_type == 'moran') 'moran_i' else 'geary_c'
  
  # Compute autocorrelation using parallel processing
  stat_list = parallel::mclapply(
    seq_along(1:length(unique(as.vector(unlist(combo[[1]]))))),
    function(i_combo){
      i = unique(as.vector(unlist(combo[[1]])))[i_combo]
      genes_tmp = unique(as.vector(unlist(combo[[2]][combo[[1]] == i])))
      stat_list_tmp = list()
      
      for(j in genes_tmp){
        gene_expr = x@tr_counts[[i]][j, ]
        
        if(stat_type == 'moran'){
          stat_est = spdep::moran(
            x = gene_expr,
            listw = x@misc[['sthet']][['listws']][[i]],
            n = length(x@misc[['sthet']][['listws']][[i]]$neighbours),
            S0 = spdep::Szero(x@misc[['sthet']][['listws']][[i]])
          )
        } else {
          stat_est = spdep::geary(
            x = gene_expr,
            listw = x@misc[['sthet']][['listws']][[i]],
            n = length(x@misc[['sthet']][['listws']][[i]]$neighbours),
            n1 = length(x@misc[['sthet']][['listws']][[i]]$neighbours) - 1,
            S0 = spdep::Szero(x@misc[['sthet']][['listws']][[i]])
          )
        }
        
        stat_list_tmp[[j]] = stat_est
      }
      
      return(stat_list_tmp)
    },
    mc.cores = cores,
    mc.preschedule = FALSE
  )
  
  names(stat_list) = unique(as.vector(unlist(combo[[1]])))
  
  # Store results in STlist gene_meta
  for(i in names(stat_list)){
    for(j in names(stat_list[[i]])){
      combo_name = c(i, j)
      
      if(overwrite | is.na(as.vector(
        x@gene_meta[[combo_name[1]]][
          x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2],
          result_col
        ]
      ))){
        result_value = if(stat_type == 'moran') {
          as.vector(stat_list[[i]][[j]][['I']])
        } else {
          as.vector(stat_list[[i]][[j]][['C']])
        }
        
        x@gene_meta[[combo_name[1]]][
          x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2],
          result_col
        ] = result_value
      }
    }
  }
  
  return(x)
}


##
# @title gene_moran_i_notest
# @description Calculates Moran's I from ST data (wrapper for calc_spatial_autocorr)
# @param x an STlist with normalized gene counts.
# @param combo a table with combinations of samples and genes to calculate statistics
# @return x a STlist including the calculated Moran's I
# @note Wrapper for backward compatibility
#
gene_moran_i_notest = function(x=NULL, combo=NULL, overwrite=T, cores=NULL){
  calc_spatial_autocorr(x = x, combo = combo, stat_type = 'moran', overwrite = overwrite, cores = cores)
}


##
# @title gene_geary_c_notest
# @description Calculates Geary's C from ST data (wrapper for calc_spatial_autocorr)
# @param x an STlist with normalized gene counts.
# @param combo a table with combinations of samples and genes to calculate statistics
# @return x a STlist including the calculated Geary's C
# @note Wrapper for backward compatibility
#
gene_geary_c_notest = function(x=NULL, combo=NULL, overwrite=T, cores=NULL){
  calc_spatial_autocorr(x = x, combo = combo, stat_type = 'geary', overwrite = overwrite, cores = cores)
}

