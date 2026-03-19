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

  # Validate inputs and call core
  if (is.null(x) || !methods::is(x, "STlist")) {
    stop('Input x must be an STlist object')
  }
  
  if (is.null(genes)) {
    stop('Please enter one or more genes to calculate statistics.')
  }
  
  # Convert numeric sample indices to names
  if (!is.null(samples) && is.numeric(samples)) {
    samples = names(x@tr_counts)[samples]
  }
  
  # Report genes not present in samples
  if (!is.null(samples)) {
    for (i in samples) {
      notgenes = genes[!(genes %in% rownames(x@tr_counts[[i]]))]
      if (!rlang::is_empty(notgenes)) {
        message(paste0(paste(notgenes, collapse=', '), ": Not present in the transformed counts for sample ", i), ".\n")
      }
    }
  }
  
  # Call core implementation
  result = SThet_core(
    x = x,
    genes = genes,
    samples = samples,
    method = method,
    k = k,
    overwrite = overwrite,
    cores = cores
  )

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('SThet completed in ', round(end_t, 2), ' min.\n'))
  }

  return(result)
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

