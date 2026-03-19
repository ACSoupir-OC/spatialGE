##
# SThet_invdist_test.R - Modular Implementation with Statistical Tests
#
# This file contains the alternative SThet implementation using spdep::moran.test 
# and spdep::geary.test for statistical significance testing.
#
# Created: 2026-03-19
# Pattern: validation → prepare → calculate_with_tests → format
#

# ============================================================================
# SThet_invdist_test - Public interface (modular version)
# ============================================================================

##
#' @title SThet_invdist_test: Computes spatial autocorrelation with statistical tests
#' @description Alternative implementation using spdep::moran.test and spdep::geary.test
#'
#' @param x an STlist
#' @param genes a vector of gene names to compute statistics
#' @param samples the samples to compute statistics
#' @param method The spatial statistic(s) to estimate. It can be set to 'moran',
#' 'geary' or both. Default is 'moran'
#' @param k the number of neighbors to estimate weights
#' @param overwrite logical indicating if previous statistics should be overwritten
#' @param cores the number of cores to use during computations
#' @return an STlist containing spatial statistics
#'
#' @export
#
SThet_invdist_test = function(x=NULL, genes=NULL, samples=NULL, method='moran', k=NULL, overwrite=T, cores=NULL, verbose=TRUE){
  # Record time
  zero_t = Sys.time()
  if(verbose){
    cat(paste0('SThet_invdist_test started.\n'))
  }
  
  # Validate inputs
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
  result = SThet_invdist_test_core(
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
    cat(paste0('SThet_invdist_test completed in ', round(end_t, 2), ' min.\n'))
  }
  
  return(result)
}


# ============================================================================
# SThet_invdist_test_core - Core implementation with statistical tests
# ============================================================================

##
# @title SThet_invdist_test_core
# @description Core implementation with statistical significance testing
#
# @param x an STlist
# @param genes a vector of gene names
# @param samples samples to compute statistics
# @param method spatial statistic(s): 'moran', 'geary', or both
# @param k number of neighbors for k-nearest neighbors
# @param overwrite logical indicating if previous statistics should be overwritten
# @param cores number of cores for parallelization
# @return an STlist with spatial statistics and p-values in gene_meta
#
SThet_invdist_test_core = function(x, genes, samples, method = 'moran', k = NULL,
                                    overwrite = TRUE, cores = NULL) {
  
  # Validate input parameters (reuse SThet validation)
  validated = SThet_validate_input(
    x = x,
    genes = genes,
    samples = samples,
    method = method,
    k = k,
    overwrite = overwrite,
    cores = cores
  )
  
  # Prepare data (reuse SThet preparation)
  prepared = SThet_prepare_data(
    x = x,
    genes = validated$genes,
    samples = validated$samples
  )
  
  # Calculate spatial weights if needed
  if (overwrite || is.null(x@misc[['sthet']][['listws']])) {
    if (!is.null(k)) {
      x@misc[['sthet']][['listws']] = create_listw_from_knn(x, ks = k)
    } else {
      x@misc[['sthet']][['listws']] = create_listw_from_dist(x, cores = cores)
    }
  }
  
  # Perform calculations with statistical tests based on method
  if ('moran' %in% method) {
    x = SThet_calculate_moran_test(
      x = x,
      combo = prepared$combo,
      overwrite = overwrite,
      cores = cores
    )
  }
  
  if ('geary' %in% method) {
    x = SThet_calculate_geary_test(
      x = x,
      combo = prepared$combo,
      overwrite = overwrite,
      cores = cores
    )
  }
  
  # Format and return results
  x = SThet_format_results(x, validated$genes, validated$samples, method)
  
  return(x)
}


# ============================================================================
# Calculation Functions with Statistical Tests - Moran's I
# ============================================================================

##
# @title SThet_calculate_moran_test
# @description Calculate Moran's I with statistical significance testing
#
# @param x an STlist with spatial weights already calculated
# @param combo tibble with sample-gene combinations
# @param overwrite logical indicating if previous statistics should be overwritten
# @param cores number of cores for parallelization
# @return x an STlist with Moran's I values in gene_meta
#
SThet_calculate_moran_test = function(x, combo, overwrite = TRUE, cores = NULL) {
  
  # Define cores available - Windows-compatible parallelization
  if (.Platform$OS.type == 'windows') {
    cores = 1
  }
  if (is.null(cores)) {
    cores = count_cores(length(x@spatial_meta))
  } else {
    cores = as.integer(cores)
    if (is.na(cores)) {
      stop('Could not recognize number of cores requested')
    }
  }
  
  # Compute Moran's I using parallel processing with statistical test
  stat_list = parallel::mclapply(
    seq_along(1:length(unique(as.vector(unlist(combo[[1]]))))),
    function(i_combo) {
      i = unique(as.vector(unlist(combo[[1]])))[i_combo]
      genes_tmp = unique(as.vector(unlist(combo[[2]][combo[[1]] == i])))
      stat_list_tmp = list()
      
      for (j in genes_tmp) {
        gene_expr = x@tr_counts[[i]][j, ]
        
        # Use moran.test for statistical significance
        stat_est = spdep::moran.test(
          x = gene_expr,
          listw = x@misc[['sthet']][['listws']][[i]]
        )
        
        stat_list_tmp[[j]] = stat_est
      }
      
      return(stat_list_tmp)
    },
    mc.cores = cores,
    mc.preschedule = FALSE
  )
  
  names(stat_list) = unique(as.vector(unlist(combo[[1]])))
  
  # Store results in STlist gene_meta
  for (i in names(stat_list)) {
    for (j in names(stat_list[[i]])) {
      combo_name = c(i, j)
      
      if (overwrite || is.na(as.vector(
        x@gene_meta[[combo_name[1]]][
          x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2],
          'moran_i'
        ]
      ))) {
        # Extract estimate from test result
        x@gene_meta[[combo_name[1]]][
          x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2],
          'moran_i'
        ] = as.vector(stat_list[[i]][[j]][['estimate']][1])
      }
    }
  }
  
  return(x)
}


# ============================================================================
# Calculation Functions with Statistical Tests - Geary's C
# ============================================================================

##
# @title SThet_calculate_geary_test
# @description Calculate Geary's C with statistical significance testing
#
# @param x an STlist with spatial weights already calculated
# @param combo tibble with sample-gene combinations
# @param overwrite logical indicating if previous statistics should be overwritten
# @param cores number of cores for parallelization
# @return x an STlist with Geary's C values in gene_meta
#
SThet_calculate_geary_test = function(x, combo, overwrite = TRUE, cores = NULL) {
  
  # Define cores available - Windows-compatible parallelization
  if (.Platform$OS.type == 'windows') {
    cores = 1
  }
  if (is.null(cores)) {
    cores = count_cores(length(x@spatial_meta))
  } else {
    cores = as.integer(cores)
    if (is.na(cores)) {
      stop('Could not recognize number of cores requested')
    }
  }
  
  # Compute Geary's C using parallel processing with statistical test
  stat_list = parallel::mclapply(
    seq_along(1:length(unique(as.vector(unlist(combo[[1]]))))),
    function(i_combo) {
      i = unique(as.vector(unlist(combo[[1]])))[i_combo]
      genes_tmp = unique(as.vector(unlist(combo[[2]][combo[[1]] == i])))
      stat_list_tmp = list()
      
      for (j in genes_tmp) {
        gene_expr = x@tr_counts[[i]][j, ]
        
        # Use geary.test for statistical significance
        stat_est = spdep::geary.test(
          x = gene_expr,
          listw = x@misc[['sthet']][['listws']][[i]]
        )
        
        stat_list_tmp[[j]] = stat_est
      }
      
      return(stat_list_tmp)
    },
    mc.cores = cores,
    mc.preschedule = FALSE
  )
  
  names(stat_list) = unique(as.vector(unlist(combo[[1]])))
  
  # Store results in STlist gene_meta
  for (i in names(stat_list)) {
    for (j in names(stat_list[[i]])) {
      combo_name = c(i, j)
      
      if (overwrite || is.na(as.vector(
        x@gene_meta[[combo_name[1]]][
          x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2],
          'geary_c'
        ]
      ))) {
        # Extract estimate from test result
        x@gene_meta[[combo_name[1]]][
          x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2],
          'geary_c'
        ] = as.vector(stat_list[[i]][[j]][['estimate']][1])
      }
    }
  }
  
  return(x)
}
