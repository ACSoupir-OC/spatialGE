##
# SThet_core.R - Core Implementation for Spatial Heterogeneity Analysis
#
# This file contains the core logic for SThet (spatial heterogeneity analysis)
# Computes Moran's I and Geary's C spatial autocorrelation statistics
#
# Created: 2026-03-19
# Pattern: validation → prepare → calculate → format
#

# ============================================================================
# SThet_core - Main core function
# ============================================================================

##
# @title SThet_core
# @description Core implementation of spatial heterogeneity analysis
# @details Computes Moran's I and/or Geary's C for specified genes and samples
#          This is the internal core function - use SThet() for public interface
#
# @param x an STlist
# @param genes a vector of gene names to compute statistics
# @param samples a vector of sample names to compute statistics
# @param method The spatial statistic(s) to estimate: 'moran', 'geary', or both
# @param k the number of neighbors for k-nearest neighbors (NULL = distance-based)
# @param overwrite logical indicating if previous statistics should be overwritten
# @param cores integer indicating the number of cores to use
# @return an STlist with spatial statistics stored in gene_meta
#
SThet_core = function(x, genes, samples, method = 'moran', k = NULL, 
                      overwrite = TRUE, cores = NULL) {
  
  # Validate input parameters
  validated = SThet_validate_input(
    x = x,
    genes = genes,
    samples = samples,
    method = method,
    k = k,
    overwrite = overwrite,
    cores = cores
  )
  
  # Prepare data (extract samples, coordinates, match genes)
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
  
  # Perform calculations based on method
  if ('moran' %in% method) {
    x = SThet_calculate_moran(
      x = x,
      combo = prepared$combo,
      overwrite = overwrite,
      cores = cores
    )
  }
  
  if ('geary' %in% method) {
    x = SThet_calculate_geary(
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
# Validation Functions
# ============================================================================

##
# @title SThet_validate_input
# @description Validate all input parameters for SThet
#
# @param x an STlist
# @param genes a vector of gene names
# @param samples samples to compute statistics (can be numeric indices)
# @param method spatial statistic(s): 'moran', 'geary', or both
# @param k number of neighbors for k-nearest neighbors
# @param overwrite logical indicating if previous statistics should be overwritten
# @param cores number of cores for parallelization
# @return list with validated parameters
#
SThet_validate_input = function(x, genes, samples, method, k, overwrite, cores) {
  
  # Validate STlist
  if (is.null(x) || !methods::is(x, "STlist")) {
    stop('Input x must be an STlist object')
  }
  
  # Validate genes
  if (is.null(genes)) {
    stop('Please enter one or more genes to calculate statistics.')
  }
  
  if (!is.character(genes) || length(genes) == 0) {
    stop('genes must be a non-empty character vector')
  }
  
  # Validate samples - convert numeric indices to names
  if (is.null(samples)) {
    samples = names(x@tr_counts)
  } else if (is.numeric(samples)) {
    samples = names(x@tr_counts)[samples]
  }
  
  if (!is.character(samples) || length(samples) == 0) {
    stop('samples must be NULL, numeric indices, or character vector')
  }
  
  # Validate method
  if (!all(method %in% c('moran', 'geary'))) {
    stop("method must be 'moran', 'geary', or c('moran', 'geary')")
  }
  
  # Validate k (if provided)
  if (!is.null(k)) {
    k = as.integer(k)
    if (is.na(k) || k <= 0) {
      stop("If using k nearest-neighbors, please input a positive integer for k.")
    }
  }
  
  # Validate overwrite
  if (!is.logical(overwrite) || length(overwrite) != 1) {
    stop('overwrite must be a single logical value')
  }
  
  # Validate cores
  if (!is.null(cores)) {
    cores = as.integer(cores)
    if (is.na(cores) || cores < 1) {
      stop('cores must be NULL or a positive integer')
    }
  }
  
  return(list(
    genes = genes,
    samples = samples,
    method = method,
    k = k,
    overwrite = overwrite,
    cores = cores
  ))
}


# ============================================================================
# Data Preparation Functions
# ============================================================================

##
# @title SThet_prepare_data
# @description Extract tissue spots, coordinates, and match genes to samples
#
# @param x an STlist
# @param genes validated gene names
# @param samples validated sample names
# @return list with combo tibble and gene availability info
#
SThet_prepare_data = function(x, genes, samples) {
  
  # Generate combination of sample x gene
  combo = tibble::tibble()
  gene_availability = list()
  
  for (i in samples) {
    # Check if gene names are in the data set
    subsetgenes = genes[genes %in% rownames(x@tr_counts[[i]])]
    combo = dplyr::bind_rows(combo, expand.grid(i, subsetgenes))
    
    # Track genes not present
    notgenes = genes[!(genes %in% rownames(x@tr_counts[[i]]))]
    gene_availability[[i]] = list(
      present = subsetgenes,
      not_present = notgenes
    )
    
    # Add columns in gene meta data if not already present
    if (!('moran_i' %in% colnames(x@gene_meta[[i]]))) {
      x@gene_meta[[i]][['moran_i']] = NA
    }
    if (!('geary_c' %in% colnames(x@gene_meta[[i]]))) {
      x@gene_meta[[i]][['geary_c']] = NA
    }
  }
  
  return(list(
    combo = combo,
    gene_availability = gene_availability
  ))
}


# ============================================================================
# Calculation Functions - Moran's I
# ============================================================================

##
# @title SThet_calculate_moran
# @description Calculate Moran's I spatial autocorrelation for all gene-sample combinations
#
# @param x an STlist with spatial weights already calculated
# @param combo tibble with sample-gene combinations
# @param overwrite logical indicating if previous statistics should be overwritten
# @param cores number of cores for parallelization
# @return x an STlist with Moran's I values in gene_meta
#
SThet_calculate_moran = function(x, combo, overwrite = TRUE, cores = NULL) {
  
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
  
  # Compute Moran's I using parallel processing
  stat_list = parallel::mclapply(
    seq_along(1:length(unique(as.vector(unlist(combo[[1]]))))),
    function(i_combo) {
      i = unique(as.vector(unlist(combo[[1]])))[i_combo]
      genes_tmp = unique(as.vector(unlist(combo[[2]][combo[[1]] == i])))
      stat_list_tmp = list()
      
      for (j in genes_tmp) {
        gene_expr = x@tr_counts[[i]][j, ]
        
        stat_est = spdep::moran(
          x = gene_expr,
          listw = x@misc[['sthet']][['listws']][[i]],
          n = length(x@misc[['sthet']][['listws']][[i]]$neighbours),
          S0 = spdep::Szero(x@misc[['sthet']][['listws']][[i]])
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
        x@gene_meta[[combo_name[1]]][
          x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2],
          'moran_i'
        ] = as.vector(stat_list[[i]][[j]][['I']])
      }
    }
  }
  
  return(x)
}


# ============================================================================
# Calculation Functions - Geary's C
# ============================================================================

##
# @title SThet_calculate_geary
# @description Calculate Geary's C spatial autocorrelation for all gene-sample combinations
#
# @param x an STlist with spatial weights already calculated
# @param combo tibble with sample-gene combinations
# @param overwrite logical indicating if previous statistics should be overwritten
# @param cores number of cores for parallelization
# @return x an STlist with Geary's C values in gene_meta
#
SThet_calculate_geary = function(x, combo, overwrite = TRUE, cores = NULL) {
  
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
  
  # Compute Geary's C using parallel processing
  stat_list = parallel::mclapply(
    seq_along(1:length(unique(as.vector(unlist(combo[[1]]))))),
    function(i_combo) {
      i = unique(as.vector(unlist(combo[[1]])))[i_combo]
      genes_tmp = unique(as.vector(unlist(combo[[2]][combo[[1]] == i])))
      stat_list_tmp = list()
      
      for (j in genes_tmp) {
        gene_expr = x@tr_counts[[i]][j, ]
        
        stat_est = spdep::geary(
          x = gene_expr,
          listw = x@misc[['sthet']][['listws']][[i]],
          n = length(x@misc[['sthet']][['listws']][[i]]$neighbours),
          n1 = length(x@misc[['sthet']][['listws']][[i]]$neighbours) - 1,
          S0 = spdep::Szero(x@misc[['sthet']][['listws']][[i]])
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
        x@gene_meta[[combo_name[1]]][
          x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2],
          'geary_c'
        ] = as.vector(stat_list[[i]][[j]][['C']])
      }
    }
  }
  
  return(x)
}


# ============================================================================
# Results Formatting
# ============================================================================

##
# @title SThet_format_results
# @description Format results and ensure gene_meta is properly structured
#
# @param x an STlist with calculated statistics
# @param genes original gene list
# @param samples original sample list
# @param method methods used
# @return x an STlist with formatted results
#
SThet_format_results = function(x, genes, samples, method) {
  # Results are already stored in gene_meta by calculate functions
  # This function ensures proper structure and can add metadata if needed
  
  # Add metadata about the analysis
  x@misc[['sthet']][['last_run']] = list(
    timestamp = Sys.time(),
    genes = genes,
    samples = samples,
    method = method
  )
  
  return(x)
}
