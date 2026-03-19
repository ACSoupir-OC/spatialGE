##
# STgradient: Spatial gradient analysis - Core implementation
#
# This file contains the core STgradient algorithm
#
##

#' @title STgradient_core: Core algorithm for spatial gradient analysis
#' @description Core implementation of STgradient - calculates Spearman correlations between gene expression and distances to reference domain
#' @details Internal function - use STgradient() for public interface
#' @param x STlist object with transformed gene expression
#' @param samples Vector of sample names to process
#' @param annot Name of annotation column in @spatial_meta
#' @param ref Reference tissue domain
#' @param exclude Optional domain to exclude
#' @param out_rm Remove gene expression outliers (IQR method)
#' @param limit Limit analysis to spots within this distance threshold
#' @param distsumm Distance summary metric: "min" or "avg"
#' @param min_nb Minimum number of neighbors required
#' @param robust Use robust regression
#' @param nb_dist_thr Neighborhood distance threshold
#' @param log_dist Apply log transform to distances
#' @param topgenes Number of high-variance genes to test
#' @param cores Number of cores for parallelization
#' @param verbose Print progress messages
#' @return List of data frames with correlation results
#'
#' @importFrom magrittr %>%
#' @importFrom stats IQR cor.test p.adjust quantile
#'
STgradient_core = function(x, samples, annot, ref, exclude, out_rm, limit, distsumm, min_nb, robust, nb_dist_thr, log_dist, topgenes, cores, verbose){

  # Record time
  zero_t = Sys.time()

  # Parallel processing across samples
  results_ls = parallel::mclapply(samples, function(sample_name){
    # Prepare distances and categorize spots
    dist_data = STgradient_prepare_distances(x, sample_name, annot, ref, exclude)

    # Filter reference spots by neighbor count
    nbs_keep = STgradient_filter_neighbors(dist_data$dist_tmp, dist_data$ref_tmp, min_nb, nb_dist_thr)

    # Summarize distances from reference
    dists_summ = STgradient_summarize_distances(dist_data$dist_tmp, dist_data$nonref_tmp, dist_data$ref_tmp,
                                                 nbs_keep, distsumm, limit)

    if(nrow(dists_summ) > 1){
      # Identify variable genes
      vargenes = STgradient_identify_variable_genes(x, sample_name, dists_summ, topgenes)

      # Extract expression data
      vargenes_expr = STgradient_extract_expression(x, sample_name, vargenes, dists_summ)

      if(nrow(vargenes_expr) > 0){
        # Detect outliers if requested
        if(out_rm & !robust){
          outliers = STgradient_detect_outliers(vargenes_expr, out_rm)
        } else {
          outliers = list()
        }

        # Calculate correlations
        results = STgradient_calculate_correlations(vargenes_expr, outliers, robust, log_dist, sample_name)

        # Format results
        results = STgradient_format_results(results, distsumm)
      } else{
        results = data.frame()
      }

    } else{
      results = data.frame()
    }

    return(results)
  }, mc.cores=cores)

  names(results_ls) = samples

  # Clean up and report timing
  results_ls = STgradient_cleanup(results_ls, zero_t, verbose)

  return(results_ls)
}
