##
#' @title STgradient: Tests of gene expression spatial gradients
#' @description Calculates Spearman's coefficients to detect genes showing expression spatial gradients
#' @details
#' The `STgradient` function fits linear models and calculates Spearman coefficients
#' between the expression of a gene and the minimum or average distance of spots or
#' cells to a reference tissue domain. In other wordsm the `STgradient` function
#' can be used to investigate if a gene is expressed higher in spots/cells closer to
#' a specific reference tissue domain, compared to spots/cells farther from the
#' reference domain (or viceversa as indicated by the Spearman's cofficient).
#'
#' @param x an STlist with transformed gene expression
#' @param samples the samples on which the test should be executed
#' @param topgenes the number of high-variance genes to be tested. These genes are
#' selected in descending order of variance as caclulated using Seurat's vst method
#' @param annot the name of a column in `@spatial_meta` containing the tissue domain
#' assignments for each spot or cell. These assignments can be generated using the
#' `STclust` function
#' @param ref one of the tissue domains in the column specified in `annot`,
#' corresponding to the "reference" cluster or domain. Spearman's correlations will
#' be calculated using spots assigned to domains other than this reference domain
#' (or domains specified in `exclude`).
#' @param exclude optional, a cluster/domain to exclude from the analysis
#' @param out_rm logical (optional), remove gene expression outliers defined by
#' the interquartile method. This option is only valid when `robust=F`
#' @param limit limite the analysis to spots/cells with distances to `ref` shorther
#' than the value specified here. Useful when gradients might occur at smaller scales
#' or when the domain in `ref` is scattered through the tissue. Caution must be used
#' due to difficult interpretation of imposed limits. It is suggested to run analysis
#' without restricted distances in addition for comparison.
#' @param distsumm the distance summary metric to use in correlations. One of `min` or `avg`
#' @param min_nb the minimum number of immediate neighbors a spot or cell has to
#' have in order to be included in the analysis. This parameter seeks to reduce the
#' effect of isolated `ref` spots on the correlation
#' @param robust logical, whether to use robust regression (`MASS` and `sfsmisc` packages)
#' @param nb_dist_thr a numeric vector of length two indicating the tolerance interval to assign
#' spots/cells to neighborhoods. The wider the range of the interval, the more likely
#' distinct neighbors to be considered. If NULL, `c(0.75, 1.25)` and `c(0.25, 3)` is assigned
#' for Visium and CosMx respectively.
#' @param log_dist logical, whether to apply the natural logarithm to the spot/cell
#' distances. It applies to all distances a constant (1e-200) to avoid log(0)
#' @param cores the number of cores used during parallelization. If NULL (default),
#' the number of cores is defined automatically
#' @param verbose logical, whether to print text to console
#' @return a list of data frames with the results of the test
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom stats IQR lm cor.test p.adjust quantile
#'
STgradient = function(x=NULL, samples=NULL, topgenes=2000, annot=NULL, ref=NULL, exclude=NULL,
                      out_rm=FALSE, limit=NULL, distsumm='min', min_nb=3, robust=TRUE,
                      nb_dist_thr=NULL, log_dist=FALSE, cores=NULL, verbose=TRUE){

  # To prevent NOTES in R CMD check
  . = NULL

  # Record time
  zero_t = Sys.time()

  # Validate inputs
  validation = STgradient_validate_input(x, samples, topgenes, annot, ref, exclude, min_nb, nb_dist_thr, cores, verbose)

  # Parallel processing across samples
  results_ls = parallel::mclapply(validation$samples, function(sample_name){
    # Prepare distances and categorize spots
    dist_data = STgradient_prepare_distances(x, sample_name, validation$annot, validation$ref, validation$exclude)

    # Filter reference spots by neighbor count
    nbs_keep = STgradient_filter_neighbors(dist_data$dist_tmp, dist_data$ref_tmp, validation$min_nb, validation$nb_dist_thr)

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
  }, mc.cores=validation$cores)

  names(results_ls) = validation$samples

  # Clean up and report timing
  results_ls = STgradient_cleanup(results_ls, zero_t, verbose)

  return(results_ls)
}