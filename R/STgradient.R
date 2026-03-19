##
# STgradient: Spatial gradient analysis - Public interface
#
# This file contains the public STgradient function
# and imports from STgradient_helpers.R for helper functions
#

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
#' @seealso STgradient_legacy for reproducibility with original implementation
#'
#' @importFrom magrittr %>%
#' @importFrom stats p.adjust
#' @importFrom tibble add_column column_to_rownames
#' @importFrom dplyr mutate select
#' @importFrom parallel mclapply count
STgradient = function(x=NULL, samples=NULL, topgenes=2000, annot=NULL, ref=NULL, exclude=NULL,
                      out_rm=FALSE, limit=NULL, distsumm='min', min_nb=3, robust=TRUE,
                      nb_dist_thr=NULL, log_dist=FALSE, cores=NULL, verbose=TRUE){

  # Validate inputs
  validation = STgradient_validate_input(x, samples, topgenes, annot, ref, exclude, min_nb, nb_dist_thr, cores, verbose)

  # Call core function
  results_ls = STgradient_core(
    x = x,
    samples = validation$samples,
    annot = validation$annot,
    ref = validation$ref,
    exclude = validation$exclude,
    out_rm = out_rm,
    limit = limit,
    distsumm = distsumm,
    min_nb = min_nb,
    robust = robust,
    nb_dist_thr = validation$nb_dist_thr,
    log_dist = log_dist,
    topgenes = topgenes,
    cores = validation$cores,
    verbose = verbose
  )

  return(results_ls)
}
