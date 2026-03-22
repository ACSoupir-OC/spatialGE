##
# STgradient: Spatial gradient analysis - Public interface
#
# This file contains the public STgradient function
# and imports from STgradient_helpers.R for helper functions
#

#' @title STgradient: Tests of gene expression spatial gradients
#' @description
#' `r lifecycle::badge("stable")`
#' 
#' Calculates Spearman's coefficients to detect genes showing expression spatial gradients
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
#' @examples
#' \dontrun{
#' # ============================================================================
#' # STEP 1: Load TNBC test data and create STlist
#' # ============================================================================
#' # Load TNBC test data from package
#' data_path <- system.file("tests/testthat/data/tnbc_bassiouni", package = "spatialGE")
#' 
#' # Get file paths for all 8 samples
#' samples_dirs <- list.dirs(data_path, recursive = FALSE)
#' count_files <- character(length(samples_dirs))
#' coord_files <- character(length(samples_dirs))
#' 
#' for (i in seq_along(samples_dirs)) {
#'   count_files[i] <- file.path(samples_dirs[i], 
#'                               "GSM6433585_092A_filtered_feature_bc_matrix.h5")
#'   coord_files[i] <- file.path(samples_dirs[i], "spatial", 
#'                               "GSM6433585_092A_tissue_positions_list.csv")
#' }
#' 
#' # Load clinical data
#' clin_file <- file.path(data_path, "bassiouni_clinical.csv")
#' 
#' # Create STlist object
#' st_obj <- STlist(rnacounts = count_files,
#'                  spotcoords = coord_files,
#'                  samples = clin_file)
#'
#' # Transform data (variance stabilizing transformation)
#' st_obj <- transform_data(st_obj)
#' # Expected output: STlist object with 8 samples, transformed counts
#'
#' # ============================================================================
#' # STEP 2: Run spatial clustering to get tissue domain annotations
#' # ============================================================================
#' # Run STclust to identify tissue domains
#' st_obj <- STclust(st_obj, ws = 0.025, ks = c(2, 3), cores = 1, verbose = FALSE)
#' # Expected output: Adds 'stclust_spw0.025_k2' and 'stclust_spw0.025_k3' columns
#'
#' # ============================================================================
#' # STEP 3: Define direction vector for spatial gradient analysis
#' # ============================================================================
#' # Extract spatial coordinates
#' coords <- st_obj@spotcoords[[1]]
#' 
#' # Calculate gradient direction (conceptual vector from one edge to another)
#' dir_x <- diff(range(coords$x))  # Total x range ~5000 microns
#' dir_y <- diff(range(coords$y))  # Total y range ~5000 microns
#' 
#' # Direction vector represents gradient from tissue edge to edge
#' # Example: gradient from lower-left (0,0) to upper-right (max_x, max_y)
#' # This defines the spatial axis along which we test gene expression gradients
#'
#' # ============================================================================
#' # STEP 4: Run STgradient to detect gradient genes
#' # ============================================================================
#' # Analyze sample 1, using cluster 1 as reference domain
#' result <- STgradient(
#'   x = st_obj,
#'   samples = 1,
#'   topgenes = 100,  # Use fewer genes for quick demonstration
#'   annot = "stclust_spw0.025_k2",  # Column name from STclust output
#'   ref = "1",  # Reference cluster (tissue domain)
#'   distsumm = "min",  # Use minimum distance metric
#'   cores = 1,
#'   verbose = TRUE
#' )
#' # Expected output: List with 2 data frames (k=2 and k=3 clustering results)
#'
#' # ============================================================================
#' # STEP 5: Identify and visualize top gradient genes
#' # ============================================================================
#' # View first 10 genes (sorted by Spearman correlation rho)
#' print(head(result[[1]], 10))
#' # Expected output: Data frame with columns: gene, rho, p.value, adjusted.p.value
#'
#' # Identify top positive gradient genes (expressed higher closer to reference)
#' top_positive <- result[[1]][order(-result[[1]]$rho), ][1:5, ]
#' print(top_positive)
#' # Expected output: Top 5 genes with strongest positive correlation to proximity
#'
#' # Identify top negative gradient genes (expressed higher farther from reference)
#' top_negative <- result[[1]][order(result[[1]]$rho), ][1:5, ]
#' print(top_negative)
#' # Expected output: Top 5 genes with strongest negative correlation to proximity
#'
#' # ============================================================================
#' # STEP 6: Visualize gradient gene expression patterns
#' # ============================================================================
#' # Plot expression of top positive gradient gene
#' top_gene <- top_positive$gene[1]
#' gene_expr <- st_obj@counts[[1]][top_gene, ]
#' 
#' # Create spatial plot data
#' plot_data <- data.frame(
#'   x = coords$x,
#'   y = coords$y,
#'   expression = gene_expr,
#'   cluster = st_obj@spatial_meta[[1]]$stclust_spw0.025_k2
#' )
#' 
#' # Create spatial scatter plot with ggplot2
#' library(ggplot2)
#' ggplot(plot_data, aes(x = x, y = y, color = expression)) +
#'   geom_point(size = 0.5) +
#'   scale_color_viridis_c() +
#'   labs(title = paste("Expression of", top_gene),
#'        color = "Expression") +
#'   theme_minimal()
#' # Expected output: Spatial heatmap showing gene expression gradient
#'
#' # ============================================================================
#' # STEP 7: Compare different distance metrics
#' # ============================================================================
#' # Run with average distance metric
#' result_avg <- STgradient(
#'   x = st_obj,
#'   samples = 1,
#'   topgenes = 100,
#'   annot = "stclust_spw0.025_k2",
#'   ref = "1",
#'   distsumm = "avg",
#'   cores = 1,
#'   verbose = FALSE
#' )
#' 
#' # Compare correlations between min and avg distance metrics
#' cor(result[[1]]$rho, result_avg[[1]]$rho)
#' # Expected output: Correlation coefficient (~0.95-1.0 for highly correlated genes)
#' }
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
