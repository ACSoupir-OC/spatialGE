##
# STdiff: Differential expression analysis for spatial transcriptomics data
# 
# Main STdiff() function as dispatcher calling modular functions for
# non-spatial tests, spatial models, and result compilation.
#
# @importFrom magrittr %>%
# @importFrom rlang :=
# @importFrom stats p.adjust

# Source modular functions
helper_path <- system.file('R', 'STdiff_helpers.R', package='spatialGE', mustWork=FALSE)
if(helper_path == '' || !file.exists(helper_path)){ source('R/STdiff_helpers.R') } else { source(helper_path) }

nonspatial_path <- system.file('R', 'STdiff_nonspatial.R', package='spatialGE', mustWork=FALSE)
if(nonspatial_path == '' || !file.exists(nonspatial_path)){ source('R/STdiff_nonspatial.R') } else { source(nonspatial_path) }

spatial_path <- system.file('R', 'STdiff_spatial.R', package='spatialGE', mustWork=FALSE)
if(spatial_path == '' || !file.exists(spatial_path)){ source('R/STdiff_spatial.R') } else { source(spatial_path) }

results_path <- system.file('R', 'STdiff_results.R', package='spatialGE', mustWork=FALSE)
if(results_path == '' || !file.exists(results_path)){ source('R/STdiff_results.R') } else { source(results_path) }

# Legacy function for reproducibility
legacy_path <- system.file('R', 'STdiff_legacy.R', package='spatialGE', mustWork=FALSE)
if(legacy_path == '' || !file.exists(legacy_path)){
  if(file.exists('R/STdiff_legacy.R')){ source('R/STdiff_legacy.R') }
} else { source(legacy_path) }

##
#' @title STdiff: Differential gene expression analysis for spatial transcriptomics
#' @description
#' `r lifecycle::badge("stable")`
#' 
#' Tests for differentially expressed genes between groups of spots/cells
#' in spatial transcriptomics data. First runs non-spatial tests (linear models, t-tests,
#' or Wilcoxon tests) to detect DE genes, then optionally fits spatial mixed models
#' with Matern covariance structure on a subset of DE genes.
#'
#' @details
#' Tests for genes with significantly higher or lower expression in one group of
#' spots/cells with respect to the rest. Users can specify cluster assignments via
#' metadata columns or via STclust parameters (w, k, deepSplit). Uses spaMM::fitme
#' and is computationally expensive. Set sp_topgenes=0 to run non-spatial tests only.
#'
#' @param x an STlist object containing spatial transcriptomics data
#' @param samples vector of sample names or indices to test. If NULL, uses all samples
#' @param annot column name in spatial_meta for cluster annotations. Required if w/k not specified
#' @param w spatial weight parameter for STclust (required if annot is NULL)
#' @param k number of clusters for STclust, or 'dtc' for dynamicTreeCut (required if annot is NULL)
#' @param deepSplit deepSplit parameter for dynamicTreeCut clusters (required if k='dtc')
#' @param topgenes number of top variable genes to select (default=5000). If NULL, all genes used
#' @param pval_thr p-value threshold for selecting DE genes from non-spatial tests (default=0.05)
#' @param pval_adj p-value adjustment method (default='fdr', passed to stats::p.adjust)
#' @param test_type type of test: 'mm' (linear models), 't_test', or 'wilcoxon'
#' @param sp_topgenes proportion (0-1) of DE genes to fit spatial models on. If 0, no spatial tests
#' @param clusters optional vector of specific cluster names to test (vs all clusters)
#' @param pairwise whether to perform pairwise tests (TRUE) or reference-based (FALSE)
#' @param verbose verbosity level (0=silent, 1=progress, 2=detailed)
#' @param cores number of cores for parallelization. If NULL, auto-detected
#'
#' @return a list with one data frame per sample containing differential expression results.
#' Each data frame includes columns: sample, gene, cluster_1, cluster_2 (if pairwise),
#' avg_log2fc, mm_p_val/ttest_p_val/wilcox_p_val, adj_p_val, exp_p_val, exp_adj_p_val,
#' comments (indicating spatial model status: 'no_spatial_test', 'time_out',
#' 'no_convergence', or NA for successful spatial models)
#'
#' @export
#' @examples
#' \dontrun{
#' # Load TNBC test data
#' data_dir <- "tests/testthat/data/tnbc_bassiouni"
#' count_files <- list.files(data_dir, pattern='counts', full.names=TRUE)[1:2]
#' coord_files <- list.files(data_dir, pattern='mapping', full.names=TRUE)[1:2]
#' clin_file <- file.path(data_dir, "bassiouni_clinical.csv")
#' 
#' # Create STlist
#' tnbc <- STlist(rnacounts=count_files, spotcoords=coord_files, samples=clin_file)
#' tnbc <- transform_data(tnbc)
#' 
#' # Run non-spatial DE (Wilcoxon)
#' tnbc <- STdiff(tnbc, test='wilcoxon', group='tissue_type')
#' 
#' # Extract and view results
#' res <- STdiff_compile_results(tnbc)
#' head(res)
#' 
#' # Alternative: Run non-spatial tests only
#' result = STdiff(x=stlist_obj, annot='cluster', test_type='mm', sp_topgenes=0)
#' }
#'
#' @keywords spatial transcriptomics differential expression spaMM
STdiff = function(x=NULL, samples=NULL, annot=NULL, w=NULL, k=NULL, deepSplit=NULL,
                  topgenes=5000, pval_thr=0.05, pval_adj='fdr', test_type='mm', sp_topgenes=0.2,
                  clusters=NULL, pairwise=FALSE, verbose=1L, cores=NULL){

  . = NULL
  zero_t = Sys.time()

  # Convert user inputs to expected types
  topgenes = as.integer(ceiling(topgenes))
  pval_thr = as.double(pval_thr)
  sp_topgenes = as.double(sp_topgenes)
  
  if(sp_topgenes < 0 || sp_topgenes > 1){ stop('sp_topgenes must be between 0 and 1.') }
  
  verbose = as.integer(verbose)
  if(!is.integer(verbose) || !(verbose %in% c(0L, 1L, 2L))){ verbose = 1L }

  if(!test_type %in% c('mm', 't_test', 'wilcoxon')){
    stop('test_type must be one of "mm", "t_test", or "wilcoxon".')
  }

  if(pairwise && !is.null(clusters) && length(clusters) < 2){
    stop('If pairwise tests requested, at least two clusters are required.')
  }

  if(round(topgenes, 0) <= 0){ stop('topgenes must be a positive integer.') }

  # Define samples
  if(is.null(samples)){
    samples = names(x@spatial_meta)
  } else {
    if(is.numeric(samples)){ samples = as.vector(na.omit(names(x@spatial_meta)[samples])) }
    else { samples = samples[samples %in% names(x@spatial_meta)] }
    if(length(samples) == 0 || !any(samples %in% names(x@spatial_meta))){
      stop('None of the requested samples are present in the STlist.')
    }
  }

  # Define annotation column
  if(is.null(annot)){
    if(!is.null(w) && !is.null(k)){
      if(length(w) != 1 || length(k) != 1){ stop('Please specify a single value for w and a single value for k.') }
      annot = paste0('stclust_spw', as.character(w))
      if(k == 'dtc'){
        if(is.null(deepSplit)){ stop('If k="dtc", then specify deepSplit.') }
        else if(is.logical(deepSplit)){ annot = paste0(annot, '_dspl', ifelse(deepSplit, 'True', 'False')) }
        else { annot = paste0(annot, '_dspl', deepSplit) }
      } else {
        if(!is.numeric(k)){ stop('Specify a valid k value.') }
        annot = paste0(annot, '_k', as.character(k))
      }
    } else { stop('If no specific annotation is specified, please specify both w and k (STclust parameters).') }
  } else {
    if(length(annot) == 0){ stop('If w and k are not specified, one annotation column from @spatial_meta should be specified.') }
    else if(length(annot) > 1){ stop('Only one annotation column from @spatial_meta can be tested at a time.') }
  }

  # Validate annotation exists in samples
  samples_tmp = samples
  for(i in samples){
    if(!(annot %in% colnames(x@spatial_meta[[i]]))){
      samples_tmp = grep(i, samples_tmp, value=TRUE, invert=TRUE)
      if(verbose){ cat(paste0('Skipping ', i, '. Annotation not available for this sample.\n')) }
    }
    if(length(samples_tmp) == 0){ stop('No samples left to test. Are the requested annotations/clusters present in at least one sample?') }
  }
  samples = samples_tmp

  # Run non-spatial tests
  if(verbose){
    cat('\nRunning STdiff...\n')
    cat(paste0('\tTest type: ', test_type, '\n'))
    cat(paste0('\tSamples: ', length(samples), '\n'))
    cat(paste0('\tAnnotation: ', annot, '\n'))
    cat(paste0('\tTop genes: ', topgenes, '\n'))
    cat(paste0('\tP-value threshold: ', pval_thr, '\n'))
    cat(paste0('\tSpatial models: ', ifelse(sp_topgenes > 0, 'YES', 'NO'), '\n'))
    cat('\n')
  }

  # Step 1: Run non-spatial differential expression tests
  prep = STdiff_run_nonspatial(
    x = x, samples = samples, annot = annot, topgenes = topgenes, pval_thr = pval_thr,
    pval_adj = pval_adj, test_type = test_type, clusters = clusters, pairwise = pairwise,
    verbose = verbose, cores = cores
  )

  # Step 2: Run spatial models if requested
  spatial = NULL
  if(test_type == 'mm' && sp_topgenes > 0){
    if(verbose){ cat('\tRunning spatial tests...\n') }
    spatial = STdiff_run_spatial(x = x, prep = prep, sp_topgenes = sp_topgenes, cores = cores, verbose = verbose)
  }

  # Step 3: Compile results
  result_de = STdiff_compile_results(prep = prep, spatial = spatial)

  # Print completion message
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){ cat(paste0('\nSTdiff completed in ', round(end_t, 2), ' min.\n')) }

  return(result_de)
}
