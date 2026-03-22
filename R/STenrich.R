##
# STenrich: Spatial enrichment of gene expression sets
#
# This file contains the public STenrich function (modular implementation)
# and imports from STenrich_helpers.R for modular helper functions
#
# @importFrom magrittr %>%
# @importFrom stats p.adjust
# @importFrom tibble add_column column_to_rownames
# @importFrom dplyr mutate select
#

##
#' @title STenrich: Test for spatial enrichment of gene expression sets
#' @description
#' `r lifecycle::badge("stable")`
#' 
#' Test for spatial enrichment of gene expression sets in spatial transcriptomics data
#' @details The function performs a randomization test to assess if the sum of
#' distances between cells/spots with high expression of a gene set is lower than
#' the sum of distances among randomly selected cells/spots. The cells/spots are
#' considered as having high gene set expression if the average expression of genes in a
#' set is higher than the average expression plus `num_sds` times the standard deviation.
#' Control over the size of regions with high expression is provided by setting the
#' minimum number of cells/spots (`min_units`). This method is a modification of
#' the method devised by Hunter et al. 2021 (zebrafish melanoma study).
#'
#' @param x an STlist with transformed gene expression
#' @param samples a vector with sample names or indexes to run analysis
#' @param gene_sets a named list of gene sets to test. The names of the list should
#' identify the gene sets to be tested
#' @param score_type Controls how gene set expression is calculated. The options are
#' the average expression among genes in a set ('avg'), or a GSEA score ('gsva'). The
#' default is 'avg'
#' @param reps the number of random samples to be extracted. Default is 1000 replicates
#' @param annot name of the annotation within `x@spatial_meta` containing the spot/cell
#' categories. Needs to be used in conjunction with `domain`
#' @param domain the domain to restrict the analysis. Must exist within the spot/cell
#' categories included in the selected annotation (i.e., `annot`)
#' @param num_sds the number of standard deviations to set the minimum gene set
#' expression threshold. Default is one (1) standard deviation
#' @param min_units Minimum number of spots with high expression of a pathway for
#' that gene set to be considered in the analysis. Defaults to 20 spots or cells
#' @param min_genes the minimum number of genes of a gene set present in the data set
#' for that gene set to be included. Default is 5 genes
#' @param pval_adj_method the method for multiple comparison adjustment of p-values.
#' Options are the same as that of `p.adjust`. Default is 'BH'
#' @param seed the seed number for the selection of random samples. Default is 12345
#' @param cores the number of cores used during parallelization. If NULL (default),
#' the number of cores is defined automatically
#' @param verbose logical, whether to print text to console
#' @return a list of data frames with the results of the test
#'
#' @export
#'
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
#' # Define gene sets (Hallmark pathways)
#' gene_sets <- list(
#'   "HALLMARK_APOPTOSIS" = c("CASP3", "CASP7", "BAX"),
#'   "HALLMARK_GLYCOLYSIS" = c("HK2", "PFKP", "LDHA")
#' )
#'
#' # Run enrichment analysis
#' tnbc <- STenrich(tnbc, gene_sets=gene_sets, n_perm=100)
#' }
#'
#' @seealso STenrich_legacy for reproducibility with original implementation
#'
#' @importFrom magrittr %>%
#' @importFrom stats p.adjust
#' @importFrom tibble add_column column_to_rownames
#' @importFrom dplyr mutate select
#' @importFrom parallel mclapply count
STenrich = function(x, samples=NULL, gene_sets=NULL, score_type='avg', reps=1000,
                    annot=NULL, domain=NULL, num_sds=1, min_units=20, min_genes=5,
                    pval_adj_method='BH', seed=12345, cores=NULL, verbose=TRUE){

  # Record time
  zero_t = Sys.time()
  
  if(verbose){
    cat("Running STenrich...\n")
  }
  
  # Validate inputs
  validated = STenrich_validate_input(x, samples, gene_sets, score_type, annot, domain, num_sds, min_units, min_genes, pval_adj_method)
  x = validated$x
  samples = validated$samples
  gene_sets = validated$gene_sets
  score_type = validated$score_type
  annot = validated$annot
  domain = validated$domain
  num_sds = validated$num_sds
  min_units = validated$min_units
  min_genes = validated$min_genes
  pval_adj_method = validated$pval_adj_method
  
  # Call core workflow
  result = STenrich_core(x, samples, gene_sets, score_type, reps, annot, domain, 
                         num_sds, min_units, min_genes, pval_adj_method, seed, cores, verbose)
  
  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('STenrich completed in ', round(end_t, 2), ' min.\n'))
  }
  
  return(result)
}
