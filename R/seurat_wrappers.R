## Seurat Wrapper Functions
## Convenience functions for common Seurat-spatialGE workflows

##
#' Run spatialGE analysis on Seurat object
#' @description
#' `r lifecycle::badge("stable")`
#' 
#' Convenience wrapper that converts Seurat to STlist and runs spatial analysis.
#' @details
#' This function provides a seamless workflow from Seurat to spatialGE:
#' 1. Converts Seurat object to STlist
#' 2. Runs specified spatial analysis (SThet, STclust, STdiff, etc.)
#' 3. Optionally converts results back to Seurat format
#' 
#' @param x a Seurat object (spatial or non-spatial)
#' @param analysis type of analysis: 'SThet', 'STclust', 'STdiff', 'STenrich', 'STgradient'
#' @param genes vector of gene names to analyze (required for SThet, STgradient)
#' @param assay which assay to use for counts (default: 'RNA')
#' @param verbose logical, whether to print progress
#' @param ... additional arguments passed to the analysis function
#' @return depends on analysis type:
#'   - SThet: STlist with gene_meta updated
#'   - STclust: STlist with spatial_meta updated (cluster assignments)
#'   - STdiff: list of differential expression results
#'   - STenrich: list of enrichment results
#'   - STgradient: list of gradient analysis results
#' 
#' @export
#' @examples
#' \dontrun{
#' # Run heterogeneity analysis
#' result <- spatialGE_from_seurat(seurat_obj, analysis = 'SThet', 
#'                                  genes = c('GENE1', 'GENE2'))
#' 
#' # Run clustering
#' result <- spatialGE_from_seurat(seurat_obj, analysis = 'STclust')
#' 
#' # Run differential expression
#' result <- spatialGE_from_seurat(seurat_obj, analysis = 'STdiff',
#'                                  annot = 'cell_type')
#' }
#'
#' @seealso \code{\link{as.STlist.Seurat}}, \code{\link{as.Seurat.STlist}}
#'
spatialGE_from_seurat <- function(x, analysis = 'SThet', genes = NULL,
                                   assay = 'RNA', verbose = TRUE, ...) {
  
  # Convert to STlist
  st_obj <- as.STlist.Seurat(x, assay = assay, verbose = verbose)
  
  # Run analysis
  if (verbose) cat("Running", analysis, "analysis...\n")
  
  result <- switch(analysis,
    SThet = {
      if (is.null(genes)) {
        stop("genes parameter required for SThet analysis")
      }
      SThet(st_obj, genes = genes, verbose = verbose, ...)
    },
    STclust = {
      STclust(st_obj, verbose = verbose, ...)
    },
    STdiff = {
      STdiff(st_obj, verbose = verbose, ...)
    },
    STenrich = {
      STenrich(st_obj, verbose = verbose, ...)
    },
    STgradient = {
      STgradient(st_obj, verbose = verbose, ...)
    },
    stop("Unknown analysis type: ", analysis, 
         ". Valid options: SThet, STclust, STdiff, STenrich, STgradient")
  )
  
  return(result)
}


##
#' Add spatialGE results to Seurat object
#' @description
#' `r lifecycle::badge("stable")`
#' 
#' Adds spatial analysis results (e.g., cluster assignments, Moran's I) to a Seurat object's metadata.
#' @details
#' After running spatialGE analysis on an STlist, this function transfers the results
#' back to the original Seurat object for integrated visualization and downstream analysis.
#' 
#' @param seurat_obj original Seurat object
#' @param st_obj STlist object with analysis results
#' @param result_type type of results to add: 'clusters', 'moran', 'geary', 'custom'
#' @param custom.metadata named list of metadata columns to add
#' @param verbose logical, whether to print progress
#' @return Seurat object with added metadata
#' 
#' @export
#' @examples
#' \dontrun{
#' # Convert to STlist and run clustering
#' st_obj <- as.STlist.Seurat(seurat_obj)
#' st_obj <- STclust(st_obj)
#' 
#' # Add cluster assignments back to Seurat
#' seurat_obj <- add_spatialGE_to_seurat(seurat_obj, st_obj, result_type = 'clusters')
#' 
#' # Visualize in Seurat
#' Seurat::DimPlot(seurat_obj, group.by = 'spatialGE_cluster')
#' }
#'
#' @seealso \code{\link{as.STlist.Seurat}}, \code{\link{as.Seurat.STlist}}
#'
add_spatialGE_to_seurat <- function(seurat_obj, st_obj, result_type = 'clusters',
                                     custom.metadata = NULL, verbose = TRUE) {
  
  # Validate inputs
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  if (!methods::is(st_obj, "STlist")) {
    stop("st_obj must be an STlist object")
  }
  
  if (verbose) cat("Adding spatialGE results to Seurat object...\n")
  
  # Get cell barcodes from Seurat object
  seurat_cells <- colnames(seurat_obj)
  
  # Extract metadata from STlist
  metadata_to_add <- list()
  
  if (result_type == 'clusters') {
    # Extract cluster assignments from spatial_meta
    for (samp in names(st_obj@spatial_meta)) {
      samp_meta <- st_obj@spatial_meta[[samp]]
      if (!is.null(samp_meta) && 'cluster' %in% colnames(samp_meta)) {
        # Match barcodes
        if ('barcode' %in% colnames(samp_meta)) {
          cluster_map <- setNames(samp_meta$cluster, samp_meta$barcode)
          metadata_to_add[[samp]] <- cluster_map
        }
      }
    }
    
    # Add to Seurat
    if (length(metadata_to_add) > 0) {
      # Merge all samples
      all_clusters <- unlist(metadata_to_add)
      names(all_clusters) <- gsub("^[^_]*_", "", names(all_clusters))  # Remove sample prefix if needed
      
      # Match to Seurat cells
      matched <- all_clusters[names(all_clusters) %in% seurat_cells]
      
      if (length(matched) > 0) {
        seurat_obj$spatialGE_cluster <- NA
        seurat_obj$spatialGE_cluster[names(matched)] <- matched
        if (verbose) cat("  Added cluster assignments for", length(matched), "cells\n")
      }
    }
    
  } else if (result_type == 'moran' || result_type == 'geary') {
    # Extract gene-level statistics from gene_meta
    stat_col <- ifelse(result_type == 'moran', 'moran.I', 'geary.C')
    
    for (samp in names(st_obj@gene_meta)) {
      samp_meta <- st_obj@gene_meta[[samp]]
      if (!is.null(samp_meta) && stat_col %in% colnames(samp_meta)) {
        gene_stats <- setNames(samp_meta[[stat_col]], samp_meta$gene)
        metadata_to_add[[samp]] <- gene_stats
      }
    }
    
    # Add as feature metadata
    if (length(metadata_to_add) > 0) {
      all_stats <- unlist(metadata_to_add)
      seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = data.frame(all_stats), col.name = paste0('spatialGE_', result_type))
      if (verbose) cat("  Added", result_type, "statistics for", length(all_stats), "genes\n")
    }
    
  } else if (result_type == 'custom' && !is.null(custom.metadata)) {
    # Add custom metadata
    for (col_name in names(custom.metadata)) {
      seurat_obj[[col_name]] <- custom.metadata[[col_name]]
      if (verbose) cat("  Added custom metadata:", col_name, "\n")
    }
  }
  
  if (verbose) cat("Done!\n")
  return(seurat_obj)
}


##
#' Create Seurat-compatible gene sets from spatialGE results
#' @description
#' `r lifecycle::badge("stable")`
#' 
#' Converts STenrich or STgradient results to gene set format for Seurat GSEA.
#' @details
#' This function extracts significant gene sets from spatialGE enrichment analysis
#' and formats them for use with Seurat's AddModuleScore or similar functions.
#' 
#' @param st_result results from STenrich or STgradient analysis
#' @param pval.thr p-value threshold for significance (default: 0.05)
#' @param return.type 'list' (gene sets) or 'data.frame' (summary table)
#' @return gene sets in Seurat-compatible format
#' 
#' @export
#' @examples
#' \dontrun{
#' # Run enrichment analysis
#' enrich_result <- STenrich(st_obj, gene_sets = pathway_list)
#' 
#' # Convert to Seurat gene sets
#' gene_sets <- spatialGE_to_seurat_genesets(enrich_result, pval.thr = 0.05)
#' 
#' # Add module scores
#' seurat_obj <- Seurat::AddModuleScore(seurat_obj, features = gene_sets)
#' }
#'
#' @seealso \code{\link{STenrich}}, \code{\link{STgradient}}
#'
spatialGE_to_seurat_genesets <- function(st_result, pval.thr = 0.05, 
                                          return.type = 'list') {
  
  # Validate input
  if (!is.list(st_result)) {
    stop("st_result must be a list (output from STenrich or STgradient)")
  }
  
  # Extract significant gene sets
  significant_sets <- list()
  
  for (samp in names(st_result)) {
    samp_result <- st_result[[samp]]
    
    if (is.data.frame(samp_result) && 'pval' %in% colnames(samp_result)) {
      sig_genesets <- samp_result[samp_result$pval < pval.thr, ]
      
      if (nrow(sig_genesets) > 0) {
        if ('geneset' %in% colnames(sig_genesets)) {
          # Gene set names available
          for (i in 1:nrow(sig_genesets)) {
            gs_name <- paste0(samp, "_", sig_genesets$geneset[i])
            significant_sets[[gs_name]] <- sig_genesets$genes[i]
          }
        }
      }
    }
  }
  
  if (return.type == 'list') {
    return(significant_sets)
  } else if (return.type == 'data.frame') {
    # Create summary table
    summary_df <- data.frame(
      geneset = names(significant_sets),
      n_genes = sapply(significant_sets, length),
      stringsAsFactors = FALSE
    )
    return(summary_df)
  }
  
  return(significant_sets)
}
