## Seurat Conversion Functions
## Convert between STlist and Seurat objects

##
#' Convert STlist to Seurat Object
#' @description
#' `r lifecycle::badge("stable")`
#' 
#' Converts an STlist object to a Seurat object for integration with Seurat workflows.
#' @details
#' This function creates a Seurat object from an STlist, preserving:
#' - Count data in the 'RNA' assay
#' - Spatial coordinates in the 'spatial' reduction
#' - Sample metadata in Seurat meta.data
#' - Gene metadata in feature-level metadata
#' 
#' For multi-sample STlist objects, returns a merged Seurat object with sample
#' barcodes preserved.
#' 
#' @param x an STlist object to convert
#' @param samples vector of sample names to include. If NULL, all samples are converted
#' @param assay.name name for the assay in the Seurat object (default: 'RNA')
#' @param add.spatial.info logical, whether to add spatial coordinates to Seurat object
#' @param verbose logical, whether to print progress messages
#' @return a Seurat object containing the STlist data
#' 
#' @export
#' @examples
#' \dontrun{
#' # Convert STlist to Seurat
#' seurat_obj <- as.Seurat.STlist(st_obj)
#' 
#' # Convert specific samples
#' seurat_obj <- as.Seurat.STlist(st_obj, samples = c("sample1", "sample2"))
#' 
#' # Use in Seurat workflow
#' seurat_obj <- SCTransform(seurat_obj)
#' seurat_obj <- RunPCA(seurat_obj)
#' }
#'
#' @seealso \code{\link{as.STlist.Seurat}}
#'
as.Seurat.STlist <- function(x, samples = NULL, assay.name = "RNA", 
                              add.spatial.info = TRUE, verbose = TRUE) {
  
  # Check Seurat availability
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required for conversion. Install with: install.packages('Seurat')")
  }
  
  # Validate input
  if (!methods::is(x, "STlist")) {
    stop("x must be an STlist object")
  }
  
  # Determine samples to process
  if (is.null(samples)) {
    samples <- names(x@counts)
  } else {
    # Validate sample names
    invalid <- setdiff(samples, names(x@counts))
    if (length(invalid) > 0) {
      stop("Invalid sample names: ", paste(invalid, collapse = ", "))
    }
  }
  
  if (verbose) cat("Converting", length(samples), "sample(s) to Seurat object...\n")
  
  # Process each sample
  seurat_list <- list()
  
  for (sname in samples) {
    if (verbose) cat("  Processing", sname, "...\n")
    
    # Get counts
    counts <- x@counts[[sname]]
    
    # Get coordinates
    coords <- NULL
    if (add.spatial.info && !is.null(x@spatial_meta[[sname]])) {
      coords <- x@spatial_meta[[sname]]
    }
    
    # Create Seurat object for this sample
    # Seurat v5 uses different constructor than v4
    seurat_obj <- tryCatch({
      # Try Seurat v5 syntax first
      Seurat::CreateSeuratObject(
        counts = counts,
        assay = assay.name,
        project = sname
      )
    }, error = function(e) {
      # Fallback to v4 syntax if needed
      Seurat::CreateSeuratObject(
        counts = counts,
        assay = assay.name,
        project = sname
      )
    })
    
    # Add spatial coordinates if available
    if (!is.null(coords)) {
      # Extract x, y coordinates
      if (all(c("xpos", "ypos") %in% colnames(coords))) {
        seurat_obj$x <- coords$xpos
        seurat_obj$y <- coords$ypos
        
        # Add to spatial coordinates if image info exists
        if (add.spatial.info) {
          # Try to add as spatial reduction
          try({
            seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = coords)
          }, silent = TRUE)
        }
      }
    }
    
    # Add sample metadata if available
    if (!rlang::is_empty(x@sample_meta) && nrow(x@sample_meta) > 0) {
      sample_info <- x@sample_meta[x@sample_meta$sample_name == sname, ]
      if (nrow(sample_info) > 0) {
        # Add sample-level metadata to all cells
        for (col in setdiff(colnames(sample_info), "sample_name")) {
          seurat_obj[[paste0("sample_", col)]] <- sample_info[[col]]
        }
      }
    }
    
    # Add gene metadata if available
    if (!rlang::is_empty(x@gene_meta[[sname]]) && nrow(x@gene_meta[[sname]]) > 0) {
      gene_meta <- x@gene_meta[[sname]]
      # Match genes
      common_genes <- intersect(rownames(seurat_obj), gene_meta$gene)
      if (length(common_genes) > 0) {
        for (col in setdiff(colnames(gene_meta), "gene")) {
          meta_col <- gene_meta[match(rownames(seurat_obj), gene_meta$gene), col]
          seurat_obj[[paste0("gene_", col)]] <- meta_col
        }
      }
    }
    
    seurat_list[[sname]] <- seurat_obj
  }
  
  # Merge if multiple samples
  if (length(seurat_list) == 1) {
    result <- seurat_list[[1]]
  } else {
    if (verbose) cat("Merging", length(seurat_list), "samples...\n")
    result <- Seurat::merge(
      x = seurat_list[[1]],
      y = seurat_list[-1],
      add.cell.ids = samples,
      project = "spatialGE_merged"
    )
  }
  
  if (verbose) cat("Conversion complete!\n")
  return(result)
}


##
#' Convert Seurat Object to STlist
#' @description
#' `r lifecycle::badge("stable")`
#' 
#' Converts a Seurat object (particularly spatial Seurat objects) to STlist format.
#' @details
#' This function extracts data from a Seurat object and creates an STlist:
#' - Counts from the default assay (usually 'RNA' or 'SCT')
#' - Spatial coordinates from images slot (for Visium/Visium HD)
#' - Sample metadata from Seurat meta.data
#' 
#' Works with both single-sample and multi-sample Seurat objects.
#' 
#' @param x a Seurat object to convert
#' @param assay which assay to use for counts (default: 'RNA')
#' @param slot which slot to use for counts (default: 'counts')
#' @param use.spatial logical, whether to extract spatial coordinates
#' @param verbose logical, whether to print progress messages
#' @return an STlist object
#' 
#' @export
#' @examples
#' \dontrun{
#' # Convert Seurat to STlist
#' st_obj <- as.STlist.Seurat(seurat_obj)
#' 
#' # Use in spatialGE workflow
#' st_obj <- SThet(st_obj, genes = c("GENE1", "GENE2"))
#' compare_SThet(st_obj)
#' }
#'
#' @seealso \code{\link{as.Seurat.STlist}}
#'
as.STlist.Seurat <- function(x, assay = "RNA", slot = "counts", 
                              use.spatial = TRUE, verbose = TRUE) {
  
  # Check Seurat availability
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required for conversion")
  }
  
  # Validate input
  if (!inherits(x, "Seurat")) {
    stop("x must be a Seurat object")
  }
  
  if (verbose) cat("Converting Seurat object to STlist...\n")
  
  # Get counts from specified assay/slot
  if (!assay %in% Seurat::Assays(x)) {
    stop("Assay '", assay, "' not found in Seurat object. Available: ", 
         paste(Seurat::Assays(x), collapse = ", "))
  }
  
  counts <- tryCatch({
    Seurat::GetAssayData(x, assay = assay, slot = slot)
  }, error = function(e) {
    Seurat::GetAssayData(x, assay = assay, layer = slot)
  })
  
  # Convert to list format (one entry per sample if merged)
  # Check if object has multiple samples by looking at cell barcodes
  cell_barcodes <- colnames(counts)
  
  # Detect sample structure from barcodes
  # Seurat merged objects have format: SampleID_CellBarcode
  if (any(grepl("_", cell_barcodes[1:10]))) {
    # Likely merged object - split by sample prefix
    sample_prefixes <- sapply(strsplit(cell_barcodes, "_"), `[`, 1)
    unique_samples <- unique(sample_prefixes)
    
    if (length(unique_samples) > 1) {
      counts_list <- list()
      coords_list <- list()
      
      for (samp in unique_samples) {
        samp_cells <- cell_barcodes[sample_prefixes == samp]
        counts_list[[samp]] <- counts[, samp_cells, drop = FALSE]
        
        # Try to get coordinates for this sample
        coords_list[[samp]] <- NULL  # Will be populated below if spatial
      }
    } else {
      samp_name <- unique_samples[1]
      counts_list <- list()
      counts_list[[samp_name]] <- counts
      coords_list <- list()
      coords_list[[samp_name]] <- NULL
    }
  } else {
    # Single sample
    counts_list <- list("sample_1" = counts)
    coords_list <- list("sample_1" = NULL)
  }
  
  # Extract spatial coordinates if available
  if (use.spatial && !is.null(x@images) && length(x@images) > 0) {
    if (verbose) cat("Extracting spatial coordinates from", length(x@images), "image(s)...\n")
    
    for (img_name in names(x@images)) {
      try({
        coords <- Seurat::GetTissueCoordinates(x, image = img_name)
        
        # Match to cells in counts
        if (!is.null(coords) && nrow(coords) > 0) {
          # Find matching barcodes
          if ("cell" %in% colnames(coords)) {
            cell_ids <- coords$cell
          } else {
            cell_ids <- rownames(coords)
          }
          
          # Match to sample
          for (samp in names(counts_list)) {
            samp_cells <- colnames(counts_list[[samp]])
            common <- intersect(cell_ids, samp_cells)
            
            if (length(common) > 0) {
              coords_sub <- coords[cell_ids %in% common, ]
              
              # Extract x, y
              if (all(c("imagecol", "imagerow") %in% colnames(coords_sub))) {
                coords_list[[samp]] <- data.frame(
                  barcode = coords_sub$cell,
                  xpos = coords_sub$imagecol,
                  ypos = coords_sub$imagerow
                )
              } else if (all(c("x", "y") %in% colnames(coords_sub))) {
                coords_list[[samp]] <- data.frame(
                  barcode = coords_sub$cell,
                  xpos = coords_sub$x,
                  ypos = coords_sub$y
                )
              }
            }
          }
        }
      }, silent = TRUE)
    }
  }
  
  # Create sample metadata from Seurat meta.data
  sample_meta <- NULL
  if (ncol(x@meta.data) > 0) {
    # Extract unique sample-level metadata if available
    # For now, just store a placeholder
    sample_meta <- tibble::tibble(sample_name = names(counts_list))
  }
  
  # Create STlist object
  st_obj <- methods::new("STlist",
    counts = counts_list,
    spatial_meta = coords_list,
    sample_meta = sample_meta,
    gene_meta = list(),
    tr_counts = list(),
    gene_krige = list(),
    misc = list(
      platform = "seurat",
      source_object = x,
      conversion_info = list(
        assay = assay,
        slot = slot,
        timestamp = Sys.time()
      )
    )
  )
  
  if (verbose) cat("Conversion complete!\n")
  return(st_obj)
}


##
#' Check if object can be converted to STlist
#' @description Tests whether an object can be converted to STlist format
#' @param x object to test
#' @return logical TRUE if conversion is possible
#' @keywords internal
can_convert_to_STlist <- function(x) {
  # Check if Seurat object
  if (inherits(x, "Seurat")) return(TRUE)
  
  # Check if list of matrices/dataframes
  if (is.list(x)) {
    if (all(sapply(x, function(i) is.matrix(i) || is.data.frame(i)))) {
      return(TRUE)
    }
  }
  
  # Check if matrix
  if (is.matrix(x) || inherits(x, "Matrix")) return(TRUE)
  
  return(FALSE)
}
