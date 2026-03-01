##
# STclust: Spatial clustering (modular implementation)
# 
# This file contains the public STclust function (modular implementation)
# and the legacy STclust_legacy function for reproducibility.
#
# @importFrom magrittr %>%
# @importFrom methods as is new

# Source helper functions
source('R/STclust_helpers.R')

# Source seurat helpers (for Seurat_FindVariableFeatures)
source('R/seurat_helpers.R')

# Source core modular functions
source('R/STclust_core.R')

# Legacy function is defined in STclust_legacy.R
# Import it for reproducibility
if(file.exists('R/STclust_legacy.R')){
  source('R/STclust_legacy.R')
}