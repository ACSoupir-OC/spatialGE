##
# STdiff: Differential expression analysis for spatial transcriptomics data
# 
# This file contains the public STdiff function (modular implementation)
# and the legacy STdiff_legacy function for reproducibility.
#
# @importFrom magrittr %>%
# @importFrom rlang :=
# @importFrom stats p.adjust
#

# Source helper functions - use absolute path for reproducibility
helper_path <- system.file('R', 'STdiff_helpers.R', package='spatialGE', mustWork=FALSE)
if(helper_path == '' || !file.exists(helper_path)){
  # If not installed, use relative path
  source('R/STdiff_helpers.R')
} else {
  source(helper_path)
}

# Source core modular functions
core_path <- system.file('R', 'STdiff_core.R', package='spatialGE', mustWork=FALSE)
if(core_path == '' || !file.exists(core_path)){
  source('R/STdiff_core.R')
} else {
  source(core_path)
}

# Legacy function is defined in STdiff_legacy.R
# Import it for reproducibility
legacy_path <- system.file('R', 'STdiff_legacy.R', package='spatialGE', mustWork=FALSE)
if(legacy_path == '' || !file.exists(legacy_path)){
  if(file.exists('R/STdiff_legacy.R')){
    source('R/STdiff_legacy.R')
  }
} else {
  source(legacy_path)
}