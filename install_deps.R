#!/usr/bin/env Rscript

# Install spatialGE dependencies
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Core CRAN packages
cran_packages <- c(
    "arrow", "BiocManager", "concaveman", "data.table", "dynamicTreeCut",
    "dplyr", "fastcluster", "ggforce", "ggplot2", "ggpolypath", "ggrepel",
    "gstat", "hdf5r", "jpeg", "jsonlite", "khroma", "lifecycle", "magrittr",
    "Matrix", "MASS", "png", "RColorBrewer", "Rcpp", "RcppEigen", "RcppProgress",
    "readr", "readxl", "rlang", "scales", "sctransform", "sfsmisc", "sp",
    "spaMM", "spdep", "stringr", "tibble", "tidyr", "uwot", "wordspace",
    "testthat", "curl", "geoR", "ggpubr", "httr", "janitor", "kableExtra",
    "knitr", "msigdbr", "progress", "rmarkdown", "rpart", "R.utils", "spatstat"
)

# Install CRAN packages
cat("Installing CRAN packages...\n")
install.packages(cran_packages, dependencies = TRUE, quiet = TRUE)

# Bioconductor packages
bioc_packages <- c("ComplexHeatmap", "DelayedArray", "DelayedMatrixStats", 
                   "EBImage", "GSVA", "SeuratObject", "scSpatialSIM")

cat("Installing Bioconductor packages...\n")
BiocManager::install(bioc_packages, ask = FALSE, quiet = TRUE)

cat("All packages installed!\n")
