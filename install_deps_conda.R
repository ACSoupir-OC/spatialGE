#!/usr/bin/env Rscript

# Install spatialGE dependencies to conda environment
lib_path <- "/home/node/miniconda3/envs/spatialge-r/lib/R/library"
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Core CRAN packages for spatialGE
cran_packages <- c(
    "arrow", "concaveman", "data.table", "dynamicTreeCut",
    "dplyr", "fastcluster", "ggforce", "ggplot2", "ggrepel",
    "gstat", "hdf5r", "jpeg", "jsonlite", "khroma", "lifecycle", "magrittr",
    "Matrix", "MASS", "png", "RColorBrewer", "Rcpp", "RcppEigen", "RcppProgress",
    "readr", "readxl", "rlang", "scales", "sctransform", "sp",
    "spaMM", "spdep", "stringr", "tibble", "tidyr", "uwot", "wordspace",
    "testthat"
)

# Install CRAN packages to conda lib
cat("Installing CRAN packages to conda environment...\n")
for (pkg in cran_packages) {
    if (!require(pkg, character.only = TRUE, lib.loc = lib_path, quietly = TRUE)) {
        cat("Installing", pkg, "...\n")
        tryCatch({
            install.packages(pkg, lib = lib_path, dependencies = TRUE, quiet = TRUE)
        }, error = function(e) {
            cat("Failed to install", pkg, ":", e$message, "\n")
        })
    }
}

# Bioconductor packages
bioc_packages <- c("BiocManager", "ComplexHeatmap", "DelayedArray", "DelayedMatrixStats", 
                   "EBImage", "GSVA", "BiocParallel", "Biobase", "SummarizedExperiment")

cat("\nInstalling Bioconductor packages...\n")
BiocManager::install(bioc_packages, lib = lib_path, ask = FALSE, quiet = TRUE, update = FALSE)

cat("\nAll packages installed!\n")
