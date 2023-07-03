## ----install, eval = FALSE-----------------------------------
## ## For installing Bioconductor packages
## if (!requireNamespace("BiocManager", quietly = TRUE)) {
##     install.packages("BiocManager")
## }
## 
## ## Install required packages
## BiocManager::install(
##     c(
##         "usethis", ## Utilities
##         "BiocFileCache",
##         "RefManageR",
##         "here",
##         "Hmisc",
##         "biocthis",
##         "lobstr",
##         "postcards",
##         "scater",
##         "sessioninfo",
##         "smokingMouse",
##         "stringr",
##         "SummarizedExperiment", ## Main containers / vis
##         "iSEE",
##         "edgeR", ## RNA-seq
##         "ExploreModelMatrix",
##         "limma",
##         "recount3",
##         "rlang",
##         "pheatmap", ## Visualization
##         "ggplot2",
##         "ggrepel",
##         "patchwork",
##         "RColorBrewer",
##         "ComplexHeatmap",
##         "cowplot",
##         "spatialLIBD", ## Advanced
##         "variancePartition"
##     )
## )
## 
## ## Install smokingMouse, which is currently under review at Bioconductor
## ## at https://github.com/Bioconductor/Contributions/issues/3024.
## BiocManager::install("LieberInstitute/smokingMouse")


## ----session_packages, eval = TRUE, message = FALSE----------
## Load the package at the top of your script
library("sessioninfo")

## Utilities
library("BiocFileCache")
library("BiocStyle")
library("biocthis")
library("here")
library("lobstr")
library("postcards")
library("usethis")
library("sessioninfo")

## Data
library("smokingMouse")

## Main containers / vis
library("SummarizedExperiment")
library("iSEE")

## RNA-seq
library("edgeR")
library("ExploreModelMatrix")
library("limma")
library("recount3")

## QCA
library("scater")

## Variance Partition
library("variancePartition")

## Visualization: plots & text
library("ComplexHeatmap")
library("ggplot2")
library("patchwork")
library("pheatmap")
library("RColorBrewer")
library("Hmisc")
library("stringr")
library("cowplot")
library("rlang")
library("ggrepel")

## Spatial transcriptomics
library("spatialLIBD")


## ----session_info--------------------------------------------
## Reproducibility information
options(width = 120)
session_info()
proc.time()

