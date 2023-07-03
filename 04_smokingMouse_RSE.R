## ----download_data_biocfilecache-------------------------------------
## Load the container package for this type of data
library("SummarizedExperiment")

## Download and cache the file
library("BiocFileCache")
bfc <- BiocFileCache::BiocFileCache()
cached_rse_gene <- BiocFileCache::bfcrpath(
    x = bfc,
    "https://github.com/LieberInstitute/SPEAQeasyWorkshop2023/raw/devel/provisional_data/rse_gene_mouse_RNAseq_nic-smo.Rdata"
)

## Check the local path on our cache
cached_rse_gene

## Load the rse_gene object
load(cached_rse_gene, verbose = TRUE)

## General overview of the object
rse_gene


## ----explore_assays--------------------------------------------------
## Explore main assay (of raw counts)
assay(rse_gene)[1:3, 1:3] ## counts for first 3 genes and 3 samples
## Access the same raw data with assays()
assays(rse_gene)$counts[1:3, 1:3]
## Access lognorm counts
assays(rse_gene)$logcounts[1:3, 1:3]


## ----explore_colData-------------------------------------------------
## Data for first 3 samples and 5 variables
colData(rse_gene)[1:3, 1:5]


## ----explore_rowData-------------------------------------------------
## Data for first 3 genes and 5 variables
rowData(rse_gene)[1:3, 1:5]


## ----exercise1_data, eval=FALSE, echo=FALSE--------------------------
## ## Solution
## head(rse_gene_nic$flowcell)


## ----extract_nicotine_data-------------------------------------------
## Original dimensions of the data
dim(rse_gene)
rse_gene_nic <- rse_gene[, which(rse_gene$Expt == "Nicotine")]
## New dimensions
dim(rse_gene_nic)


## ----exercise2_data, eval=FALSE, echo=FALSE--------------------------
## ## Solution
## table(rse_gene_nic$Expt)


## ----exercise3_data, eval=FALSE, echo=FALSE--------------------------
## ## Solution
## table(rse_gene$Age)
## pup_samples <- rse_gene[, which(rse_gene$Age == "Pup")]
## table(pup_samples$Sex)

