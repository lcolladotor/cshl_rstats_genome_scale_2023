## ----download_data_biocfilecache_speaqeasy_example---------------
## Load the container package for this type of data
library("SummarizedExperiment")

## Download and cache the file
library("BiocFileCache")
bfc <- BiocFileCache::BiocFileCache()
cached_rse_gene_example <- BiocFileCache::bfcrpath(
    x = bfc,
    "https://github.com/LieberInstitute/SPEAQeasy-example/raw/master/rse_speaqeasy.RData"
)

## Check the local path on our cache
cached_rse_gene_example

## Load the rse_gene object
load(cached_rse_gene_example, verbose = TRUE)

## General overview of the object
rse_gene

## We can check how big the object is with lobstr
lobstr::obj_size(rse_gene)


## ----------------------------------------------------------------
class(rse_gene$trimmed)

## logical vectors can take 2 values (plus the third `NA` if it's missing)

