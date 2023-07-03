## ------------------------------------------------------------
## Lets build a simple SummarizedExperiment object following information
## from the documentation
library("SummarizedExperiment")
## ?SummarizedExperiment

## Adapted from the official documentation:

## First we create the data pieces that we'll use to build our
## SummarizedExperiment object. In this case, we'll have 200 genes
## measured in 6 samples.
nrows <- 200
ncols <- 6

## Let's make up some count numbers at random
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)

## Then some basic infomratino for our genes
rowRanges <- GRanges(
    rep(c("chr1", "chr2"), c(50, 150)),
    IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
    strand = sample(c("+", "-"), 200, TRUE),
    feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))

## Next, we create some information about samples
colData <- DataFrame(
    Treatment = rep(c("ChIP", "Input"), 3),
    row.names = LETTERS[1:6]
)

## Finally we put all these pieces together in a single R object
rse <- SummarizedExperiment(
    assays = SimpleList(counts = counts),
    rowRanges = rowRanges,
    colData = colData
)

## Overview
rse


## ----isee_basic, eval = FALSE--------------------------------
## ## Let's explore the `rse` object interactively
## library("iSEE")
## iSEE::iSEE(rse)


## ----download_sce_layer--------------------------------------
## Lets get some data using spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")
sce_layer

## We can check how big the object is with lobstr
lobstr::obj_size(sce_layer)


## ----explore_sce_layer, eval = FALSE-------------------------
## iSEE::iSEE(sce_layer)

