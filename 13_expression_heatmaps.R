## ---- include = FALSE--------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)


## ----vignetteSetup_expheatmap, echo=FALSE, message=FALSE, warning = FALSE----
## For links
library(BiocStyle)

## Bib setup
library(RefManageR)

## Write bibliography information
bib <- c(
    smokingMouse = citation("smokingMouse")[1],
    SummarizedExperiment = citation("SummarizedExperiment")[1],
    ComplexHeatmap = citation("ComplexHeatmap")[1]
)

options(max.print = 50)


## ---- echo=FALSE-------------------
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(data(airway, package = "airway"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))


## ----------------------------------
library("SummarizedExperiment")
library("ComplexHeatmap")
library("circlize")


## ----download_data_biocfilecache_expheatmap----
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


## ----------------------------------
## Extract genes and samples of interest
rse_gene_pup_nic <- rse_gene[
    rowData(rse_gene)$DE_in_pup_brain_nicotine == TRUE,
    rse_gene$Age == "Pup" & rse_gene$Expt == "Nicotine"
]

## Extract logcounts and add name columns
logs_pup_nic <- assay(rse_gene_pup_nic, 2)
colnames(logs_pup_nic) <- paste0("Pup_", 1:dim(rse_gene_pup_nic)[2])

## Create function to remove technical variables' contributions,
## this is from jaffelab package (https://github.com/LieberInstitute/jaffelab)
cleaningY <- function(y, mod, P) {
    stopifnot(P <= ncol(mod))
    stopifnot(`Input matrix is not full rank` = qr(mod)$rank ==
        ncol(mod))
    Hat <- solve(t(mod) %*% mod) %*% t(mod)
    ty <- t(y)
    ty[is.na(ty)] <- 0
    beta <- (Hat %*% ty)
    cleany <- y - t(as.matrix(mod[, -c(seq_len(P))]) %*% beta[-seq_len(P), ])
    return(cleany)
}

## Remove contribution of technical variables
formula <- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr +
    overallMapRate + mitoRate
model <- model.matrix(formula, data = colData(rse_gene_pup_nic))
logs_pup_nic <- cleaningY(logs_pup_nic, model, P = 2)

## Center the data to make differences more evident
logs_pup_nic <- (logs_pup_nic - rowMeans(logs_pup_nic)) / rowSds(logs_pup_nic)


## ----------------------------------
## Prepare annotation for our heatmap
## For this heatmap I want to be able to see the Group to which each sample belongs
## as well as the Sex of the pup

top_ans <- HeatmapAnnotation(
    df = as.data.frame(colData(rse_gene_pup_nic)[, c("Sex", "Group")]),
    col = list(
        "Sex" = c("F" = "hotpink1", "M" = "dodgerblue"),
        "Group" = c("Control" = "brown2", "Experimental" = "deepskyblue3")
    )
)

## Also, I want to add the FDR associated to each gene
## Even though we do have that data, for this particular exercise we are gonna simulate it

rowData(rse_gene_pup_nic)$FDR <- sample(
    x = 0:5000,
    size = dim(rse_gene_pup_nic)[1],
    prob = c(rep(0.0012, 1001), rep(0.00005, 4000))
) / 10000

left_ans <- rowAnnotation(
    FDR = rowData(rse_gene_pup_nic)$FDR,
    col = list(FDR = colorRamp2(c(0, 0.049), c("#ecf39e", "#4f772d"))),
    annotation_legend_param = list(FDR = list(at = c(0, 0.01, 0.02, 0.03, 0.04, 0.05)))
)


## ---- out.width = "1000px"---------
## Finally, let's plot!
Heatmap(logs_pup_nic,
    name = "logcounts",
    show_row_names = FALSE,
    top_annotation = top_ans,
    left_annotation = left_ans,
    row_km = 2,
    column_km = 2,
    col = colorRamp2(c(-4, -0.0001, 00001, 4), c("darkblue", "lightblue", "lightsalmon", "darkred")),
    row_title = NULL,
    column_title = NULL,
    column_names_gp = gpar(fontsize = 7)
)


## ----------------------------------
## For a) we need to use he function anno_barplot() to generate the barplot and gpar() to fill the bars

top_ans <- HeatmapAnnotation(
    df = as.data.frame(colData(rse_gene_pup_nic)[, c("Sex", "Group")]),
    library_size = anno_barplot(colData(rse_gene_pup_nic)$sum, gp = gpar(fill = "#ffd500")),
    col = list(
        "Sex" = c("F" = "hotpink1", "M" = "dodgerblue"),
        "Group" = c("Control" = "brown2", "Experimental" = "deepskyblue3")
    )
)

## For b), I want to know if a gene is protein_coding or not_protein_coding.

## First, we can explore the rowData of our object to see if we have a column with
## that information

names(rowData(rse_gene_pup_nic))

## Maybe gene_type has what we need
unique(rowData(rse_gene_pup_nic)$gene_type)

## It does but not in the format we need it, so we need to create a new column to complete the task
## This columns is going to have protein_coding and not_protein_coding as elements
rowData(rse_gene_pup_nic)$pc_status <- "protein_coding"
rowData(rse_gene_pup_nic)$pc_status[rowData(rse_gene_pup_nic)$gene_type != "protein_coding"] <- "not_protein_coding"

## Now we can add the information and specify the color
left_ans <- rowAnnotation(
    FDR = sample(x = 0:5000, size = dim(rse_gene_pup_nic)[1], prob = c(rep(0.0012, 1001), rep(0.00005, 4000))) / 10000,
    pc_status = rowData(rse_gene_pup_nic)$pc_status,
    col = list(FDR = colorRamp2(c(0, 0.049), c("#ecf39e", "#4f772d")), pc_status = c("not_protein_coding" = "#ffe5d9", "protein_coding" = "#9d8189")),
    annotation_legend_param = list(FDR = list(at = c(0, 0.01, 0.02, 0.03, 0.04, 0.05)))
)


## ---- out.width = "1000px"---------
Heatmap(logs_pup_nic,
    name = " ",
    show_row_names = FALSE,
    top_annotation = top_ans,
    left_annotation = left_ans,
    row_km = 2,
    column_km = 2,
    col = colorRamp2(c(-4, -0.0001, 00001, 4), c("darkblue", "lightblue", "lightsalmon", "darkred")),
    row_title = NULL,
    column_title = NULL,
    column_names_gp = gpar(fontsize = 7)
)

