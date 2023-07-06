## ----"download_10x_data"----------------------------------------------------
## Download and save a local cache of the data provided by 10x Genomics
bfc <- BiocFileCache::BiocFileCache()
lymph.url <-
    paste0(
        "https://cf.10xgenomics.com/samples/spatial-exp/",
        "1.1.0/V1_Human_Lymph_Node/",
        c(
            "V1_Human_Lymph_Node_filtered_feature_bc_matrix.tar.gz",
            "V1_Human_Lymph_Node_spatial.tar.gz",
            "V1_Human_Lymph_Node_analysis.tar.gz"
        )
    )
lymph.data <- sapply(lymph.url, BiocFileCache::bfcrpath, x = bfc)


## ----"extract_files"--------------------------------------------------------
## Extract the files to a temporary location
## (they'll be deleted once you close your R session)
xx <- sapply(lymph.data, utils::untar, exdir = file.path(tempdir(), "outs"))
## The names are the URLs, which are long and thus too wide to be shown here,
## so we shorten them to only show the file name prior to displaying the
## utils::untar() output status
names(xx) <- basename(names(xx))
xx

## List the files we downloaded and extracted
## These files are typically SpaceRanger outputs
lymph.dirs <- file.path(
    tempdir(), "outs",
    c("filtered_feature_bc_matrix", "spatial", "raw_feature_bc_matrix", "analysis")
)
list.files(lymph.dirs)


## $ cd /dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A

## $ du -sh --apparent-size genes/genes.gtf

## 1.4G	genes/genes.gtf


## ----"use_gencode_gtf"------------------------------------------------------
## Download the Gencode v32 GTF file and cache it
gtf_cache <- BiocFileCache::bfcrpath(
    bfc,
    paste0(
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
        "release_32/gencode.v32.annotation.gtf.gz"
    )
)

## Show the GTF cache location
gtf_cache


## ----wrapper_functions------------------------------------------------------
## Import the data as a SpatialExperiment object using wrapper functions
## provided by spatialLIBD
library("spatialLIBD")
spe_wrapper <- read10xVisiumWrapper(
    samples = file.path(tempdir(), "outs"),
    sample_id = "lymph",
    type = "sparse", data = "filtered",
    images = c("lowres", "hires", "detected", "aligned"), load = TRUE,
    reference_gtf = gtf_cache
)

## Explore the resulting SpatialExperiment object
spe_wrapper

## Size of the data
lobstr::obj_size(spe_wrapper)


## ----"run_shiny_app_wrapper"------------------------------------------------
## Run our shiny app
if (interactive()) {
    vars <- colnames(colData(spe_wrapper))

    run_app(
        spe_wrapper,
        sce_layer = NULL,
        modeling_results = NULL,
        sig_genes = NULL,
        title = "spatialLIBD: human lymph node by 10x Genomics (made with wrapper)",
        spe_discrete_vars = c(vars[grep("10x_", vars)], "ManualAnnotation"),
        spe_continuous_vars = c("sum_umi", "sum_gene", "expr_chrM", "expr_chrM_ratio"),
        default_cluster = "10x_graphclust"
    )
}


## ---------------------------------------------------------------------------
## Basic spatial graph visualizing a discrete variable that we provided to
## the argument "clustervar". Here we chose to visualize
## "10x_kmeans_10_clusters" which contains the clustering results from the
## K-means algorithm when k = 10.
## We are plotting the one sample we have called "lymph". This is the same
## name we chose earlier when we ran spatialLIBD::read10xVisiumWrapper().
vis_clus(
    spe = spe_wrapper,
    sampleid = "lymph",
    clustervar = "10x_kmeans_10_clusters"
)

## Next we ignore the histology image and stop plotting it by setting
## "spatial = FALSE"
vis_clus(
    spe = spe_wrapper,
    sampleid = "lymph",
    clustervar = "10x_kmeans_10_clusters",
    spatial = FALSE
)

## Finally, we use the "colors" argument to specify our own colors. We use
## the general setNames() R function to create a named vector that has as
## values the colors and as names the same names we have for our "clustervar".
vis_clus(
    spe = spe_wrapper,
    sampleid = "lymph",
    clustervar = "10x_kmeans_10_clusters",
    spatial = FALSE,
    colors = setNames(Polychrome::palette36.colors(10), 1:10)
)

## Here we can see the vector we provided as input to "colors":
setNames(Polychrome::palette36.colors(10), 1:10)

