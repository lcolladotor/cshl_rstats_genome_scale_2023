## ----'start', message=FALSE----
## Load recount3 R package
library("recount3")


## ----'quick_example'--------
## Lets download all the available projects
human_projects <- available_projects()

## Find your project of interest. Here we'll use
## SRP009615 as an example
proj_info <- subset(
    human_projects,
    project == "SRP009615" & project_type == "data_sources"
)
## Build a RangedSummarizedExperiment (RSE) object
## with the information at the gene level
rse_gene_SRP009615 <- create_rse(proj_info)
## Explore the resulting object
rse_gene_SRP009615

## How large is it?
lobstr::obj_size(rse_gene_SRP009615)


## ----"interactive_display", eval = FALSE----
## ## Explore available human projects interactively
## proj_info_interactive <- interactiveDisplayBase::display(human_projects)
## ## Choose only 1 row in the table, then click on "send".
## 
## ## Lets double check that you indeed selected only 1 row in the table
## stopifnot(nrow(proj_info_interactive) == 1)
## ## Now we can build the RSE object
## rse_gene_interactive <- create_rse(proj_info_interactive)


## ----"tranform_counts"------
## We'll compute read counts, which is what most downstream software
## uses.
## For other types of transformations such as RPKM and TPM, use
## transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)


## ----"expand_attributes"----
## Lets make it easier to use the information available for this study
## that was provided by the original authors of the study.
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)
colData(rse_gene_SRP009615)[
    ,
    grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]

