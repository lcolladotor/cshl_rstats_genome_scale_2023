## ----download_SRP045638----------------------------------------------
library("recount3")

human_projects <- available_projects()

rse_gene_SRP045638 <- create_rse(
    subset(
        human_projects,
        project == "SRP045638" & project_type == "data_sources"
    )
)
assay(rse_gene_SRP045638, "counts") <- compute_read_counts(rse_gene_SRP045638)


## ----describe_issue--------------------------------------------------
## Can you notice the problem with the sample information?
rse_gene_SRP045638$sra.sample_attributes[1:3]


## ----solve_issue-----------------------------------------------------
rse_gene_SRP045638$sra.sample_attributes <- gsub("dev_stage;;Fetal\\|", "", rse_gene_SRP045638$sra.sample_attributes)
rse_gene_SRP045638$sra.sample_attributes[1:3]


## ----attributes------------------------------------------------------
rse_gene_SRP045638 <- expand_sra_attributes(rse_gene_SRP045638)

colData(rse_gene_SRP045638)[
    ,
    grepl("^sra_attribute", colnames(colData(rse_gene_SRP045638)))
]


## ----re_cast---------------------------------------------------------
## Recast character vectors into numeric or factor ones
rse_gene_SRP045638$sra_attribute.age <- as.numeric(rse_gene_SRP045638$sra_attribute.age)
rse_gene_SRP045638$sra_attribute.disease <- factor(tolower(rse_gene_SRP045638$sra_attribute.disease))
rse_gene_SRP045638$sra_attribute.RIN <- as.numeric(rse_gene_SRP045638$sra_attribute.RIN)
rse_gene_SRP045638$sra_attribute.sex <- factor(rse_gene_SRP045638$sra_attribute.sex)

## Summary of our variables of interest
summary(as.data.frame(colData(rse_gene_SRP045638)[
    ,
    grepl("^sra_attribute.[age|disease|RIN|sex]", colnames(colData(rse_gene_SRP045638)))
]))


## ----new_variables---------------------------------------------------
## We'll want to look for differences between prenatal and postnatal samples
rse_gene_SRP045638$prenatal <- factor(ifelse(rse_gene_SRP045638$sra_attribute.age < 0, "prenatal", "postnatal"))
table(rse_gene_SRP045638$prenatal)

## http://rna.recount.bio/docs/quality-check-fields.html
rse_gene_SRP045638$assigned_gene_prop <- rse_gene_SRP045638$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP045638$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP045638$assigned_gene_prop)
with(colData(rse_gene_SRP045638), plot(assigned_gene_prop, sra_attribute.RIN))

## Hm... lets check if there is a difference between these two groups
with(colData(rse_gene_SRP045638), tapply(assigned_gene_prop, prenatal, summary))


## ----filter_rse------------------------------------------------------
## Lets save our full object for now in case we change our minds later on
rse_gene_SRP045638_unfiltered <- rse_gene_SRP045638

## Lets drop some bad samples. On a real analysis, you would likely use
## some statistical method for identifying outliers such as scuttle::isOutlier()
hist(rse_gene_SRP045638$assigned_gene_prop)
table(rse_gene_SRP045638$assigned_gene_prop < 0.3)
rse_gene_SRP045638 <- rse_gene_SRP045638[, rse_gene_SRP045638$assigned_gene_prop > 0.3]

## Lets compute the mean expression levels.
##
## Note: in a real analysis we would likely do this with RPKMs or CPMs instead
## of counts. That is, we would use one of the following options:
# edgeR::filterByExpr() https://bioconductor.org/packages/edgeR/ https://rdrr.io/bioc/edgeR/man/filterByExpr.html
# genefilter::genefilter() https://bioconductor.org/packages/genefilter/ https://rdrr.io/bioc/genefilter/man/genefilter.html
# jaffelab::expression_cutoff() http://research.libd.org/jaffelab/reference/expression_cutoff.html
#
gene_means <- rowMeans(assay(rse_gene_SRP045638, "counts"))
summary(gene_means)

## We can now drop genes with low expression levels
rse_gene_SRP045638 <- rse_gene_SRP045638[gene_means > 0.1, ]

## Final dimensions of our RSE object
dim(rse_gene_SRP045638)

## Percent of genes that we retained:
round(nrow(rse_gene_SRP045638) / nrow(rse_gene_SRP045638_unfiltered) * 100, 2)


## ----normalize-------------------------------------------------------
## Use edgeR::calcNormFactors() to address the composition bias
library("edgeR")
dge <- DGEList(
    counts = assay(rse_gene_SRP045638, "counts"),
    genes = rowData(rse_gene_SRP045638)
)
dge <- calcNormFactors(dge)


## ----explore_gene_prop_by_age----------------------------------------
library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP045638)), aes(y = assigned_gene_prop, x = prenatal)) +
    geom_boxplot() +
    theme_bw(base_size = 20) +
    ylab("Assigned Gene Prop") +
    xlab("Age Group")


## ----statiscal_model-------------------------------------------------
mod <- model.matrix(~ prenatal + sra_attribute.RIN + sra_attribute.sex + assigned_gene_prop,
    data = colData(rse_gene_SRP045638)
)
colnames(mod)


## ----run_limma-------------------------------------------------------
library("limma")
vGene <- voom(dge, mod, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
    eb_results,
    coef = 2,
    number = nrow(rse_gene_SRP045638),
    sort.by = "none"
)
dim(de_results)
head(de_results)

## Differentially expressed genes between pre and post natal with FDR < 5%
table(de_results$adj.P.Val < 0.05)

## We can now visualize the resulting differential expression results
plotMA(eb_results, coef = 2)

## We can also make a volcano plot
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)
de_results[de_results$gene_name %in% c("ZSCAN2", "VASH2", "KIAA0922"), ]


## ----pheatmap--------------------------------------------------------
## Extract the normalized expression values from our limma-voom result
## from earlier
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

## We can now build a table with information about our samples and
## then make the names a bit more friendly by making them easier to
## understand
df <- as.data.frame(colData(rse_gene_SRP045638)[, c("prenatal", "sra_attribute.RIN", "sra_attribute.sex")])
colnames(df) <- c("AgeGroup", "RIN", "Sex")

## Next, we can make a basic heatmap
library("pheatmap")
pheatmap(
    exprs_heatmap,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_col = df
)


## ----plot_mds--------------------------------------------------------
## For nicer colors
library("RColorBrewer")

## Mapping age groups into colors
col.group <- df$AgeGroup
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

## MDS by age groups
limma::plotMDS(vGene$E, labels = df$AgeGroup, col = col.group)

## Mapping Sex values into colors
col.sex <- df$Sex
levels(col.sex) <- brewer.pal(nlevels(col.sex), "Dark2")
col.sex <- as.character(col.sex)

## MDS by Sex
limma::plotMDS(vGene$E, labels = df$Sex, col = col.sex)

