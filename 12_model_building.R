## ----download_data_biocfilecache_repeat_modeling---------------------
## Load the container package for this type of data
library("SummarizedExperiment")

## Download data
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

## Nicotine data
rse_gene_nic <- rse_gene[, which(rse_gene$Expt == "Nicotine")]

## Retain genes that passed filtering step
rse_gene_filt <- rse_gene_nic[rowData(rse_gene_nic)$retained_after_feature_filtering == TRUE, ]

## Separate data by Age
rse_gene_pups <- rse_gene_filt[, which(rse_gene_filt$Age == "Pup")]
rse_gene_adults <- rse_gene_filt[, which(rse_gene_filt$Age == "Adult")]

library("scuttle")
## Filter adult samples
outliers_library_size <- isOutlier(rse_gene_adults$sum, nmads = 3, type = "lower")
outliers_detected_num <- isOutlier(rse_gene_adults$detected, nmads = 3, type = "lower")
outliers_totalAssignedGene <- isOutlier(rse_gene_adults$totalAssignedGene, nmads = 3, type = "lower")
outliers_mito <- isOutlier(rse_gene_adults$mitoRate, nmads = 3, type = "higher")
outliers_rRNArate <- isOutlier(rse_gene_adults$rRNA_rate, nmads = 3, type = "higher")

not_outliers <- which(!(outliers_library_size | outliers_detected_num | outliers_totalAssignedGene | outliers_mito | outliers_rRNArate))
rse_gene_adults_qc <- rse_gene_adults[, not_outliers]

## Filter pup samples
outliers_library_size <- isOutlier(rse_gene_pups$sum, nmads = 3, type = "lower")
outliers_detected_num <- isOutlier(rse_gene_pups$detected, nmads = 3, type = "lower")
outliers_totalAssignedGene <- isOutlier(rse_gene_pups$totalAssignedGene, nmads = 3, type = "lower")
outliers_mito <- isOutlier(rse_gene_pups$mitoRate, nmads = 3, type = "higher")
outliers_rRNArate <- isOutlier(rse_gene_pups$rRNA_rate, nmads = 3, type = "higher")

not_outliers <- which(!(outliers_library_size | outliers_detected_num | outliers_totalAssignedGene | outliers_mito | outliers_rRNArate))
rse_gene_pups_qc <- rse_gene_pups[, not_outliers]


## ----message=FALSE, warning=FALSE------------------------------------
library("variancePartition")
library("pheatmap")

#######################   Variance Partition Analysis   #######################

## Fraction of variation attributable to each variable after correcting for all other variables


## 1. Canonical Correlation Analysis (CCA)

## Asses the correlation between each pair of sample variables

## Plot heatmap of correlations
plot_CCA <- function(age) {
    ## Data
    rse_gene <- eval(parse_expr(paste0("rse_gene_", age, "_qc")))


    ## Define variables to examine: remove those with single values
    ## For adults: all are females (so we drop 'Sex' variable)

    if (age == "adults") {
        formula <- ~ Group + Pregnancy + plate + flowcell + mitoRate + overallMapRate + totalAssignedGene + rRNA_rate + sum + detected + ERCCsumLogErr
    }
    ## For pups: none is pregnant (so 'Pregnancy' variable is not considered)
    else {
        formula <- ~ Group + Sex + plate + flowcell + mitoRate + overallMapRate + totalAssignedGene + rRNA_rate + sum + detected + ERCCsumLogErr
    }

    ## Measure correlations
    C <- canCorPairs(formula, colData(rse_gene))
    ## Heatmap
    pheatmap(
        C, ## data
        color = hcl.colors(50, "YlOrRd", rev = TRUE), ## color scale
        fontsize = 8, ## text size
        border_color = "black", ## border color for heatmap cells
        cellwidth = unit(0.4, "cm"), ## height of cells
        cellheight = unit(0.4, "cm") ## width of cells
    )

    return(C)
}


## ----message=FALSE, warning=FALSE------------------------------------
## Heatmap for adult samples
CCA_adults <- plot_CCA("adults")


## ----message=FALSE, warning=FALSE------------------------------------
## Heatmap for pup samples
CCA_pups <- plot_CCA("pups")


## ----message=FALSE, warning=FALSE------------------------------------
library("ggplot2")
library("rlang")
## 1.1  Barplots/Boxplots/Scatterplots for each pair of correlated variables

corr_plots <- function(age, sample_var1, sample_var2, sample_color) {
    ## Data
    rse_gene <- eval(parse_expr(paste("rse_gene", age, "qc", sep = "_")))
    CCA <- eval(parse_expr(paste0("CCA_", age)))

    ## Sample color by one variable
    colors <- list(
        "Group" = c("Control" = "brown2", "Experimental" = "deepskyblue3"),
        "Age" = c("Adult" = "slateblue3", "Pup" = "yellow3"),
        "Sex" = c("F" = "hotpink1", "M" = "dodgerblue"),
        "Pregnancy" = c("Yes" = "darkorchid3", "No" = "darkolivegreen4"),
        "plate" = c("Plate1" = "darkorange", "Plate2" = "lightskyblue", "Plate3" = "deeppink1"),
        "flowcell" = c(
            "HKCG7DSXX" = "chartreuse2", "HKCMHDSXX" = "magenta", "HKCNKDSXX" = "turquoise3",
            "HKCTMDSXX" = "tomato"
        )
    )

    data <- colData(rse_gene)

    ## Barplots for categorical variable vs categorical variable
    if (class(data[, sample_var1]) == "character" & class(data[, sample_var2]) == "character") {
        ## y-axis label
        if (sample_var2 == "Pregnancy") {
            y_label <- paste("Number of samples from each ", sample_var2, " group", sep = "")
        } else {
            y_label <- paste("Number of samples from each ", sample_var2, sep = "")
        }


        # Stacked barplot with counts for 2nd variable
        plot <- ggplot(data = as.data.frame(data), aes(
            x = !!rlang::sym(sample_var1),
            fill = !!rlang::sym(sample_var2)
        )) +
            geom_bar(position = "stack") +
            ## Colors by 2nd variable
            scale_fill_manual(values = colors[[sample_var2]]) +
            ## Show sample counts on stacked bars
            geom_text(aes(label = after_stat(count)),
                stat = "count",
                position = position_stack(vjust = 0.5), colour = "gray20", size = 3
            ) +
            theme_bw() +
            labs(
                subtitle = paste0("Corr: ", signif(CCA[sample_var1, sample_var2], digits = 3)),
                y = y_label
            ) +
            theme(
                axis.title = element_text(size = (7)),
                axis.text = element_text(size = (6)),
                plot.subtitle = element_text(size = 7, color = "gray40"),
                legend.text = element_text(size = 6),
                legend.title = element_text(size = 7)
            )
    }


    ## Boxplots for categorical variable vs continuous variable
    else if (class(data[, sample_var1]) == "character" & class(data[, sample_var2]) == "numeric") {
        plot <- ggplot(data = as.data.frame(data), mapping = aes(
            x = !!rlang::sym(sample_var1),
            y = !!rlang::sym(sample_var2),
            color = !!rlang::sym(sample_var1)
        )) +
            geom_boxplot(size = 0.25, width = 0.32, color = "black", outlier.color = "#FFFFFFFF") +
            geom_jitter(width = 0.15, alpha = 1, size = 1) +
            stat_smooth(geom = "line", alpha = 0.6, size = 0.4, span = 0.3, method = lm, aes(group = 1), color = "orangered3") +
            scale_color_manual(values = colors[[sample_var1]]) +
            theme_bw() +
            guides(color = "none") +
            labs(
                subtitle = paste0("Corr: ", signif(CCA[sample_var1, sample_var2], digits = 3)), y = gsub("_", " ", sample_var2),
                x = sample_var1
            ) +
            theme(
                axis.title = element_text(size = (7)),
                axis.text = element_text(size = (6)),
                plot.subtitle = element_text(size = 7, color = "gray40"),
                legend.text = element_text(size = 6),
                legend.title = element_text(size = 7)
            )
    }


    ## Scatterplots for continuous variable vs continuous variable
    else if (class(data[, sample_var1]) == "numeric" & class(data[, sample_var2]) == "numeric") {
        plot <- ggplot(as.data.frame(data), aes(
            x = !!rlang::sym(sample_var1),
            y = !!rlang::sym(sample_var2),
            color = !!rlang::sym(sample_color)
        )) +
            geom_point(size = 2) +
            stat_smooth(geom = "line", alpha = 0.4, size = 0.4, span = 0.25, method = lm, color = "orangered3") +
            ## Color by sample_color variale
            scale_color_manual(name = sample_color, values = colors[[sample_color]]) +
            theme_bw() +
            labs(subtitle = paste0("Corr: ", signif(CCA[sample_var1, sample_var2], digits = 3)), y = gsub("_", " ", sample_var2), x = gsub("_", " ", sample_var1)) +
            theme(
                axis.title = element_text(size = (7)),
                axis.text = element_text(size = (6)),
                plot.subtitle = element_text(size = 7, color = "gray40"),
                legend.text = element_text(size = 6),
                legend.title = element_text(size = 7)
            )
    }

    return(plot)
}


## ----message=FALSE, warning=FALSE------------------------------------
## Correlation plot for adults
p <- corr_plots("adults", "mitoRate", "totalAssignedGene", "Group")
p + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))


## ----message=FALSE, warning=FALSE------------------------------------
p <- corr_plots("adults", "flowcell", "plate", NULL)
p + theme(plot.margin = unit(c(1.5, 4.5, 1.5, 4.5), "cm"))


## ----message=FALSE, warning=FALSE------------------------------------
p <- corr_plots("adults", "plate", "overallMapRate", NULL)
p + theme(plot.margin = unit(c(2, 5.3, 2, 5.3), "cm"))


## ----message=FALSE, warning=FALSE------------------------------------
## Correlation plots
p <- corr_plots("adults", "sum", "detected", "Group")
p + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))

p <- corr_plots("pups", "sum", "detected", "Group")
p + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))


## ----message=FALSE, warning=FALSE------------------------------------
## ## Correlation plot for pups
p <- corr_plots("pups", "rRNA_rate", "overallMapRate", "Group")
p + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))


## ----message=FALSE, warning=FALSE------------------------------------
p <- corr_plots("pups", "plate", "overallMapRate", NULL)
p + theme(plot.margin = unit(c(2, 5.3, 2, 5.3), "cm"))


## ----message=FALSE, warning=FALSE------------------------------------
p <- corr_plots("pups", "flowcell", "overallMapRate", NULL)
p + theme(plot.margin = unit(c(2, 5.3, 2, 5.3), "cm"))


## ----message=FALSE, warning=FALSE------------------------------------
p1 <- corr_plots("adults", "Group", "plate", NULL)
p2 <- corr_plots("pups", "Group", "plate", NULL)
p3 <- corr_plots("adults", "Group", "flowcell", NULL)
p4 <- corr_plots("pups", "Group", "flowcell", NULL)


plots <- plot_grid(p1, p2, p3, p4, ncol = 2)
plots + theme(plot.margin = unit(c(1, 2.5, 1, 2.5), "cm"))


## ----message=FALSE, warning=FALSE------------------------------------
## 2. Fit model

## Fit a linear mixed model (LMM) that takes continuous variables as fixed effects and categorical variables as random effects

varPartAnalysis <- function(age, formula) {
    RSE <- eval(parse_expr(paste("rse_gene", age, "qc", sep = "_")))

    ## Ignore genes with variance 0
    genes_var_zero <- which(apply(assays(RSE)$logcounts, 1, var) == 0)
    if (length(genes_var_zero) > 0) {
        RSE <- RSE[-genes_var_zero, ]
    }

    ## Loop over each gene to fit model and extract variance explained by each variable
    varPart <- fitExtractVarPartModel(assays(RSE)$logcounts, formula, colData(RSE))

    # Sort variables by median fraction of variance explained
    vp <- sortCols(varPart)
    p <- plotVarPart(vp)

    return(list(p, vp))
}


## ----message=FALSE, warning=FALSE, eval=FALSE------------------------
## ## Violin plots
## 
## #####  Model with all variables  #####
## 
## ## Adults
## ## Define variables; random effects indicated with (1| )
## formula <- ~ (1 | Group) + (1 | Pregnancy) + (1 | plate) + (1 | flowcell) + mitoRate + overallMapRate +
##     totalAssignedGene + rRNA_rate + sum + detected + ERCCsumLogErr
## plot <- varPartAnalysis("adults", formula)[[1]]
## plot + theme(
##     plot.margin = unit(c(1, 1, 1, 1), "cm"),
##     axis.text.x = element_text(size = (7)),
##     axis.text.y = element_text(size = (7.5))
## )


## ----message=FALSE, warning=FALSE, eval=FALSE------------------------
## #####  Model without correlated variables  #####
## 
## ## Adult plots without mitoRate, plate and sum
## formula <- ~ (1 | Group) + (1 | Pregnancy) + (1 | flowcell) + overallMapRate + totalAssignedGene + rRNA_rate + detected + ERCCsumLogErr
## varPart <- varPartAnalysis("adults", formula)
## varPart_data_adults <- varPart[[2]]
## plot <- varPart[[1]]
## plot + theme(
##     plot.margin = unit(c(1, 1, 1, 1), "cm"),
##     axis.text.x = element_text(size = (7)),
##     axis.text.y = element_text(size = (7.5))
## )


## ----message=FALSE, warning=FALSE, eval=FALSE------------------------
## #####  Model with all variables  #####
## 
## ## Pups
## formula <- ~ (1 | Group) + (1 | Sex) + (1 | plate) + (1 | flowcell) + mitoRate + overallMapRate +
##     totalAssignedGene + rRNA_rate + sum + detected + ERCCsumLogErr
## plot <- varPartAnalysis("pups", formula)[[1]]
## plot + theme(
##     plot.margin = unit(c(1, 1, 1, 1), "cm"),
##     axis.text.x = element_text(size = (7)),
##     axis.text.y = element_text(size = (7.5))
## )


## ----message=FALSE, warning=FALSE, eval=FALSE------------------------
## #####  Model without correlated variables  #####
## 
## ## Pup plots without sum, rRNA_rate and plate
## formula <- ~ (1 | Group) + (1 | Sex) + (1 | flowcell) + mitoRate + overallMapRate + totalAssignedGene + detected + ERCCsumLogErr
## varPart <- varPartAnalysis("pups", formula)
## varPart_data_pups <- varPart[[2]]
## plot <- varPart[[1]]
## plot + theme(
##     plot.margin = unit(c(1, 1, 1, 1), "cm"),
##     axis.text.x = element_text(size = (7)),
##     axis.text.y = element_text(size = (7.5))
## )


## ----message=FALSE, warning=FALSE------------------------------------
## Plot of gene expression lognorm counts vs. sample variable
plot_gene_expr <- function(age, sample_var, gene_id) {
    rse_gene <- eval(parse_expr(paste("rse_gene", age, "qc", sep = "_")))
    varPart_data <- eval(parse_expr(paste0("varPart_data_", age)))

    colors <- list(
        "Group" = c("Control" = "brown2", "Experimental" = "deepskyblue3"),
        "Age" = c("Adult" = "slateblue3", "Pup" = "yellow3"),
        "Sex" = c("F" = "hotpink1", "M" = "dodgerblue"),
        "Pregnancy" = c("Yes" = "darkorchid3", "No" = "darkolivegreen4"),
        "plate" = c("Plate1" = "darkorange", "Plate2" = "lightskyblue", "Plate3" = "deeppink1"),
        "flowcell" = c(
            "HKCG7DSXX" = "chartreuse2", "HKCMHDSXX" = "magenta", "HKCNKDSXX" = "turquoise3",
            "HKCTMDSXX" = "tomato"
        )
    )

    ## Lognorm counts of the gene across samples
    data <- colData(rse_gene)
    data$gene_expr <- assays(rse_gene)$logcounts[gene_id, ]

    ## Percentage of variance explained by the variable
    percentage <- 100 * signif(varPart_data[gene_id, sample_var], digits = 3)

    ## Boxplots for discrete variables
    if (class(data[, sample_var]) == "character") {
        plot <- ggplot(data = as.data.frame(data), mapping = aes(
            x = !!rlang::sym(sample_var),
            y = gene_expr, color = !!rlang::sym(sample_var)
        )) +
            geom_boxplot(size = 0.25, width = 0.32, color = "black", outlier.color = "#FFFFFFFF") +
            geom_jitter(width = 0.15, alpha = 1, size = 1) +
            stat_smooth(geom = "line", alpha = 0.6, size = 0.4, span = 0.3, method = lm, aes(group = 1), color = "orangered3") +
            scale_color_manual(values = colors[[sample_var]]) +
            theme_bw() +
            guides(color = "none") +
            labs(
                title = gene_id,
                subtitle = paste0("Variance explained: ", percentage, "%"),
                y = "lognorm counts", x = sample_var
            ) +
            theme(
                axis.title = element_text(size = (7)),
                axis.text = element_text(size = (6)),
                plot.title = element_text(hjust = 0.5, size = 7.5, face = "bold"),
                plot.subtitle = element_text(size = 7, color = "gray40"),
                legend.text = element_text(size = 6),
                legend.title = element_text(size = 7)
            )
    }

    ## Scatterplots for continuous variables
    else {
        colors <- c(
            "mitoRate" = "khaki3", "overallMapRate" = "turquoise", "totalAssignedGene" = "plum2", "rRNA_rate" = "orange3",
            "sum" = "palegreen3", "detected" = "skyblue2", "ERCCsumLogErr" = "slateblue1"
        )

        plot <- ggplot(as.data.frame(data), aes(x = eval(parse_expr(sample_var)), y = gene_expr)) +
            geom_point(color = colors[[sample_var]], size = 2) +
            stat_smooth(geom = "line", alpha = 0.4, size = 0.4, span = 0.25, method = lm, color = "orangered3") +
            theme_bw() +
            guides(color = "none") +
            labs(
                title = gene_id,
                subtitle = paste0("Variance explained: ", percentage, "%"),
                y = "lognorm counts", x = gsub("_", " ", sample_var)
            ) +
            theme(
                plot.margin = unit(c(0.4, 0.1, 0.4, 0.1), "cm"),
                axis.title = element_text(size = (7)),
                axis.text = element_text(size = (6)),
                plot.title = element_text(hjust = 0.5, size = 7.5, face = "bold"),
                plot.subtitle = element_text(size = 7, color = "gray40"),
                legend.text = element_text(size = 6),
                legend.title = element_text(size = 7)
            )
    }

    return(plot)
}


## ----message=FALSE, warning=FALSE------------------------------------
## Function to plot gene expression vs sample variable data for top 3 most affected genes

plot_gene_expr_sample <- function(age, sample_var) {
    rse_gene <- eval(parse_expr(paste("rse_gene", age, "qc", sep = "_")))
    varPart_data <- eval(parse_expr(paste0("varPart_data_", age)))

    ## Top 3 genes most affected by sample variable
    affected_genes <- rownames(varPart_data[order(varPart_data[, "Group"], decreasing = TRUE), ][1:3, ])

    ## Plots
    plots <- list()
    for (i in 1:length(affected_genes)) {
        plots[[i]] <- plot_gene_expr(age, sample_var, affected_genes[i])
    }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 3)
}


## ----message=FALSE, warning=FALSE, eval=FALSE------------------------
## ## Adults
## 
## ## Plots for top affected genes by 'totalAssignedGene'
## plots <- plot_gene_expr_sample("adults", "totalAssignedGene")
## plots + theme(plot.margin = unit(c(3, 1, 2, 3), "cm"))
## 
## ## Plots for top affected genes by 'overallMapRate'
## plots <- plot_gene_expr_sample("adults", "overallMapRate")
## plots + theme(plot.margin = unit(c(3, 1, 2, 3), "cm"))
## 
## ## Plots for top affected genes by 'Group'
## plots <- plot_gene_expr_sample("adults", "Group")
## plots + theme(plot.margin = unit(c(3, 1, 2, 3), "cm"))


## ----message=FALSE, warning=FALSE, eval=FALSE------------------------
## ## Pups
## 
## ## Plots for top affected genes by 'overallMapRate'
## plots <- plot_gene_expr_sample("pups", "overallMapRate")
## plots + theme(plot.margin = unit(c(3, 1, 2, 3), "cm"))
## 
## ## Plots for top affected genes by 'totalAssignedGene'
## plots <- plot_gene_expr_sample("pups", "totalAssignedGene")
## plots + theme(plot.margin = unit(c(3, 1, 2, 3), "cm"))
## 
## ## Plots for top affected genes by 'Group'
## plots <- plot_gene_expr_sample("pups", "Group")
## plots + theme(plot.margin = unit(c(3, 1, 2, 3), "cm"))


## ----exercise1_varPart, message=FALSE, warning=FALSE,  eval=FALSE, echo=FALSE----
## ## Solution
## 
## ## Gene ID
## gene_id <- "ENSMUSG00000042348.10"
## ## % of variance explained by Group
## percentage <- 100 * signif(varPart_data_pups[gene_id, "Group"], digits = 3)
## ## Sample colors
## colors <- c("Control" = "brown2", "Experimental" = "deepskyblue3")
## ## Gene expression logcounts
## rse_gene_pups_qc$gene_expr <- assays(rse_gene_pups_qc)$logcounts[gene_id, ]
## 
## ## Plot
## plot <- ggplot(
##     data = as.data.frame(colData(rse_gene_pups_qc)),
##     mapping = aes(x = Group, y = gene_expr, color = Group)
## ) +
##     geom_boxplot(size = 0.25, width = 0.32, color = "black", outlier.color = "#FFFFFFFF") +
##     geom_jitter(width = 0.15, alpha = 1, size = 1) +
##     scale_color_manual(values = colors) +
##     theme_bw() +
##     guides(color = "none") +
##     labs(
##         title = gene_id,
##         subtitle = paste0("Variance explained: ", percentage, "%"),
##         y = "lognorm counts"
##     ) +
##     theme(
##         plot.margin = unit(c(2, 6, 2, 6), "cm"),
##         axis.title = element_text(size = (7)),
##         axis.text = element_text(size = (6)),
##         plot.title = element_text(hjust = 0.5, size = 7.5, face = "bold"),
##         plot.subtitle = element_text(size = 7, color = "gray40"),
##         legend.text = element_text(size = 6),
##         legend.title = element_text(size = 7)
##     )
## 
## plot


## ----exercise2_varPart, message=FALSE, warning=FALSE,  eval=FALSE, echo=FALSE----
## ## Solution
## 
## ## Gene ID
## gene_id <- "ENSMUSG00000064372.1"
## ## % of variance explained by Group
## percentage <- 100 * signif(varPart_data_pups[gene_id, "Group"], digits = 3)
## ## Sample colors
## colors <- c("Control" = "brown2", "Experimental" = "deepskyblue3")
## ## Gene expression logcounts
## rse_gene_pups_qc$gene_expr <- assays(rse_gene_pups_qc)$logcounts[gene_id, ]
## 
## ## Plot
## plot <- ggplot(
##     data = as.data.frame(colData(rse_gene_pups_qc)),
##     mapping = aes(x = Group, y = gene_expr, color = Group)
## ) +
##     geom_boxplot(size = 0.25, width = 0.32, color = "black", outlier.color = "#FFFFFFFF") +
##     geom_jitter(width = 0.15, alpha = 1, size = 1) +
##     scale_color_manual(values = colors) +
##     theme_bw() +
##     guides(color = "none") +
##     labs(
##         title = gene_id,
##         subtitle = paste0("Variance explained: ", percentage, "%"),
##         y = "lognorm counts"
##     ) +
##     theme(
##         plot.margin = unit(c(2, 6, 2, 6), "cm"),
##         axis.title = element_text(size = (7)),
##         axis.text = element_text(size = (6)),
##         plot.title = element_text(hjust = 0.5, size = 7.5, face = "bold"),
##         plot.subtitle = element_text(size = 7, color = "gray40"),
##         legend.text = element_text(size = 6),
##         legend.title = element_text(size = 7)
##     )
## 
## plot

