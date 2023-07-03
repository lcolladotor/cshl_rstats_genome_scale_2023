## ----download_data_biocfilecache_repeat--------------------------
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


## ----Data preparation, message=FALSE, warning=FALSE--------------
library("ggplot2")

## Histogram and density plot of read counts before and after normalization

## Raw counts
counts_data <- data.frame(counts = as.vector(assays(rse_gene_nic)$counts))
plot <- ggplot(counts_data, aes(x = counts)) +
    geom_histogram(colour = "black", fill = "lightgray") +
    labs(x = "read counts", y = "Frecuency") +
    theme_classic()
plot + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))

## Normalized counts
logcounts_data <- data.frame(logcounts = as.vector(assays(rse_gene_nic)$logcounts))
plot <- ggplot(logcounts_data, aes(x = logcounts)) +
    geom_histogram(aes(y = ..density..), colour = "darkgray", fill = "lightgray") +
    theme_classic() +
    geom_density(fill = "#69b3a2", alpha = 0.3) +
    labs(x = "log(CPM+0.5)", y = "Frecuency")
plot + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))


## ----  message=FALSE, warning=FALSE------------------------------
## Retain genes that passed filtering step
rse_gene_filt <- rse_gene_nic[rowData(rse_gene_nic)$retained_after_feature_filtering == TRUE, ]

## Normalized counts and filtered genes
filt_logcounts_data <- data.frame(logcounts = as.vector(assays(rse_gene_filt)$logcounts))

## Plot
plot <- ggplot(filt_logcounts_data, aes(x = logcounts)) +
    geom_histogram(aes(y = ..density..), colour = "darkgray", fill = "lightgray") +
    theme_classic() +
    geom_density(fill = "#69b3a2", alpha = 0.3) +
    labs(x = "log(CPM+0.5)", y = "Frecuency")
plot + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))


## ----QC_boxplots,  message=FALSE, warning=FALSE------------------
library("Hmisc")
library("stringr")
library("cowplot")

## Define QC metrics of interest
qc_metrics <- c("mitoRate", "overallMapRate", "totalAssignedGene", "rRNA_rate", "sum", "detected", "ERCCsumLogErr")

## Define sample variables of interest
sample_variables <- c("Group", "Age", "Sex", "Pregnancy", "plate", "flowcell")


## Function to create boxplots of QC metrics for groups of samples

QC_boxplots <- function(qc_metric, sample_var) {
    ## Define sample colors depending on the sample variable
    if (sample_var == "Group") {
        colors <- c("Control" = "brown2", "Experimental" = "deepskyblue3")
    } else if (sample_var == "Age") {
        colors <- c("Adult" = "slateblue3", "Pup" = "yellow3")
    } else if (sample_var == "Sex") {
        colors <- c("F" = "hotpink1", "M" = "dodgerblue")
    } else if (sample_var == "Pregnancy") {
        colors <- c("Yes" = "darkorchid3", "No" = "darkolivegreen4")
    } else if (sample_var == "plate") {
        colors <- c("Plate1" = "darkorange", "Plate2" = "lightskyblue", "Plate3" = "deeppink1")
    } else if (sample_var == "flowcell") {
        colors <- c(
            "HKCG7DSXX" = "chartreuse2", "HKCMHDSXX" = "magenta", "HKCNKDSXX" = "turquoise3",
            "HKCTMDSXX" = "tomato"
        )
    }

    ## Axis labels
    x_label <- capitalize(sample_var)
    y_label <- str_replace_all(qc_metric, c("_" = ""))

    ## x-axis text angle and position
    if (sample_var == "flowcell") {
        x_axis_angle <- 18
        x_axis_hjust <- 0.5
        x_axis_vjust <- 0.7
        x_axis_size <- 4
    } else {
        x_axis_angle <- 0
        x_axis_hjust <- 0.5
        x_axis_vjust <- 0.5
        x_axis_size <- 6
    }

    ## Extract sample data in colData(rse_gene_filt)
    data <- data.frame(colData(rse_gene_filt))

    ## Sample variable separating samples in x-axis and QC metric in y-axis
    ## (Coloring by sample variable)
    plot <- ggplot(data = data, mapping = aes(x = !!rlang::sym(sample_var), y = !!rlang::sym(qc_metric), color = !!rlang::sym(sample_var))) +
        ## Add violin plots
        geom_violin(alpha = 0, size = 0.4, color = "black", width = 0.7) +
        ## Spread dots
        geom_jitter(width = 0.08, alpha = 0.7, size = 1.3) +
        ## Add boxplots
        geom_boxplot(alpha = 0, size = 0.4, width = 0.1, color = "black") +
        ## Set colors
        scale_color_manual(values = colors) +
        ## Define axis labels
        labs(y = y_label, x = x_label) +
        ## Get rid of the background
        theme_bw() +
        ## Hide legend and define plot margins and size of axis title and text
        theme(
            legend.position = "none",
            plot.margin = unit(c(0.5, 0.4, 0.5, 0.4), "cm"),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = x_axis_size),
            axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, vjust = x_axis_vjust)
        )


    return(plot)
}



## Plots of all QC metrics for each sample variable
multiple_QC_boxplots <- function(sample_var) {
    i <- 1
    plots <- list()
    for (qc_metric in qc_metrics) {
        ## Call function to create each individual plot
        plots[[i]] <- QC_boxplots(qc_metric, sample_var)
        i <- i + 1
    }
    ## Arrange multiple plots into a grid
    print(plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], nrow = 2))
}


## ----message=FALSE, warning=FALSE--------------------------------
multiple_QC_boxplots("Age")


## ----message=FALSE, warning=FALSE--------------------------------
multiple_QC_boxplots("Sex")


## ----message=FALSE, warning=FALSE--------------------------------
multiple_QC_boxplots("Group")


## ----message=FALSE, warning=FALSE--------------------------------
multiple_QC_boxplots("Pregnancy")


## ----message=FALSE, warning=FALSE--------------------------------
multiple_QC_boxplots("plate")


## ----message=FALSE, warning=FALSE--------------------------------
multiple_QC_boxplots("flowcell")


## ----"QC scatterplots", message=FALSE, warning=FALSE-------------
## Scatterplots for a pair of QC metrics

QC_scatterplots <- function(sample_var, qc_metric1, qc_metric2) {
    ## Define sample colors
    if (sample_var == "Group") {
        colors <- c("Control" = "brown2", "Experimental" = "deepskyblue3")
    } else if (sample_var == "Age") {
        colors <- c("Adult" = "slateblue3", "Pup" = "yellow3")
    } else if (sample_var == "Sex") {
        colors <- c("F" = "hotpink1", "M" = "dodgerblue")
    } else if (sample_var == "Pregnancy") {
        colors <- c("Yes" = "darkorchid3", "No" = "darkolivegreen4")
    } else if (sample_var == "plate") {
        colors <- c("Plate1" = "darkorange", "Plate2" = "lightskyblue", "Plate3" = "deeppink1")
    } else if (sample_var == "flowcell") {
        colors <- c(
            "HKCG7DSXX" = "chartreuse2", "HKCMHDSXX" = "magenta", "HKCNKDSXX" = "turquoise3",
            "HKCTMDSXX" = "tomato"
        )
    }

    data <- colData(rse_gene_filt)

    ## Scatterplots for continuous variable vs continuous variable
    ## First QC metric in x-axis and second QC metric in y-axis
    plot <- ggplot(as.data.frame(data), aes(
        x = !!rlang::sym(qc_metric1),
        y = !!rlang::sym(qc_metric2),
        ## Color samples by a variable
        color = !!rlang::sym(sample_var)
    )) +
        ## Add scatterplot
        geom_point(size = 1) +
        ## Add regression line
        stat_smooth(geom = "line", alpha = 0.4, size = 0.4, span = 0.25, method = lm, color = "orangered3") +
        ## Colors
        scale_color_manual(name = sample_var, values = colors) +
        theme_bw() +
        ## Add Pearson correlation coefficient between the metrics as subtitle
        labs(
            subtitle = paste0("Corr: ", signif(cor(data[, qc_metric1], data[, qc_metric2], method = "pearson"), digits = 3)),
            ## Add axis labels
            y = gsub("_", " ", qc_metric2),
            x = gsub("_", " ", qc_metric1)
        ) +
        ## Plot margins and text size
        theme(
            plot.margin = unit(c(0.1, 1.2, 0.1, 1.2), "cm"),
            axis.title = element_text(size = (7)),
            axis.text = element_text(size = (6)),
            plot.subtitle = element_text(size = 7, color = "gray40"),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 7)
        )
    return(plot)
}



## QC scatterplots coloring by all sample variables
multiple_QC_scatterplots <- function(qc_metric1, qc_metric2) {
    sample_variables <- c("Age", "Sex", "plate", "Pregnancy", "Group", "flowcell")

    i <- 1
    plots <- list()
    for (sample_var in sample_variables) {
        plots[[i]] <- QC_scatterplots(sample_var, qc_metric1, qc_metric2)
        i <- i + 1
    }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], nrow = 3, rel_widths = c(1, 1))
}


## ----message=FALSE, warning=FALSE--------------------------------
multiple_QC_scatterplots("mitoRate", "rRNA_rate")


## ----message=FALSE, warning=FALSE--------------------------------
multiple_QC_scatterplots("mitoRate", "totalAssignedGene")


## ----message=FALSE, warning=FALSE--------------------------------
multiple_QC_scatterplots("sum", "detected")


## ----message=FALSE, warning=FALSE--------------------------------
multiple_QC_scatterplots("sum", "totalAssignedGene")


## ----message=FALSE, warning=FALSE--------------------------------
multiple_QC_scatterplots("detected", "totalAssignedGene")


## ----exercise1_EDA, message=FALSE, warning=FALSE,  eval=FALSE, echo=FALSE----
## ## Solution
## multiple_QC_scatterplots("subsets_Mito_sum", "mitoMapped")
## ## Because in mitoMapped you take reads that mapped to the whole mt chr, in subsets_Mito_sum only reads that were aligned to mt genes. But there's almost a perfect correlation between these two metrics.


## ----"QC sample filtering", message=FALSE, warning=FALSE---------
library("scater")
library("rlang")
library("ggrepel")

## Separate data by Age
rse_gene_pups <- rse_gene_filt[, which(rse_gene_filt$Age == "Pup")]
rse_gene_adults <- rse_gene_filt[, which(rse_gene_filt$Age == "Adult")]


## Find outlier samples based on their QC metrics (samples that are 3 median-absolute-deviations away from the median)


## Filter all samples together

## Drop samples with lower library sizes (sum), detected number of genes and totalAssignedGene
outliers_library_size <- isOutlier(rse_gene_filt$sum, nmads = 3, type = "lower")
outliers_detected_num <- isOutlier(rse_gene_filt$detected, nmads = 3, type = "lower")
outliers_totalAssignedGene <- isOutlier(rse_gene_filt$totalAssignedGene, nmads = 3, type = "lower")
## Drop samples with higher mitoRates and rRNA rates
outliers_mito <- isOutlier(rse_gene_filt$mitoRate, nmads = 3, type = "higher")
outliers_rRNArate <- isOutlier(rse_gene_filt$rRNA_rate, nmads = 3, type = "higher")

## Keep not outlier samples
not_outliers <- which(!(outliers_library_size | outliers_detected_num | outliers_totalAssignedGene | outliers_mito | outliers_rRNArate))
rse_gene_qc <- rse_gene_filt[, not_outliers]

## Number of samples retained
dim(rse_gene_qc)[2]

## Add new variables to rse_gene_filt with info of samples retained/dropped
rse_gene_filt$Retention_after_QC_filtering <- as.vector(sapply(rse_gene_filt$SAMPLE_ID, function(x) {
    if (x %in% rse_gene_qc$SAMPLE_ID) {
        "Retained"
    } else {
        "Dropped"
    }
}))



## Filter adult samples

outliers_library_size <- isOutlier(rse_gene_adults$sum, nmads = 3, type = "lower")
outliers_detected_num <- isOutlier(rse_gene_adults$detected, nmads = 3, type = "lower")
outliers_totalAssignedGene <- isOutlier(rse_gene_adults$totalAssignedGene, nmads = 3, type = "lower")
outliers_mito <- isOutlier(rse_gene_adults$mitoRate, nmads = 3, type = "higher")
outliers_rRNArate <- isOutlier(rse_gene_adults$rRNA_rate, nmads = 3, type = "higher")

not_outliers <- which(!(outliers_library_size | outliers_detected_num | outliers_totalAssignedGene | outliers_mito | outliers_rRNArate))
rse_gene_adults_qc <- rse_gene_adults[, not_outliers]

## Number of samples retained
dim(rse_gene_adults_qc)[2]

rse_gene_adults$Retention_after_QC_filtering <- as.vector(sapply(rse_gene_adults$SAMPLE_ID, function(x) {
    if (x %in% rse_gene_adults_qc$SAMPLE_ID) {
        "Retained"
    } else {
        "Dropped"
    }
}))



## Filter pup samples

outliers_library_size <- isOutlier(rse_gene_pups$sum, nmads = 3, type = "lower")
outliers_detected_num <- isOutlier(rse_gene_pups$detected, nmads = 3, type = "lower")
outliers_totalAssignedGene <- isOutlier(rse_gene_pups$totalAssignedGene, nmads = 3, type = "lower")
outliers_mito <- isOutlier(rse_gene_pups$mitoRate, nmads = 3, type = "higher")
outliers_rRNArate <- isOutlier(rse_gene_pups$rRNA_rate, nmads = 3, type = "higher")

not_outliers <- which(!(outliers_library_size | outliers_detected_num | outliers_totalAssignedGene | outliers_mito | outliers_rRNArate))
rse_gene_pups_qc <- rse_gene_pups[, not_outliers]

## Number of samples retained
dim(rse_gene_pups_qc)[2]

rse_gene_pups$Retention_after_QC_filtering <- as.vector(sapply(rse_gene_pups$SAMPLE_ID, function(x) {
    if (x %in% rse_gene_pups_qc$SAMPLE_ID) {
        "Retained"
    } else {
        "Dropped"
    }
}))


## ----message=FALSE, warning=FALSE--------------------------------
## Boxplots of QC metrics after sample filtering

## Boxplots
boxplots_after_QC_filtering <- function(rse_gene, qc_metric, sample_var) {
    ## Color samples
    colors <- c("Retained" = "deepskyblue", "Dropped" = "brown2")

    ## Sample shape by sample variables
    if (sample_var == "Group") {
        shapes <- c("Control" = 0, "Experimental" = 15)
    } else if (sample_var == "Age") {
        shapes <- c("Adult" = 16, "Pup" = 1)
    } else if (sample_var == "Sex") {
        shapes <- c("F" = 11, "M" = 19)
    } else if (sample_var == "Pregnancy") {
        shapes <- c("Yes" = 10, "No" = 1)
    } else if (sample_var == "plate") {
        shapes <- c("Plate1" = 12, "Plate2" = 5, "Plate3" = 4)
    } else if (sample_var == "flowcell") {
        shapes <- c(
            "HKCG7DSXX" = 3, "HKCMHDSXX" = 8, "HKCNKDSXX" = 14,
            "HKCTMDSXX" = 17
        )
    }


    y_label <- str_replace_all(qc_metric, c("_" = " "))

    data <- data.frame(colData(rse_gene))

    ## Median of the QC var values
    median <- median(eval(parse_expr(paste("rse_gene$", qc_metric, sep = ""))))
    ## Median-absolute-deviation of the QC var values
    mad <- mad(eval(parse_expr(paste("rse_gene$", qc_metric, sep = ""))))

    plot <- ggplot(data = data, mapping = aes(
        x = "", y = !!rlang::sym(qc_metric),
        color = !!rlang::sym("Retention_after_QC_filtering")
    )) +
        geom_jitter(alpha = 1, size = 2, aes(shape = eval(parse_expr((sample_var))))) +
        geom_boxplot(alpha = 0, size = 0.15, color = "black") +
        scale_color_manual(values = colors) +
        scale_shape_manual(values = shapes) +
        labs(x = "", y = y_label, color = "Retention after QC filtering", shape = sample_var) +
        theme_classic() +
        ## Median line
        geom_hline(yintercept = median, size = 0.5) +
        ## Line of median + 3 MADs
        geom_hline(yintercept = median + (3 * mad), size = 0.5, linetype = 2) +
        ## Line of median - 3 MADs
        geom_hline(yintercept = median - (3 * mad), size = 0.5, linetype = 2) +
        theme(
            axis.title = element_text(size = (9)),
            axis.text = element_text(size = (8)),
            legend.position = "right",
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 9)
        )

    return(plot)
}


## ----message=FALSE, warning=FALSE--------------------------------
## Plots

## All samples together
p <- boxplots_after_QC_filtering(rse_gene_filt, "mitoRate", "Age")
p + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))


## ----message=FALSE, warning=FALSE--------------------------------
## Adult samples
p <- boxplots_after_QC_filtering(rse_gene_adults, "mitoRate", "Group")
p + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))


## ----message=FALSE, warning=FALSE--------------------------------
## Pup samples
p <- boxplots_after_QC_filtering(rse_gene_pups, "rRNA_rate", "Group")
p + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))


## ----exercise2_EDA, message=FALSE, warning=FALSE, eval=FALSE, echo=FALSE----
## ## Solution
## p <- boxplots_after_QC_filtering(rse_gene_adults, "mitoRate", "Pregnancy")
## p + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))
## 
## p <- boxplots_after_QC_filtering(rse_gene_adults, "mitoRate", "plate")
## p + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))
## 
## p <- boxplots_after_QC_filtering(rse_gene_adults, "mitoRate", "flowcell")
## p + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))
## 
## p <- boxplots_after_QC_filtering(rse_gene_adults, "sum", "Group")
## p + theme(plot.margin = unit(c(2, 4, 2, 4), "cm"))

