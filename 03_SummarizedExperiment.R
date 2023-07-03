## ---- include = FALSE------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)


## ----vignetteSetup_SEreview, echo=FALSE, message=FALSE, warning = FALSE----
## For links
library(BiocStyle)

## Bib setup
library(RefManageR)

## Write bibliography information
bib <- c(
    smokingMouse = citation("smokingMouse")[1],
    SummarizedExperiment = citation("SummarizedExperiment")[1]
)

options(max.print = 50)


## ---- echo=FALSE-----------------------------------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(data(airway, package = "airway"))


## --------------------------------------------------------------------
library(SummarizedExperiment)
library(airway)

data(airway, package = "airway")
se <- airway


## --------------------------------------------------------------------
## For a) you could only print the summary of the object but since the idea is to understand
## how to explore the object find other function that gives you the answer.
se

## Same thing for b, you could just print the colData and count the samples, but this is not
## efficient when our data consists in hundreds of samples. Find the answer using other tools.

colData(se)


## --------------------------------------------------------------------
## In our object, if you look at the part that says assays, we can see that at the moment
## we only have one with the name "counts"

se

## To see the data that's stored in that assay you can do either one of the next commands
assay(se)
assays(se)$counts

## Note that assay() does not support $ operator
# assay(se)$counts

## We would have to do:
assay(se, 1)
assay(se, "counts")

## If you use assays() without specifying the element you want to see it shows you the length
## of the list and the name of each element
assays(se)

## To obtain a list of names as a vector you can do:
assayNames(se)

## Which can also be use to change the name of the assays
assayNames(se)[1] <- "foo"
assayNames(se)
assayNames(se)[1] <- "counts"


## --------------------------------------------------------------------
## To calculate the library size use

apply(assay(se), 2, sum)

