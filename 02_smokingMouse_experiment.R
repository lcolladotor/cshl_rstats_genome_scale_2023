## ---- include = FALSE------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)


## ----vignetteSetup, echo=FALSE, message=FALSE, warning = FALSE-------
## For links
library(BiocStyle)

## Bib setup
library(RefManageR)

## Write bibliography information
bib <- c(
    smokingMouse = citation("smokingMouse")[1],
    SummarizedExperiment = citation("SummarizedExperiment")[1]
)

