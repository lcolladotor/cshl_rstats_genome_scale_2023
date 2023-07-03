styler::style_dir(filetype = "Rmd", transformers = biocthis::bioc_style())

rmds <- dir(pattern = "\\.Rmd$")
sapply(rmds, knitr::knit, tangle = TRUE)
