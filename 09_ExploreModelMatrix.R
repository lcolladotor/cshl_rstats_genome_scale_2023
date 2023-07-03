## ----model.matrix----------------------------------------------------
## ?model.matrix
mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat
colnames(mat)


## ----lm_example------------------------------------------------------
summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))


## ----EMM_example1----------------------------------------------------
## Example data
(sampleData <- data.frame(
    genotype = rep(c("A", "B"), each = 4),
    treatment = rep(c("ctrl", "trt"), 4)
))

## Let's make the visual aids provided by ExploreModelMatrix
vd <- ExploreModelMatrix::VisualizeDesign(
    sampleData = sampleData,
    designFormula = ~ genotype + treatment,
    textSizeFitted = 4
)

## Now lets plot these images
cowplot::plot_grid(plotlist = vd$plotlist)


## ----EMM_example1_interactive, eval = FALSE--------------------------
## ## We are using shiny again here
## app <- ExploreModelMatrix(
##     sampleData = sampleData,
##     designFormula = ~ genotype + treatment
## )
## if (interactive()) shiny::runApp(app)

