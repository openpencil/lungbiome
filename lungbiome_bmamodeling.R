##*****************************************************************************************************
##  For complete background and details, please refer to:
##  Shankar, J. et al. Microbiome and cytokine signatures of bacterial pneumonia
##  and tracheobronchitis.(Manuscript under review) (2015).
##
##  Modeling was performed under the following environment:
##
## > sessionInfo()
## R version 3.2.0 (2015-04-16)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.10.3 (Yosemite)
##
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
##
## attached base packages:
## [1] parallel  grid      stats     graphics  grDevices utils     datasets  methods   base
##
## other attached packages:
## [1] BoomSpikeSlab_0.5.2 Boom_0.2            MASS_7.3-40         doMC_1.3.3          iterators_1.0.7
## [6] foreach_1.4.2       data.table_1.9.4    RColorBrewer_1.1-2  plyr_1.8.2          scales_0.2.4
## [11] reshape_0.8.5       ggplot2_1.0.1       rlecuyer_0.3-3
##
## loaded via a namespace (and not attached):
## [1] Rcpp_0.11.6      magrittr_1.5     munsell_0.4.2    colorspace_1.2-6 stringr_1.0.0    tools_3.2.0
## [7] gtable_0.1.2     digest_0.6.8     reshape2_1.4.1   codetools_0.2-11 stringi_0.5-5    chron_2.3-45
## [13] proto_0.3-10
##
## This script provides the code for estimating the BMA models.
## ******************************************************************************************************

##### 1. Load required libraries #####
## clone the GitHub directory and set it as the working directory
## setwd('/path_to_GitHub_clone')
## Create a directory 'output' for saving output from the analysis
dir.create("./output", showWarnings = T, recursive = T)

##### 2. Load algorithms and project data #####
source(file = "./lungbiome_algorithms.R")
source(file = "./lungbiome_loadprojectdata.R")

##### 3. Run BMA MLR classification model for pneumonia, tracheobronchitis and colonization #####
runmlrbmamodel <- function(responsevar, bmaiterations, regdata) {
    yvar <- regdata[, responsevar]
    xvar <- regdata[, setdiff(colnames(regdata), responsevar)]
    ## eliminate variables that are zero throughout or have a constant value in all columns
    zeroes <- which(colSums(xvar) == 0)
    if (length(zeroes) > 0) {
        xvar <- xvar[, -zeroes]
    }
    ## setting ems to the lowest that works. 8 for 3 levels. (addtomodelsize= 5 + 3 levels of outcome)
    bmaresults <- runmlrbma(x = xvar, y = yvar, addtomodelsize = 5, iter = bmaiterations,
        seed = 101)
    ## Save the model results
    saveRDS(object = bmaresults, file = "./output/lungbiome_MLR_diagnoses.RDS")
    return(outlist)
}
mlr_results <- runmlrbmamodel(responsevar = "diagnosis_simple_code",
                              bmaiterations = 80000,
                              regdata = regmatrix_diagnoses)


##### 5. Run BMA linear regression model for checking association of microbiome diversity with therapy #####
runlinearbma <- function(regmatrix, bmaiterations) {
    ## Eliminate variables that are zero throughout or have a constant value in all columns.
    zeroes <- which(sapply(colnames(regmatrix), function(cname) {
        out <- max(regmatrix[, cname]) - min(regmatrix[, cname]) == 0
        return(out)
    }))
    if (length(zeroes) > 0) {
        regmatrix <- regmatrix[, -zeroes]
    }
    yvector <- regmatrix[, "invsimpson"]
    xvars <- setdiff(colnames(regmatrix), "invsimpson")
    xmatrix <- regmatrix[, xvars]
    bmaresults <- runbmac(x = xmatrix, y = yvector, defaultmode = "gaussian",
                          modelsizearray = c(1, 3, 5), iter = bmaiterations, seed = 101)
    ## Save the model results
    saveRDS(object = bmaresults, file = "./output/lungbiome_linearreg_therapy.RDS")
    return(bmaout)
}
linearregresults <- runlinearbma(regmatrix = regmatrix_therapy, bmaiterations = 80000)


## At the end of this script, you should have the following files:
## ./output/lungbiome_MLR_diagnoses.RDS
## ./output/lungbiome_linearreg_therapy.RDS

##### Next scripts: Extract model findings and evaluate model #####
## Extract model findings: lungbiome_modelextraction.R
## Evaluate model performance: lungbiome_evaluation.R