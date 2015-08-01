##*****************************************************************************************************
##  For complete background and details, please refer to:
##  Shankar, J. et al. Microbiome and cytokine signatures of bacterial pneumonia
##  and tracheobronchitis.(Manuscript under review) (2015).
##
##  Analysis was performed under the following environment:
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
## This script provides the code for performing univariate comparisons.
## ***********************************************************************************

##### 0. Load required libraries and files #####
## clone the GitHub directory and set it as the working directory
## setwd('/path_to_GitHub_clone')
## Create a directory 'output' for saving output from the analysis
dir.create("./output", showWarnings = T, recursive = T)

##### 1. Load libraries and project data #####
source("./lungbiome_utilities.R")
source("./lungbiome_loadprojectdata.R")

##### 2. Welch's statistic and bootstrap CI function #####
generatebootwelch <- function(whichdata, comparevar, splitonvar, bootstrapit, numberofbootstraps) {
    datasub <- whichdata[, c(comparevar, splitonvar)]
    ## get these rows from dataset
    getwelchstats <- function(bootstrapit) {
        if (bootstrapit == T) {
            ## resample indices of the rows with replacement
            samplerows <- sample(1:nrow(datasub), size = nrow(datasub), replace = T)
        } else {
            samplerows <- 1:nrow(datasub)
        }
        sampled_data <- datasub[samplerows, ]
        refdata <- data.frame(sampled_data[which(sampled_data[, splitonvar] == 0),
            comparevar])
        testdata <- data.frame(sampled_data[which(sampled_data[, splitonvar] == 1),
            comparevar])
        colnames(testdata) <- comparevar
        colnames(refdata) <- comparevar
        ## Welch's statistic
        welchvector <- sapply(comparevar, function(variablename) {
            testvector <- testdata[, variablename]
            refvector <- refdata[, variablename]
            if (length(testvector) == 0 | length(refvector) == 0) {
                return(NA)
            } else {
                welchstat <- getwelch(vectofnum1 = refvector, vectofnum2 = testvector)
                return(welchstat)
            }
        }, simplify = T)
        return(welchvector)
    }
    if (bootstrapit == T) {
          welchboot <- mclapply(1:numberofbootstraps, function(num) {
            out <- getwelchstats(bootstrapit = T)
            return(out)
        }, mc.cores = 8)
        welchstat <- ldply(welchboot)
    } else {
        welchstat <- getwelchstats(bootstrapit = F)
    }
    return(welchstat)
}

diagabbr <- c("PNA", "COL", "TBr")
names(diagabbr) <- c(2, 0, 1)

getbootci <- function(whichdata, whichvar, acrosswhat,
                      numreplicates, numcomparisons, filename) {
    if (filename %in% cytokines) {
        whichdata[, acrosswhat] <- diagabbr[as.character(whichdata$diagnosis_simple_code)]
    }
    combos <- combn(x = unique(whichdata[, acrosswhat]), m = 2)
    bootcilist <- sapply(1:ncol(combos), function(combocol) {
        dsub <- whichdata[whichdata[, acrosswhat] %in% combos[, combocol], ]
        contrastname <- paste(unique(dsub[, acrosswhat]), collapse = "_")
        dsub[, acrosswhat] <- if ("COL" %in% dsub[, acrosswhat]) {
            ifelse(dsub[, acrosswhat] == "COL", 0, 1)
        } else {
            ifelse(dsub[, acrosswhat] == "TBr", 0, 1)
        }
        welchvector <- generatebootwelch(whichdata = dsub, comparevar = whichvar,
                                         splitonvar = acrosswhat, bootstrapit = T,
                                         numberofbootstraps = numreplicates)
        welchreal <- generatebootwelch(whichdata = dsub, comparevar = whichvar,
                                       splitonvar = acrosswhat, bootstrapit = F,
                                       numberofbootstraps = numreplicates)
        bootci <- quantile(welchvector[, whichvar],
                           prob = c(0.025/numcomparisons, 1 - (0.025/numcomparisons)),
                           na.rm = T)
        return(c(contrastname, welchreal, bootci))
    }, simplify = F)
    bootply <- ldply(bootcilist)
    bootply$significant <- ifelse((as.numeric(bootply[, 3]) * as.numeric(bootply[, 4])) > 0, "yes", "no")
    whichci <- (1 - (2 * 0.025/numcomparisons)) * 100
    colnames(bootply) <- c("contrast",
                           whichvar,
                           paste("lower", whichci, sep = "_"),
                           paste("upper", whichci, sep = "_"),
                           "significant")
    saveRDS(bootply, sprintf("./output/lungbiome_univariate_bootstraptests_%s.RDS", filename))
}

##### 4. Univariate analysis of diversity #####

## Are diversity values significantly different across colonization, pneumonia and tracheobronchitis?
## Table 3 in the manuscript
getbootci(whichdata = countmatrix, whichvar = "invsimpson", acrosswhat = "diagnosis_simple",
          numreplicates = 20000, numcomparisons = 1, filename = "diversity")


##### 5. Univariate analysis of cytokines #####

## Which cytokines are significantly different across colonization, pneumonia and tracheobronchitis?
## Table 4 in the manuscript
sapply(cytokines, function(ckine) {
    cat(ckine, "\n")
    getbootci(whichdata = regmatrix_diagnoses, whichvar = ckine, acrosswhat = "diagnosis_simple",
        numreplicates = 20, numcomparisons = 1, filename = ckine)
}, simplify = F)


## At the end of this script, you should have the following files:
## ./output/lungbiome_univariate_bootstraptests_diversity.RDS (Table 3)
## ./output/lungbiome_univariate_bootstraptests_{cytokine_name}.RDS (Table 4)