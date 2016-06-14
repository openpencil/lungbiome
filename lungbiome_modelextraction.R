##*****************************************************************************************************
##  For complete background and details, please refer to:
##  Shankar J. et al. Looking beyond respiratory cultures: Microbiome-cytokine signatures of bacterial
##  pneumonia and tracheobronchitis in lung transplant recipients. Am J Transplant. Wiley Online Library;
##  2015; Available from: http://dx.doi.org/10.1111/ajt.13676
##
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
## This script provides the code for extracting summaries from the MLR and linear
## regression BMA models.
## ******************************************************************************************************

##### 1. Load required libraries and files #####
## clone the GitHub directory and set it as the working directory
## setwd('/path_to_GitHub_clone')
## Create a directory 'output' for saving output from the analysis
dir.create("./output", showWarnings = T, recursive = T)

##### 2. Load model results from the output directory #####
mlr_results <- readRDS(file = "./output/lungbiome_MLR_diagnoses.RDS")
linearreg_results <- readRDS(file = "./output/lungbiome_linearreg_therapy.RDS")

##### 3. Extract MLR model results #####
extract_mlr_results <- function(modelinfo, vector_of_continuousvars) {
    xmatrix <- modelinfo$xmatrix
    yvector <- modelinfo$yvector
    allbetas <- modelinfo$model$beta
    ## discard burn-in betas
    betanoburn <- allbetas[-(1:ceiling(0.1 * nrow(allbetas))), ]
    ## get betas for each contrast
    beta10 <- betanoburn[, grep(":1$", colnames(betanoburn))]
    beta20 <- betanoburn[, grep(":2$", colnames(betanoburn))]
    colnames(beta20) <- gsub(":2$", ":1", colnames(beta20))
    beta21 <- beta20[, colnames(beta10)] - beta10
    betalist <- list(b10 = beta10, b20 = beta20, b21 = beta21)
    allstats <- sapply(names(betalist), function(bname) {
        betavals <- betalist[[bname]]
        zeroval <- gsub("b(\\d{1})(\\d{1})", "\\2", bname)
        diseaseval <- gsub("b(\\d{1})(\\d{1})", "\\1", bname)
        l1index <- which(yvector == zeroval)
        l2index <- which(yvector == diseaseval)
        ## remove the intercept
        intindex <- grep("ntercept", colnames(betavals))
        betaci <- sapply(colnames(betavals)[-intindex], function(bcol) {
            ## include spike at 0 -- data after discarding burn-in is included
            cival <- c(quantile(x = betavals[, bcol], probs = c(0.025, 0.975)),
                       meanbeta = mean(betavals[, bcol]))
            names(cival) <- c("lowerci", "upperci", "meanbeta")
            bcolmod <- gsub(":\\d{1}$", "", bcol)
            if (bcolmod %in% continuousvars) {
                x1 <- mean(xmatrix[l1index, bcolmod])
                x2 <- mean(xmatrix[l2index, bcolmod])
                x1vect <- xmatrix[l1index, bcolmod]
                x2vect <- xmatrix[l2index, bcolmod]
                ## get welch statistics for continuous variable
                wnum <- x2 - x1
                wdenom <- sqrt((sd(x2vect))^2/length(x2vect) + (sd(x1vect))^2/length(x1vect))
                wtest <- wnum/wdenom
            } else {
                ## for categorical variables
                x2 <- median(xmatrix[l2index, bcolmod])
                x1 <- median(xmatrix[l1index, bcolmod])
                wtest <- NA
            }
            ## calculate average effect of the change in the microbe/cytokine level
            ##  on response probability = beta x (mean(Xp) - mean(Xc))
            betax <- cival["meanbeta"] * (x2 - x1)
            names(betax) <- "betax"
            ## compute pip across 72000 combinations
            ip <- length(which(betavals[, bcol] > 0))/nrow(betavals)
            out <- c(cival, betax, pip = ip, welch = wtest)
            return(out)
        })
        outstats <- data.frame(t(betaci))
        outstats$variable <- gsub(":\\d{1}$", "", rownames(outstats))
        return(outstats)
    }, simplify = F)
    return(allstats)
}

modelresults <- extract_mlr_results(modelinfo = mlr_results, vector_of_continuousvars = continuousvars)

processresultsforplotting <- function(contrastdata) {
    ## Reduce to only significant variables
    significant <- which(contrastdata$lowerci * contrastdata$upperci >= 0)
    ## upper or lower CI are not zero
    signonzero <- which(contrastdata$lowerci != 0 | contrastdata$upperci != 0)
    sigintersect <- rownames(contrastdata)[intersect(significant, signonzero)]
    contrastdata[sigintersect, "sig"] <- "sig"
    contrastdata[setdiff(rownames(contrastdata), sigintersect), "sig"] <- "signot"
    contrastdata$group <- ""
    contrastdata$group[contrastdata$variable %in% cytokines] <- "Cytokine"
    contrastdata$group[contrastdata$variable %in% microbiome] <- "Microbiome"
    contrastdata$group[grep("acr", contrastdata$variable)] <- "ACR"
    contrastdata$group[grep("ltx", contrastdata$variable)] <- "Time since LTx"
    contrastdata$group[grep("age", contrastdata$variable)] <- "Age"
    contrastdata$combo <- paste(contrastdata$group, contrastdata$sig, sep = "")
    contrastdata$ipcent <- contrastdata$pip * 100
    return(contrastdata)
}

modeldataforplotting <- sapply(names(modelresults), function(contrastname) {
    contrastdata <- modelresults[[contrastname]]
    out <- processresultsforplotting(contrastdata = contrastdata)
    return(out)
}, simplify = F)

saveRDS(object = modeldataforplotting, file = "./output/lungbiome_MLR_modeldataforplotting.RDS")

#### 4. Extract linear regression results #####
extract_linearreg_results <- function(listofmodels, getmodeldiagnostics) {
    modelinfo <- sapply(names(listofmodels), function(filenamestring) {
        modeldata <- listofmodels[[filenamestring]]$model$model
        cvems <- listofmodels[[filenamestring]]$ems
        ## cat(filenamestring, '\n')
        burntenpercent <- 0.1 * nrow(modeldata$beta)
        iprob <- summary(modeldata, burn = burntenpercent)
        resid <- iprob$residual.sd
        rsq <- iprob$rsquare
        if (getmodeldiagnostics == T) {
            summarydata <- list(residualsd = resid, rsquare = rsq)
        } else {
            probcoeff <- iprob$coefficients
            allvars <- setdiff(rownames(probcoeff), "(Intercept)")
            ## CIs of the posterior draws for all variables.
            getci <- sapply(allvars, function(var) {
                ## get remaining betas
                betavector <- modeldata$beta[-c(1:burntenpercent), var]
                ## eliminate zeroes
                betaslab <- betavector[which(betavector != 0)]
                betaci <- quantile(x = betaslab, probs = c(0.025, 0.5, 0.975))
                ## with zeroes
                betacizero <- quantile(x = betavector, probs = c(0.025, 0.5, 0.975))
                inclprob <- probcoeff[var, "inc.prob"]
                out <- c(betaci, inclprob, betacizero)
                names(out) <- c("lowerci", "median", "upperci", "pip", "lowerci_z",
                  "median_z", "upperci_z")
                return(out)
            })
            summarydata <- data.table(t(getci), colnames(getci), rep(filenamestring,
                ncol(getci)))
            setnames(summarydata, colnames(summarydata)[c(8, 9)], c("variable", "filename"))
            summarydata$ems <- cvems
        }
        return(summarydata)
    }, simplify = F)
    return(modelinfo)
}

linearreg_therapy_diagnostics <- extract_linearreg_results(listofmodels = list(linearregresults = linearreg_results),
                                                           getmodeldiagnostics = T)
saveRDS(object = linearreg_therapy_diagnostics, file = "./output/lungbiome_linearreg_therapy_diagnostics.RDS")

linearreg_therapy_results <- extract_linearreg_results(listofmodels = list(linearregresults = linearreg_results),
                                                       getmodeldiagnostics = F)
saveRDS(object = linearreg_therapy_results, file = "./output/lungbiome_linearreg_therapy_results.RDS")

## At the end of this script, you should have the following files:
## ./output/lungbiome_MLR_modeldataforplotting.RDS
## ./output/lungbiome_linearreg_therapy_diagnostics.RDS
## ./output/lungbiome_linearreg_therapy_results.RDS

##### Next script: lungbiome_modelvisualization.R ######
