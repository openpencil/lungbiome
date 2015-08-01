#*****************************************************************************************************
##  For complete background and details, please refer to:
##  Shankar, J. et al. Microbiome and cytokine signatures of bacterial pneumonia
##  and tracheobronchitis.(Manuscript under review) (2015).
##
##  Modeling and visualization was performed under the following environment:
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
## This script provides the code for evaluating the classification performance of the MLR-BMA model.
## ****************************************************************************************************

##### 0. Creates output directory #####
dir.create("./output", showWarnings = T, recursive = T)

##### 1. Load project data from the output directory #####
source(file = "./lungbiome_loadprojectdata.R")
source(file = "./lungbiome_algorithms.R")

##### 2. Estimate MLR BMA on over 1200 bootstrap resamples of the study data. #####
mlrbootstrap <- function(xvars, responsevar, iternum, regdata, bootstrapnum) {
    yvar <- regdata[, responsevar]
    xvar <- regdata[, xvars]
    ## Eliminate variables that are zero throughout or have a constant value in all columns.
    zeroes <- which(colSums(xvar) == 0)
    if (length(zeroes) > 0) {
        xvar <- xvar[, -zeroes]
    }
    ## if bootstrapnum is set to 0, then the function estimates the model on the data
    ## without bootstrap resampling
    if (bootstrapnum == 0) {
        resampledrows <- 1:nrow(xvar)
    } else {
        resampledrows <- sample(x = 1:nrow(xvar), size = nrow(xvar), replace = T)
    }
    ## Check to see that the response has the requisite number of categories and observations
    ycheck <- as.data.frame.table(table(yvar[resampledrows]), stringsAsFactors = F)
    if ((any(unique(ycheck$Freq) < 2) | length(unique(ycheck$Var1)) < 3) == F) {
        ## Set expected model size to the lowest that works. 8 for 3 levels.
        ## (addtomodelsize = 5 + 3 levels of outcome)
        xvarboot <- xvar[resampledrows, ]
        yvarboot <- yvar[resampledrows]
        zeroboot <- which(colSums(xvarboot) == 0)
        if (length(zeroboot) > 0) {
            xvarboot <- xvarboot[, -zeroboot]
        }
        bmaresults <- runmlrbma(x = xvarboot, y = yvarboot, addtomodelsize = 5,
                                iter = iternum, seed = 101)
        saveRDS(object = bmaresults,
                file = sprintf("./output/lungbiome_MLR_diagnoses_bootstrapnum%s.RDS", bootstrapnum))
    }
}

covariables_in_mlrmodel <- c(microbiome, cytokines, ltxindications,
                             therapyvars, clinical, patients)

## Run MLR-BMA model on all 1200 bootstrap resamples, and on the original data (at
## bootstrapnum=0) Note: This is a time-consuming step. Each model takes around 7-8
## minutes to run.
sapply(c(0:1200), function(bnum) {
    mlrbootstrap(xvars = covariables_in_mlrmodel,
                 responsevar = "diagnosis_simple_code",
                 iternum = 80000,
                 regdata = regmatrix_diagnoses,
                 bootstrapnum = bnum)
    }, simplify = F)

## At the end of this computation, you will have 1200 MLR-BMA model results in the output folder:
## lungbiome_MLR_diagnoses_bootstrapnum1.RDS
## lungbiome_MLR_diagnoses_bootstrapnum2.RDS
## lungbiome_MLR_diagnoses_bootstrapnum3.RDS ... and so on.
## You will also have the MLR-BMA model run on original data (no bootstrap resample)
## lungbiome_MLR_diagnoses_bootstrapnum0.RDS

## Save the coefficients generated from the bootstrap resamples ##
bootstrapmodelfilenames <- list.files(path = "./output", pattern = ".*bootstrapnum.*.RDS",
                                      all.files = T, full.names = T, recursive = T)
originalmodelname <- "./output/lungbiome_MLR_diagnoses_bootstrapnum0.RDS"
bootstrapmodelfiles <- setdiff(bootstrapmodelfilenames, originalmodelname)
bootstrapcoefficients <- sapply(bootstrapmodelfiles, function(fname) {
    out <- readRDS(file = fname)
    ## Discard the model. We need just the model coefficients
    out$model <- NULL
    return(out)
    }, simplify = F)
names(bootstrapcoefficients) <- sprintf("bootstrapcoefficients_%04d",
                                        as.numeric(gsub(".*_bootstrapnum(\\d{1,}).RDS",
                                                        "\\1", names(bootstrapcoefficients))))
saveRDS(object = bootstrapcoefficients, file = "./output/bootstrapcoefficientsforauc.RDS")


##### 3. Compute Optimism #####
calculateauc <- function(mlr_model_data, betacoefficients, thresholdinterval) {
    ## covariate matrix
    covariables_x <- mlr_model_data$xmatrix
    ## response vector
    response_y <- mlr_model_data$yvector
    ## assign intercept to 1 (originally this is NULL)
    covariables_x$`(Intercept)` <- 1
    ## isolate regression coefficients for the 1st contrast
    colint1 <- grep(":1$", rownames(betacoefficients))
    betaint1 <- betacoefficients[colint1, "beta"]
    ## isolate regression coefficients for the 1st contrast
    colint2 <- grep(":2$", rownames(betacoefficients))
    betaint2 <- betacoefficients[colint2, "beta"]
    ## strip out the numbers and make names similar
    names(betaint1) <- gsub(":\\d{1}", "", names(betaint1))
    names(betaint2) <- gsub(":\\d{1}", "", names(betaint2))
    
    ## compute auc for each contrast in the models
    sapply(c("pc", "tc", "pt"), function(whichcontrast) {
        if (whichcontrast == "pc") {
            rowint <- which(response_y == 2 | response_y == 0)
            betaint <- betaint2
        } else if (whichcontrast == "tc") {
            rowint <- which(response_y == 1 | response_y == 0)
            betaint <- betaint1
        } else {
            rowint <- which(response_y == 1 | response_y == 2)
            betaint3 <- betaint2[names(betaint1)] - betaint1
            betaint <- betaint3
        }
        ## compare covariable names in the regression dataset and beta coefficients
        if (all(names(betaint) %in% colnames(covariables_x))) {
            xdata <- covariables_x[rowint, names(betaint)]
            ydata <- response_y[rowint]
            diagnosis_interest <- which(ydata == max(ydata))
            diagnosis_reference <- which(ydata == min(ydata))
            
            ## calculate log odds of the response
            betax <- as.matrix(xdata[, names(betaint)]) %*% as.matrix(betaint)
            
            ## Convert logit to probability
            probresponse <- plogis(betax)
            equallyspacedxaxis <- seq(from = 0, to = 1, by = thresholdinterval)
            
            ## calculate metrics at each give probability threshold
            metrics <- sapply(rev(equallyspacedxaxis), function(thresholdp) {
                ## falsepositive belongs to the diagnosis_reference group but has been falsely
                ## classified to the diagnosis_interest group
                falsepositives <- which(which(probresponse > thresholdp) %in% diagnosis_reference)
                ## compute false positive rate
                fpr <- length(falsepositives)/length(diagnosis_reference)
                ## truenegative belongs to the diagnosis_reference group and has been rightly
                ## classified to this group.
                truenegatives <- which(which(probresponse <= thresholdp) %in% diagnosis_reference)
                ## compute true positive rate
                tnr <- length(truenegatives)/length(diagnosis_reference)
                ## truepositive belongs to the diagnosis_interest group and has been rightly
                ## classified to this group.  tpr is the same as sensitivity and recall
                truepositives <- which(which(probresponse > thresholdp) %in% diagnosis_interest)
                tpr <- length(truepositives)/length(diagnosis_interest)
                out <- c(TPR = tpr, FPR = fpr, TNR = tnr)
                return(out)
            }, simplify = T)
            calcauc <- sapply(2:ncol(metrics), function(colindex) {
                ## compute auc area of trapezium = 1/2 * base * sum of || sides
                ## x=fpr, y=tpr
                ## since we are moving in the opposite direction of the x-axis,
                ## using TNR instead of FPR (note: TNR is the complement of FPR)
                auc <- 0.5 * abs((metrics["TNR", colindex] - metrics["TNR", colindex - 
                  1]) * (metrics["TPR", colindex] + metrics["TPR", colindex - 1]))
            })
            aucvalue <- sum(calcauc)
            return(data.table(t(metrics), auc = aucvalue))
        } else {
            stop("For contrast:", whichcontrast, ",",
                 "covariable names in regression dataset do not\n
                 match with the names of the beta coefficients!")
        }
    }, simplify = F)
}

original_data <- readRDS("./output/lungbiome_MLR_diagnoses_bootstrapnum0.RDS")

getoptimism <- sapply(names(bootstrapcoefficients), function(bootstrapiteration, tnum = 0.01) {
    ## Compute the AUC using the regression coefficients estimated in Step 2,3 on the
    ## data from the bootstrap resamples.
    auc_on_bootstrap_resamples <- calculateauc(mlr_model_data = bootstrapcoefficients[[bootstrapiteration]],
                                               betacoefficients = bootstrapcoefficients[[bootstrapiteration]]$modelsummary,
                                               thresholdinterval = tnum)
    ## cat('AUC on bootstrap resamples computed.\n')
    ## Compute the AUC using the regression coefficients estimated in Steps 2,3 on the original data
    ## cat(bootstrapiteration,'\n')
    auc_on_original_data <- calculateauc(mlr_model_data = original_data,
                                         betacoefficients = bootstrapcoefficients[[bootstrapiteration]]$modelsummary,
                                         thresholdinterval = tnum)
    ## cat('AUC on original data computed.\n')
    ## Calculate the optimism for this model by subtracting the AUC computed on original data from
    ## the AUC computed from bootstrap resamples
    optimism <- sapply(c("pc", "tc", "pt"), function(contrast) {
        out <- unique(auc_on_bootstrap_resamples[[contrast]]$auc - auc_on_original_data[[contrast]]$auc)
        return(out)
        })
    return(optimism)
    }, simplify = F)

##### 4. Compute average optimism over bootstrap resamples #####
optimismdata <- ldply(getoptimism)
average_optimism <- sapply(c("pc", "tc", "pt"), function(cname) {
    out <- mean(optimismdata[, cname])
    return(out)
    })

##### 5. Compute orginal AUC #####
auc_original <- calculateauc(mlr_model_data = original_data,
                             betacoefficients = original_data$modelsummary,
                             thresholdinterval = 0.01)
saveRDS(object = auc_original, file = "./output/lungbiome_originalauc.RDS")

##### 6. Estimate final optimism-corrected AUC #####
auc_final_values <- sapply(c("pc", "tc", "pt"), function(cname) {
    auc_original_value <- unique(auc_original[[cname]]$auc)
    auc_optimism_corrected <- auc_original_value - average_optimism[cname]
    return(c(auc_original_value, auc_optimism_corrected))
    })

##### 7. Plot the ROC curves ##### Design and labels ##
linetypecodes <- c("dotdash", "F1")
names(linetypecodes) <- c("pc", "tc")

contrastnames <- c("Pneumonia vs. colonization", "Tracheobronchitis vs. colonization")
contrastvalues <- c("pc", "tc")
names(contrastnames) <- contrastvalues

renamelabels <- function(inlabels) {
    outlabels <- contrastnames[inlabels]
    return(outlabels)
}

aucdata <- ldply(auc_original)
plotroc <- function(aucplottingdata) {
    aucplot <- aucdata[which(aucdata$.id != "pt"), ]
    p <- ggplot(aucplot, aes(x = FPR, y = TPR, linetype = as.factor(.id)))
    p <- p + geom_abline(aes(intercept = 0, slope = 1), colour = "grey80")
    p <- p + geom_path()
    p <- p + lightertheme
    p <- p + scale_linetype_manual(values = linetypecodes,
                                   name = "Multivariable signatures",
                                   labels = renamelabels)
    p <- p + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")
    p <- p + theme(axis.title.x = element_text(size = 10, family = "Helvetica", colour = "black"),
                   axis.title.y = element_text(size = 10, family = "Helvetica", colour = "black"),
                   axis.text.x = element_text(size = 10, family = "Helvetica", colour = "black"),
                   axis.text.y = element_text(size = 10, family = "Helvetica", colour = "black"),
                   legend.position = c(0.7, 0.25))
    p <- p + guides(linetype = guide_legend(name = "Multivariable signatures", keyheight = 2.5),
                    label.theme = element_text(size = 10, angle = 0),
                    title.theme = element_text(size = 10, angle = 0, face = "bold"), ncol = 1)
    ggsave("./output/lungbiome_roc_curve.pdf", plot = p,
           width = 5, height = 5, units = "in", limitsize = F)
}

plotroc(aucplottingdata = aucdata) 
