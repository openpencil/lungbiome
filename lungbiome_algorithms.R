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
## This script provides the underlying code for building and crossvalidating the BMA models.
## ******************************************************************************************************

##### 1. Load required libraries #####
## clone the GitHub directory and set it as the working directory
## setwd('/path_to_GitHub_clone')
source("./lungbiome_utilities.R")
installloadpkgs(c("doMC", "BoomSpikeSlab"))

##### 2. Set number of cores to use in analysis #####
numcore <- 4
## register cores for parallel processing
registerDoMC(cores = numcore)

##### 3. Main BMA modeling framework #####
runspike <- function(expected_model_size,
                     y,
                     x,
                     iterations = 10000,
                     defaultmode = "gaussian",
                     seed = NULL) {
    ## the reference level is always 0 for categorical response.
    if (typeof(y) == "character") {
        y <- as.integer(as.factor(y)) - 1
    }
    x <- as.matrix(x)
    if (defaultmode == "gaussian") {
        ## linear regression
        spikesmodel <- lm.spike(y ~ x, niter = iterations, expected.model.size = expected_model_size,
                                seed = seed)
    } else if (defaultmode == "logistic") {
        ## logistic regression
        spikesmodel <- logit.spike(y ~ x, niter = iterations, expected.model.size = expected_model_size,
                                   seed = seed)
    }
    ## generate model summary. 10% iterations discarded as burn-in
    spikesummary <- summary(spikesmodel, burn = round(0.1 * iterations, 0))
    betas <- spikesummary$coefficients[, "mean"]
    inclprob <- spikesummary$coefficients[, "inc.prob"]
    names(betas) <- gsub("^xscaled|^x", "", names(betas))
    return(list(betas = betas, inclprob = inclprob, model = spikesmodel, summary = spikesummary,
                xmatrix = x, yfactor = y))
}

##### 4. Get crossvalidated model size for BMA #####
crossvalidate <- function(x,
                          y,
                          folds,
                          seed,
                          iterations = 10000,
                          defaultmode = "gaussian",
    modelsizearray = c(1, 3, 7, 10)) {
    set.seed(seed)
    ## assign fold to each row in x matrix
    foldid <- sample(rep(x = (1:folds), times = (nrow(x)/folds)), replace = F)
    ## first column is the response variable
    xy <- cbind(y, x)
    ## assign test/train partitions
    assigntesttrain <- function(fold) {
        testset <- xy[which(foldid == fold), ]
        trainset <- xy[which(foldid != fold), ]
        out <- list(trainset, testset)
        names(out) <- c("train", "test")
        return(out)
    }
    ## list of test/train data subsets. Number of test/train sets == number of folds
    testtrain <- lapply(1:folds, assigntesttrain)
    ## run BMA on the training set for all model sizes
    train_path <- mclapply(testtrain, function(elementlist) {
        xscaled <- as.matrix(elementlist$train[, -1])
        y <- elementlist$train[, 1]
        modelout <- mclapply(modelsizearray, function(num_vars) {
            if (defaultmode == "gaussian") {
                spikesmodel <- lm.spike(y ~ xscaled, niter = iterations, expected.model.size = num_vars)
            } else if (defaultmode == "logistic") {
                spikesmodel <- logit.spike(y ~ xscaled, niter = iterations, expected.model.size = num_vars)
            }
            return(spikesmodel)
        }, mc.cores = 3)
        return(modelout)
    }, mc.cores = 4)

    modelerrors <- sapply(1:length(modelsizearray), function(feature) {
        mseval <- sapply(1:folds, function(fold) {
            if (defaultmode == "gaussian") {
                ## get median predictions on test set using the trained linear regression model
                predicttestreg <- predict(object = train_path[[fold]][[feature]],
                  newdata = as.matrix(testtrain[[fold]]$test[, -1]), type = "response",
                  burn = round(0.1 * iterations, 0))
                medianpred <- apply(predicttestreg, 1, median)
                y_test <- testtrain[[fold]]$test[, 1]
                ## mean square error (or L2 distance) between test response and prediction
                mse_test <- mean((y_test - medianpred)^2)
            } else if (defaultmode == "logistic") {
                ## get median predictions on test set using the trained linear regression model
                predicttestreg <- predict(object = train_path[[fold]][[feature]],
                  newdata = as.matrix(testtrain[[fold]]$test[, -1]), type = "prob",
                  burn = round(0.1 * iterations, 0))
                ## Determine response prediction
                medianpred <- apply(predicttestreg, 1, function(x) {
                  class1 <- length(which(x >= 0.5))
                  class0 <- length(which(x < 0.5))
                  voting <- ifelse(class1 >= class0, 1, 0)
                  return(voting)
                })
                y_test <- testtrain[[fold]]$test[, 1]
                ## classification error: test response does not equal prediction
                mse_test <- length(which(y_test != medianpred))
            }
            return(mse_test)
        }, simplify = T)
        ## mean MSE or mean classification error over all the folds
        avgmse <- mean(mseval)
        return(avgmse)
    }, simplify = T)
    ## identify index of the model size with the minimum MSE or classification error
    best_feature_index <- which(modelerrors == min(modelerrors))[1]
    cat(sprintf("The crossvalidated model size is %s.\n", modelsizearray[best_feature_index]))
    ## for lambda, return the index
    best_feature <- modelsizearray[best_feature_index]
    return(best_feature)
}

##### 5. BMA with crossvalidation for expected model size #####
runbmac <- function(x, y, defaultmode = "gaussian", modelsizearray = c(1, 3, 7, 10),
    iter = 10000, seed = NULL) {
    if (nrow(x) < 25) {
        nfolds <- 3
    } else {
        nfolds <- 5
    }
    if (length(modelsizearray) > 0) {
        bestmodelsize_bmac <- crossvalidate(x, y, folds = nfolds, seed = seed, iterations = iter,
                                            defaultmode = defaultmode, modelsizearray = modelsizearray)
    } else {
        ## if crossvalidation is not turned on, the sparsest model will be estimated.
        bestmodelsize_bmac <- 1
    }
    ## run BMAC with crossvalidated model size
    model_bmac <- runspike(expected_model_size = bestmodelsize_bmac, y = y, x = x,
                           iterations = iter, defaultmode = defaultmode, seed = seed)
    return(list(ems = bestmodelsize_bmac, model = model_bmac, yvector = y, xmatrix = x))
}

##### 6. Multinomial logistic regression (MLR) model #####
runmlrbma <- function(x, y, addtomodelsize = 1, iter = 10000, seed = NULL) {
    dframe <- data.frame(x, y)
    ## expected.subject.model.size should be more than the number of levels in the model
    modelsize <- length(unique(y)) + addtomodelsize
    ## BoomSpikeSlab has three possible updates for MCMC sampling.
    ## 'DA', 'RWM' and 'TIM' are the probabilities for these three updates.
    ## The values for these probabilities come from Scott, 2011, doi:10.1007/s00362-009-0205-0
    model_bma <- mlm.spike(y ~ ., niter = iter, expected.subject.model.size = modelsize,
                           data = dframe, proposal.weights = c(DA = 0.8, RWM = 0.1, TIM = 0.1))
    sumbma <- summarizebma(model_bma)
    return(list(ems = modelsize, model = model_bma, modelsummary = sumbma,
                yvector = y, xmatrix = x))
}

## summary is not implemented yet for MLR models. This function is an alternative
summarizebma <- function(modelresults) {
    stats <- sapply(colnames(modelresults$beta), function(col) {
        burn_number <- ceiling(0.1 * nrow(modelresults$beta))
        ip <- mean(modelresults$beta[-c(1:burn_number), col] != 0)
        beta <- mean(modelresults$beta[-c(1:burn_number), col])
        betaq <- quantile(modelresults$beta[-c(1:burn_number), col], c(0.025, 0.975))
        return(c(ip = ip, beta = beta, betaq))
    }, simplify = T)
    stats <- t(stats)
    stats <- stats[order(stats[, 1], decreasing = T), ]
}