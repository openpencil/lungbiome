#*****************************************************************************************************
##  For complete background and details, please refer to:
##  Shankar J. et al. Looking beyond respiratory cultures: Microbiome-cytokine signatures of bacterial
##  pneumonia and tracheobronchitis in lung transplant recipients. Am J Transplant. Wiley Online Library;
##  2015; Available from: http://dx.doi.org/10.1111/ajt.13676
##
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
## This script provides the common utility functions for all the other scripts.
## **********************************************************************************

##### 1. Install and load required libraries #####
installloadpkgs <- function(pkgs) {
    ## check packages installed
    pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
    if (length(pkgs_miss) > 0) {
        install.packages(pkgs_miss)
    } else if (length(pkgs_miss) == 0) {
        message("\n ...Packages were already installed!\n")
    }
    ## load packages not already loaded:
    attached <- search()
    attached_pkgs <- attached[grepl("package", attached)]
    need_to_attach <- pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]
    if (length(need_to_attach) > 0) {
        ## alternative to library
        for (i in 1:length(need_to_attach)) {
            require(need_to_attach[i], character.only = TRUE)
        }
    }
    if (length(need_to_attach) == 0) {
        message("\n ...Packages were already loaded!\n")
    }
}

##### 2. Global library options #####
installloadpkgs(c("rlecuyer", "ggplot2", "reshape", "scales", "plyr",
                  "RColorBrewer", "data.table", "grid", "parallel"))

# Sets the random number generator
RNGkind("L'Ecuyer-CMRG")
options(useFancyQuotes = FALSE)

##### 3. Microbiome names cleanup functions #####
## Final variable names
cleanup <- function(invalue) {
    out <- unlist(sapply(1:length(invalue), function(rnum) {
        origvalue <- invalue[rnum]
        if (origvalue == "uncl_Bacteria_1") {
            genusval <- "Unclassified"
            phylaval <- "Bacteria"
            value <- sprintf("%s:{%s}", phylaval, genusval)
        } else if (origvalue %in% unclassdata$original) {
            genusval <- unclassdata[which(unclassdata$original == origvalue), "genus"]
            phylaval <- unclassdata[which(unclassdata$original == origvalue), "phylum"]
            value <- sprintf("%s:{%s}", phylaval, genusval)
        } else if (grepl("genus_incertae_sedis", origvalue)) {
            genusval <- "incertae sedis"
            phylaval <- treedata[which(treedata$genus == origvalue), "phylum"]
            value <- sprintf("%s:{%s}", phylaval, genusval)
        } else if (grepl("incertae_sedis", origvalue)) {
            genusval <- "incertae sedis"
            origvalue <- gsub("incertae_sedis", "genus_incertae_sedis", origvalue)
            phylaval <- treedata[which(treedata$genus == origvalue), "phylum"]
            value <- sprintf("%s:{%s}", phylaval, genusval)
        } else if (origvalue %in% treedata$genus) {
            genusval <- treedata[which(treedata$genus == origvalue), "genus"]
            phylaval <- treedata[which(treedata$genus == origvalue), "phylum"]
            value <- sprintf("%s:{%s}", phylaval, genusval)
        } else {
            value <- "Other taxa (<5%)"
        }
        valueout <- gsub("_Chloroplast", "", value)
        return(valueout)
    }))
    return(out)
}

getclusternum <- function(colofint) {
    colvarname <- colofint$varname
    indgen <- grep("uncl", colvarname)
    indnotuncl <- grep("uncl", colvarname, invert = T)
    colint <- colofint$cleanmicrobe
    colgen <- colint[indgen]
    if (length(indnotuncl) != 0) {
        ## there is a pre-classified genus present.
        out <- sapply(1:length(colgen), function(indcol) {
            ## start from [2]
            outformat <- gsub("}$", sprintf(" [%d]}", indcol + 1), colint[indcol])
            return(outformat)
        })
    } else {
        ## there is no pre-classified genus present.
        out <- sapply(1:length(colgen), function(indcol) {
            ## start from [1]
            outformat <- gsub("}$", sprintf(" [%d]}", indcol), colint[indcol])
            ## throw out the [1]
            outformat <- gsub("(.*)(\\s\\[1\\])(\\})$", "\\1\\3", outformat)
            return(outformat)
        })
    }
    names(out) <- as.character(indgen)
    outcol <- sapply(seq_along(colint), function(scol) {
        if (scol %in% indgen) {
            finalcol <- out[as.character(scol)]
        } else {
            finalcol <- colint[scol]
        }
    })
    return(outcol)
}

##### 4. Rescaling functions #####
ourrescale <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
    if (zero_range(from) || zero_range(to))
        return(x)
    out <- (x - from[1])/diff(from) * diff(to) + to[1]
}

oursignrescale <- function(value) {
    value <- unlist(value)
    negindex <- which(value < 0)
    posindex <- setdiff(1:length(value), negindex)
    neg <- value[negindex]
    pos <- value[posindex]
    negr <- ourrescale(neg, from = c(-max(abs(value), na.rm = T), 0), to = c(0, 0.5))
    names(negr) <- negindex
    posr <- ourrescale(pos, from = c(0, max(abs(value), na.rm = T)), to = c(0.5, 1))
    names(posr) <- posindex
    rval <- c(negr, posr)
    allval <- rval[as.character(c(1:length(value)))]
    names(allval) <- c()
    return(allval)
}

keepsignrescale <- function(value) {
    value <- unlist(value)
    negindex <- which(value < 0)
    posindex <- setdiff(1:length(value), negindex)
    neg <- value[negindex]
    pos <- value[posindex]
    negr <- ourrescale(neg, from = c(-max(abs(value), na.rm = T), 0), to = c(-1, 0))
    names(negr) <- negindex
    posr <- ourrescale(pos, from = c(0, max(abs(value), na.rm = T)), to = c(0, 1))
    names(posr) <- posindex
    rval <- c(negr, posr)
    allval <- rval[as.character(c(1:length(value)))]
    names(allval) <- c()
    return(allval)
}

##### 5. Favourite ggplot theme #####
lightertheme <- theme(panel.background = element_rect(fill = "#f5f5f4", colour = "#fbfbfc"),
                      panel.border = element_rect(colour = "grey80", linetype = "solid", fill = NA),
                      panel.grid.major = element_line(colour = "grey90", size = 0.08),
                      panel.grid.minor = element_line(colour = "grey90", size = 0.08),
                      panel.margin = unit(0.3, "lines"),
                      plot.background = element_rect(fill = "transparent"),
                      plot.margin = unit(c(1, 1, 1, 1), "mm"),
                      text = element_text(family = "Helvetica", colour = "black", size = 12),
                      plot.title = element_text(size = 12),
                      strip.background = element_rect(fill = "transparent", color = "#ffffff", size = 1),
                      strip.text.x = element_text(size = 12, colour = "black"),
                      strip.text.y = element_text(size = 12, colour = "black"),
                      axis.text.x = element_text(size = 12, colour = "black"),
                      axis.text.y = element_text(size = 12, colour = "black"),
                      legend.position = "left",
                      legend.background = element_rect(fill = "transparent"),
                      legend.background = element_blank(),
                      legend.key = element_blank(),
                      legend.key.size = unit(0.5, "cm"))
##### 6. Misc. functions #####
invsimpson <- function(row) {
    ## 1/\sum_{1}^{n}(p^2)
    props <- row/sum(row)
    denom <- sum(props^2)
    ## Reciprocal of the sum of the proportions squared.
    simpson <- 1/denom
    return(simpson)
    }

getwelch <- function(vectofnum1, vectofnum2) {
    ## remove NA
    one <- na.omit(vectofnum1)
    two <- na.omit(vectofnum2)
    wnum <- mean(two) - mean(one)
    wdenom <- sqrt((sd(two))^2/length(two) + (sd(one))^2/length(one))
    wtest <- wnum/wdenom
    return(wtest)
}
tosentencecase <- function(inputstring) {
    substr(inputstring, start = 1, 1) <- toupper(substr(inputstring, start = 1, 1))
    return(inputstring)
}

convertlogittoprob <- function(inputlogit) {
    odds <- exp(inputlogit)
    prob <- odds/(odds + 1)
    return(prob)
}
