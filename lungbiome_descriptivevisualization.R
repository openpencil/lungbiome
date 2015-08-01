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
## This script provides the code for descriptive data visualization.
## ***********************************************************************************

##### 0. This script will need an 'output' directory for saving plots #####

##### 1. Load environments and data #####
source("./lungbiome_utilities.R")
source("./lungbiome_loadprojectdata.R")

##### 2. Select taxa more than a specified threshold proportions #####
getgraphtaxa <- function(dataprop, threshold) {
    ## taxa that are in more than the set threshold relative abundance in at least one
    ## sample
    onlyprop <- dataprop[, setdiff(colnames(dataprop), "sample")]
    okbysample <- unlist(sapply(rownames(onlyprop), function(pname) {
        oktaxa <- colnames(onlyprop)[which(onlyprop[pname, ] >= threshold)]
        return(oktaxa)
    }))
    finaltaxa <- unique(okbysample)
    return(finaltaxa)
}

finaltaxa <- getgraphtaxa(dataprop = propmatrix, threshold = 0.05)

##### 3. Make stacked histograms and associated diversity plots #####

## Function for preparing data for graphing ##
preparedata <- function(dataprop, datacount, diagnosisgroup, finaltaxa) {
    dataprop$patient <- sprintf("P%02d", as.numeric(gsub("(.*)_(.*)", "\\1", dataprop$sample)))
    ## get the patient rows of interest
    diagnosispatients <- patientsets[[diagnosisgroup]]
    patientprop <- dataprop[which(dataprop$patient %in% diagnosispatients), ]
    ## identify taxa present in less than 0.05 proportion
    othertaxa <- setdiff(colnames(patientprop), c(finaltaxa, "sample", "patient"))
    ## collapse all the non threshold taxa into a single variable called others'
    others <- rowSums(patientprop[, othertaxa])
    ## put the final taxa to be plotted together with the others
    finaldata <- data.frame(patientprop[, finaltaxa], others = others[rownames(patientprop)])
    ## quick check. rowSums(finaldata) should be all 1
    finaldata$sample <- rownames(finaldata)
    finaldata$patient <- sprintf("P%02d", as.numeric(gsub("(.*)_.*", "\\1", finaldata$sample)))
    ## Get timesinceltx and diversity column from counts data
    rownames(datacount) <- datacount$sample
    ggdata <- data.frame(finaldata,
                         datacount[rownames(finaldata), c("timesinceltx", "diagnosis_simple", "invsimpson")],
                         stringsAsFactors = F)
    ggmelt <- melt(ggdata, measure.vars = c(finaltaxa, "others"))
    ## reconverts the factor variable into character
    ggmelt$variable <- as.character(ggmelt$variable)
    ggmelt$cleanvar <- cleanup(ggmelt$variable)
    ggmelt$genusonly <- gsub(".*\\{(.*)\\}", "\\1", ggmelt$cleanvar)
    ggmelt$timesinceltxround <- sprintf("%02.1f", ggmelt$timesinceltx)
    return(ggmelt)
}

## Function for making stacked histograms ##
plotstackedhistograms <- function(dataforplotting, diagnosisgroup, makelegend) {
    p <- ggplot(dataforplotting, aes(x = as.factor(timesinceltxround), y = value,
                                     fill = cleanvar, group = cleanvar), position = "stack")
    p <- p + geom_bar(stat = "identity")
    p <- p + scale_fill_manual(values = taxacolours, name = "Genera")
    p <- p + scale_y_continuous(breaks = c(0, 0.5, 1),
                                limits = c(0, 1), labels = c(0, 0.5, 1))
    p <- p + lightertheme
    p <- p + theme(axis.text.x = element_text(angle = 270,
                                              colour = "black",
                                              size = 10,
                                              family = "Helvetica",
                                              vjust = 0.5,
                                              hjust = 0),
                   axis.title.x = element_text(size = 10,
                                               family = "Helvetica",
                                               colour = "black"),
                   axis.title.y = element_text(size = 10,
                                               family = "Helvetica",
                                               colour = "black"),
                   axis.text.y = element_text(size = 10,
                                              family = "Helvetica",
                                              colour = "black"),
                   panel.margin = unit(c(7), "mm"))
    p <- p + xlab("") + ylab("")
    p <- p + facet_grid(~patient, scales = "free", space = "free")
    if (grepl("col", ignore.case = T, diagnosisgroup)) {
        calwidth <- 7.86/2
    } else {
        calwidth <- 7.86
    }
    if (makelegend == F) {
        p <- p + guides(colour = F, fill = F)
        ## fix the coordinates ratio: y / x
        p <- p + coord_equal(4.299)
        calheight <- 2.29
    } else {
        p <- p + guides(fill = guide_legend(override.aes = list(size = 5.5),
                                            label.theme = element_text(size = 10,
                                                                       angle = 0,
                                                                       face = "italic"),
                                            title.theme = element_text(size = 10,
                                                                       angle = 0,
                                                                       face = "bold"),
                                            ncol = 1))
        ## Make graph really tall
        calheight <- 4.5 * 1.4
    }
    ggsave(sprintf("./output/lungbiome_stackedhistogram_%s.pdf", diagnosisgroup), plot = p,
           width = calwidth, height = calheight, units = "in", limitsize = F)
}

## Function for making diversity panels ##
maxsimpson <- ceiling(max(countmatrix$invsimpson))

plotdiversitypanels <- function(dataforplotting, diagnosisgroup) {
    ## For the line to join points across categorical variables,
    ## x axis needs to be numeric (or continuous)
    sapply(unique(dataforplotting$patient), function(patientid) {
        patientdata <- dataforplotting[which(dataforplotting$patient == patientid), ]
        p <- ggplot(patientdata, aes(x = as.numeric(as.factor(timesinceltx)), y = invsimpson))
        p <- p + geom_point(colour = "#2C5777")
        p <- p + geom_line(guide = F, colour = "#2C5777")
        p <- p + coord_equal(0.254)
        p <- p + scale_x_discrete(labels = sort(unique(patientdata$timesinceltx)))
        p <- p + scale_y_continuous(breaks = c(5, 10),
                                    limits = c(0, cmaxsimpson),
                                    labels = c(5, 10))
        p <- p + lightertheme
        p <- p + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.y = element_text(size = 10,
                                                  family = "Helvetica",
                                                  colour = "black"),
                       plot.margin = unit(c(1, 1, -7, 0), "mm"))
        calwidth <- length(unique(patientdata$timesinceltx)) * 1.2
        calheight <- 4.85/2.5
        patient_diagnosis <- paste(patientid, diagnosisgroup, sep = "_")
        ggsave(sprintf("./output/lungbiome_diversitypanel_%s.pdf", patient_diagnosis), plot = p,
               width = calwidth, height = calheight, units = "cm", limitsize = F)
    }, simplify = F)
}

## Generate all graphs ##
sapply(names(patientsets), function(diagnosisgroupname) {
    plottingdata <- preparedata(dataprop = propmatrix,
                                datacount = countmatrix,
                                diagnosisgroup = diagnosisgroupname,
                                finaltaxa = finaltaxa)
    plotstackedhistograms(dataforplotting = plottingdata,
                          diagnosisgroup = diagnosisgroupname,
                          makelegend = F)
    plotdiversitypanels(dataforplotting = plottingdata,
                        diagnosisgroup = diagnosisgroupname)
    }, simplify = F)

## At the end of this script, you should have the following files:
## ./output/lungbiome_diversitypanel_P{patient_number}_Pneumonia.pdf (Figure 1A, diversity panels)
## ./output/lungbiome_diversitypanel_P{patient_number}_Tracheobronchitis.pdf (Figure 1B, diversity panels)
## ./output/lungbiome_diversitypanel_P{patient_number}_Colonization.pdf (Figure 1C, diversity panels)
## ./output/lungbiome_stackedhistogram_Pneumonia.pdf (Figure 1A, Microbiome stacked histograms)
## ./output/lungbiome_stackedhistogram_Tracheobronchitis.pdf (Figure 1B, Microbiome stacked histograms)
## ./output/lungbiome_stackedhistogram_Colonization.pdf (Figure 1C, Microbiome stacked histograms)