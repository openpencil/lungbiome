##*****************************************************************************************************
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
## This script provides the code for visualizing the MLR-BMA model findings.
## **********************************************************************************************

##### 1. Load libraries and data #####
source("./lungbiome_loadprojectdata.R")
source("./lungbiome_algorithms.R")
source("./lungbiome_utilities.R")
mlr_plottingdata <- readRDS(file = "./output/lungbiome_MLR_modeldataforplotting.RDS")

###### 2. Annotations for graphing ######
varlabels <- function(vars) {
    out <- varnames[vars]
    return(out)
}

diagnosescontrastlabels <- function(var, value) {
    value <- as.character(value)
    if (var == "contrast") {
        value[value == "b20"] <- "Pneumonia Vs.\nColonization"
        value[value == "b10"] <- "Tracheobronchitis Vs.\nColonization"
        value[value == "b21"] <- "Pneumonia Vs.\nTracheobronchitis"
    } else if (var == "facetcontrast") {
        value[value == "1_PC"] <- "Pneumonia Vs.\nColonization"
        value[value == "2_TC"] <- "Tracheobronchitis Vs.\nColonization"
        value[value == "3_PT"] <- "Pneumonia Vs.\nTracheobronchitis"
    }
    return(value)
}

diaglabels <- c("1_PC", "2_TC", "3_PT")
names(diaglabels) <- c("b20", "b10", "b21")

##### 3. Plot MLR-BMA findings overview #####
plotbmaoverview <- function(mlr_plottingdata) {
    allmlrdatacombined <- ldply(mlr_plottingdata)
    allmlrdatacombined$facetcontrast <- diaglabels[allmlrdatacombined$.id]
    plottingdata <- allmlrdatacombined[which(is.na(allmlrdatacombined$welch) == F),
        ]
    maxpip <- max(plottingdata$ipcent)
    maxwelch <- max(plottingdata$welch)
    minwelch <- min(plottingdata$welch)
    p <- ggplot(data = plottingdata, aes(x = ipcent, y = welch, colour = combo))
    p <- p + geom_hline(yintercept = 0, linetype = "dotted", colour = "grey30")
    p <- p + geom_point(size = 2.5)
    p <- p + scale_colour_manual(values = varcolours, labels = varlabels)
    p <- p + scale_x_continuous(breaks = seq(0, 25, 5), limits = c(0, ceiling(maxpip)))
    p <- p + facet_grid(facetcontrast ~ ., labeller = diagnosescontrastlabels)
    p <- p + lightertheme
    p <- p + xlab("Posterior Inclusion Probability (in %)") +
             ylab("Mean difference in levels across groups")
    p <- p + theme(axis.title.x = element_text(size = 10, family = "Helvetica", colour = "black"),
                   axis.title.y = element_text(size = 10, family = "Helvetica", colour = "black"),
                   axis.text.x = element_text(size = 10, family = "Helvetica", colour = "black"),
                   axis.text.y = element_text(size = 10, family = "Helvetica", colour = "black"),
                   legend.position = "top")
    p <- p + guides(colour = guide_legend(override.aes = list(size = 5.5),
                                          label.theme = element_text(size = 10, angle = 0),
                                          title.theme = element_text(size = 10, angle = 0, face = "bold"),
                                          ncol = 2, byrow = T, title = "Covariates",
                                          title.position = "top"))
    ggsave("./output/lungbiome_MLR_BMA_overview_plot.pdf", plot = p,
           width = 6.6, height = 7.7, units = "in", limitsize = F)
}

## Figure 2 in the manuscript. ##
plotbmaoverview(mlr_plottingdata = mlr_plottingdata)

##### 4. Plot microbiome signatures #####
plotmicrobiomesignature <- function(plottingdata, plainorpolar, ordermicro,
                                    pipcutoff, contrastname) {
    plottingdata <- data.table(plottingdata)
    ## get the sign to apply colour coding
    plottingdata[, `:=`(signwelchplot, sign(x = as.numeric(welch))), ]
    plottingdata[, `:=`(posneg, ifelse(signwelchplot == "-1", "neg", "pos")), ]
    plottingdata[, `:=`(signbetax, sign(betax)), ]
    ## limit the variables to all microbes
    datamerge <- merge(plottingdata, ordermicro, by = "variable")
    datamerge[, `:=`(welch, ifelse(is.na(welch), 0, welch)), ]
    datamerge[, `:=`(newvar, paste(sprintf("%03d", ordernum), variable, sep = "_")),
        ]
    ## limit to increased odds of a particular response (positive odds)
    megdata <- subset(x = datamerge, signbetax == 1 | is.na(signbetax))
    megdata$welchrescale <- keepsignrescale(megdata$welch)
    megdata$cleanmicrobe <- cleanup(megdata$variable)
    if (plainorpolar == "polar") {
        ## Create a duplicate of datamerge and bind to datamerge to reach a semicircle
        megcopy <- megdata
        megcopy$newvar <- paste(sprintf("%03d", megcopy$ordernum + max(megcopy$ordernum)),
                                megcopy$variable, sep = "_")
        megcopy$welchrescale <- 0
        ## Merge the two
        megdata <- rbind(megdata, megcopy)
    }
    ## Filtering by PIP
    megsub <- megdata[ipcent > pipcutoff]
    inlabel <- megsub$finalmicrobe
    names(inlabel) <- megsub$newvar
    renametaxa <- function(input) {
        out <- inlabel[input]
        return(out)
    }
    ## retain only the significant
    megsub <- megsub[sig == "sig"]
    p <- ggplot(megsub, aes(y = welchrescale, x = newvar, fill = posneg))
    p <- p + geom_bar(stat = "identity")
    p <- p + scale_fill_manual(values = posnegcolour)
    p <- p + ylab("") + xlab("")
    p <- p + lightertheme
    p <- p + guides(colour = F, fill = F)
    if (plainorpolar == "polar") {
        p <- p + coord_polar(theta = "x", start = 4.67) + ylim(c(-2, 1))
        p <- p + theme(axis.title = element_blank(),
                       axis.ticks = element_blank(),
                       text = element_blank(),
                       panel.background = element_blank(),
                       panel.border = element_blank(),
                       plot.margin = unit(c(-3, 0, -5, 0), "mm"),
                       panel.grid.major = element_line(colour = "grey70", size = 0.08),
                       panel.grid.minor = element_line(colour = "grey70", size = 0.08))
        calheight <- (length(megdata$newvar) + 3 + 0.5) * 0.11
    } else {
        p <- p + theme(axis.text.y = element_text(size = 10, colour = "black"),
                       axis.title.y = element_text(size = 10, colour = "black"),
                       axis.text.x = element_text(size = 10, colour = "black",
                                                  angle = 270, hjust = 0, vjust = 0.5),
                       strip.text.x = element_blank(),
                       axis.title.x = element_blank(),
                       panel.margin = unit(0.5, "lines"))
        p <- p + scale_x_discrete(labels = renametaxa)
    }
    ggsave(sprintf("./output/lungbiome_%s_microbiome_signature_%splot.pdf", contrastname, plainorpolar),
           plot = p, height = calheight, units = "cm", limitsize = F)
}

## Figure 3A and 3B in the manuscript ##
plotmicrobiomesignature(plottingdata = mlr_plottingdata$b20,
                        plainorpolar = "polar",
                        ordermicro = ordermicro,
                        pipcutoff = 1,
                        contrastname = "pneumonia_col")
plotmicrobiomesignature(plottingdata = mlr_plottingdata$b10,
                        plainorpolar = "polar",
                        ordermicro = ordermicro,
                        pipcutoff = 1,
                        contrastname = "tracheobronchitis_col")
## Ignore this message Warning message: In loop_apply(n, do.ply) : Stacking not
## well defined when ymin != 0

##### 4. Plot cytokine signatures #####
## Limits for the barplot
alldt <- data.table(ldply(mlr_plottingdata))
maxcytokine <- max(alldt[group == "Cytokine", welch, ])
mincytokine <- min(alldt[group == "Cytokine", welch, ])

plotcytobarplot <- function(plottingdata, ordercyto, contrastname) {
    plottingdata <- data.table(plottingdata)
    ## get the sign to apply colour coding
    plottingdata[, `:=`(signwelchplot, sign(x = as.numeric(welch))), ]
    plottingdata[, `:=`(posneg, ifelse(signwelchplot == "-1", "neg", "pos")), ]
    plottingdata[, `:=`(signbetax, sign(betax)), ]
    ## limit the variables to all microbes
    datamerge <- merge(plottingdata, ordercyto, by = "variable", all.y = T)
    datamerge[, `:=`(welch, ifelse(is.na(welch), 0, welch)), ]
    datamerge[, `:=`(newvar, paste(sprintf("%02d", ordernum), variable, sep = "_")),
        ]
    ## limit to increased odds of a particular response (positive odds)
    megdata <- subset(x = datamerge, signbetax == 1 | is.na(signbetax))
    megdata$welchrescale <- keepsignrescale(megdata$welch)
    calwidth <- (length(megdata$newvar) + 3 + 0.5) * 0.4
    ## retain only significant ones
    megsub <- megdata[sig == "sig"]
    p <- ggplot(megsub, aes(y = welch, x = newvar, fill = posneg))
    p <- p + geom_bar(stat = "identity")
    p <- p + scale_fill_manual(values = posnegcolour)
    p <- p + ylab("Mean difference in\nlevels across groups") + xlab("")
    p <- p + lightertheme
    p <- p + guides(colour = F, fill = F)
    p <- p + theme(axis.text.y = element_text(size = 10, colour = "black"),
                   axis.title.y = element_text(size = 10, colour = "black"),
                   axis.text.x = element_text(size = 10, colour = "black",
                                              angle = 270, hjust = 0, vjust = 0.5),
                   strip.text.x = element_blank(), axis.title.x = element_blank(),
                   panel.margin = unit(0.5, "lines"))
    p <- p + ylim(c(mincytokine, maxcytokine))
    ## Coord_equal adjusts the height of the graph according to
    ## the calculated width above
    p <- p + coord_equal(1)
    ## Note: Greek expressions will only work if run as a function.
    ## ggsave inherits the Greek labels.
    p <- p + scale_x_discrete(labels = c(`19_FGF_basic_sn` = "FGF basic",
                                         `01_Eotaxin_sn` = "Eotaxin",
                                         `17_G_CSF_sn` = "G-CSF",
                                         `08_GM_CSF_sn` = "GM-CSF",
                                         `21_IFN_g_sn` = expression(IFN~gamma),
                                         `10_IL_1b_sn` = expression(IL-1~beta),
                                         `13_IL_1ra_sn` = "IL-1ra",
                                         `22_IL_2_sn` = "IL-2",
                                         `14_IL_4_sn` = "IL-4",
                                         `23_IL_5_sn` = "IL-5",
                                         `11_IL_6_sn` = "IL-6",
                                         `24_IL_7_sn` = "IL-7",
                                         `06_IL_8_sn` = "IL-8",
                                         `25_IL_9_sn` = "IL-9",
                                         `15_IL_10_sn` = "IL-10",
                                         `26_IL_12_p70_sn` = "IL-12 (p70)",
                                         `16_IL_13_sn` = "IL-13",
                                         `27_IL_15_sn` = "IL-15",
                                         `12_IL_17A_sn` = "IL-17A",
                                         `07_IP_10_sn` = "IP-10",
                                         `02_MCP_1_MCAF_sn` = "MCP-1 (MCAF)",
                                         `03_MIP_1a_sn` = expression(MIP-1~alpha),
                                         `04_MIP_1b_sn` = expression(MIP-1~beta),
                                         `20_PDGF_bb_sn` = "PDGF-BB",
                                         `05_RANTES_sn` = "RANTES",
                                         `09_TNF_a_sn` = expression(TNF~alpha),
                                         `18_VEGF_sn` = "VEGF"))
    ggsave(sprintf("./output/lungbiome_%s_cytokine_signature.pdf", contrastname),
           plot = p, width = calwidth, units = "cm", limitsize = F)
}

## Figure 4C and 4D in the manuscript ##
plotcytobarplot(plottingdata = mlr_plottingdata$b20, ordercyto = ordercyto, contrastname = "pneumonia_col")
plotcytobarplot(plottingdata = mlr_plottingdata$b10, ordercyto = ordercyto, contrastname = "tracheobronchitis_col")
## This warning is harmless Warning messages: 1: Stacking not well defined when
## ymin != 0 2: In is.na(labels) : is.na() applied to non-(list or vector) of type
## 'expression'

## At the end of this script, you should have the following files:
## ./output/lungbiome_MLR_BMA_overview_plot.pdf (Figure 2)
## ./output/lungbiome_pneumonia_col_microbiome_signature_polarplot.pdf (Figure 3A)
## ./output/lungbiome_pneumonia_col_microbiome_signature_polarplot.pdf (Figure 3B)
## ./output/lungbiome_pneumonia_col_cytokine_signature.pdf (Figure 3C)
## ./output/lungbiome_tracheobronchitis_col_cytokine_signature.pdf (Figure 3D)
