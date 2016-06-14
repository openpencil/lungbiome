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
## This script loads the annotations and data from the lungbiome project.
## **********************************************************************************
## 0. Define the directory with the datasets and the annotations. #####
assetpath <- "./assets/"

##### 1. Load project specific data and annotations #####
projectdata <- readRDS(file = paste(assetpath, "lungbiome_projectdata.RDS", sep = ""))

## cytokines
cytokines <- projectdata$cytokines_sn
## reflevel: CE (COPD/Emphysema, the most common indication for lung transplant)
ltxindications <- projectdata$ltxindications
## reflevel: P10
patients <- projectdata$patients
## reflevel: None
therapyvars <- projectdata$therapy
## clinical covariables, patient characteristics
clinical <- projectdata$clinical
## microbiome variables with robust reads
microbiome <- projectdata$microbes
## continuous variables
continuousvars <- projectdata$continuous
## categorical variables
categoricalvars <- projectdata$categorical
## order of cytokines for plotting
ordercyto <- projectdata$ordercyto
## order of microbes for plotting
ordermicro <- projectdata$ordermicro
## BLAST information
unclassdata <- projectdata$unclassdata
treedata <- projectdata$treedata
## colours for graphing
taxacolours <- projectdata$taxacolours
taxanames <- projectdata$taxanames
varcolours <- projectdata$varcolours
varnames <- projectdata$varnames
## red <- brewer.pal(11, 'RdBu')[3:9][1]
## blue <- brewer.pal(11, 'RdBu')[3:9][7]
## posneg <- c(red, blue)
## names(posneg) <- c('pos', 'neg')
posnegcolour <- projectdata$posnegcolour

## diagnoses and patient groups
pnapatients <- sprintf("P%s", projectdata$pnapatients)
tbrpatients <- sprintf("P%s", projectdata$tbrpatients)
colpatients <- sprintf("P%s", projectdata$colpatients)
patientsets <- list(Pneumonia = pnapatients,
                    Tracheobronchitis = tbrpatients,
                    Colonization = colpatients)

##### 2. Load lungbiome dataframes #####
## For data pre-processing steps please see
## Shankar et al., 2015
regmatrix_diagnoses <- projectdata$diagnoses_matrix
regmatrix_therapy <- projectdata$therapy_matrix
countmatrix <- projectdata$counts
propmatrix <- projectdata$props
