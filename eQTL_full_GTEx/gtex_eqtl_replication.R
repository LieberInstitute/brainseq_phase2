####
### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)

load("eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)
dlpfc_gtex = allEqtl
load("eqtl_tables/mergedEqtl_output_hippo_4features.rda")
hippo_gtex = allEqtl
rm(allEqtl)



























	  