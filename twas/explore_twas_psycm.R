## A cleaner and more focused script that explore_twas.R

library('tibble')
library('sessioninfo')
library('purrr')
# library('dplyr')

load('rda/twas_exp.Rdata', verbose = TRUE)

## Andrew's exploration code that focuses on the 'all' part
tt <- twas_exp$all
## Drop TWAS NA p-values
tt <- tt[!is.na(tt$TWAS.P), ]
## Focus on CLOZUK+PGC2 (psycm) GWAS
tt <- tt[which(tt$type == "psycm"),]

## Compute FDR by region for each feature 4 features
tt <- map_dfr(split(tt, tt$region), function(reg) {
    map_dfr(split(reg, reg$feature), function(reg_feat) {
        reg_feat$TWAS.FDR <- p.adjust(reg_feat$TWAS.P, 'fdr')
        reg_feat <- reg_feat[order(reg_feat$TWAS.P), ]
        return(reg_feat)
    })
})
print(tt, width = 200)

ttSig <- map(split(tt, tt$region), ~ .x[.x$TWAS.FDR < 0.05, ])
map(ttSig, dim)
# $DLPFC
# [1] 5760   27
#
# $HIPPO
# [1] 4081   27

map_int(ttSig, ~ length(unique(.x$geneid)))
# DLPFC HIPPO
 # 1514  1255

## Continue

## Add GWAS p-value and OR from the original sumstats file

## Venn diagrams of features by region, then joint (grouped by gene id)

## Compare TWAS Z-scores across DLPFC and HIPPO

## Read in the 179 CLOZUK+PGC2 snps

## Use the raggr output to find the proxy snps


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
