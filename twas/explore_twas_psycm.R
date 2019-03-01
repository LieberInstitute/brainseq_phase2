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

## Compute FDR by region across all 4 features
tt <- map_dfr(split(tt, tt$region), function(reg) {
    reg$TWAS.FDR <- p.adjust(reg$TWAS.P, 'fdr')
    reg <- reg[order(reg$TWAS.P), ]
    return(reg)
})
print(tt, width = 200)

ttSig <- map(split(tt, tt$region), ~ .x[.x$TWAS.FDR < 0.05, ])
map(ttSig, dim)
# $DLPFC
# [1] 5762   27
#
# $HIPPO
# [1] 4090   27

map_int(ttSig, ~ length(unique(.x$geneid)))
# DLPFC HIPPO
#  1519  1256

## Continue


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
