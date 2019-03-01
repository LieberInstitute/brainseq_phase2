## A cleaner and more focused script that explore_twas.R

library('tibble')
library('sessioninfo')
# library('purrr')
# library('dplyr')

load('rda/twas_exp.Rdata', verbose = TRUE)

## Andrew's exploration code that focuses on the 'all' part
tt = as.data.frame(twas_exp$all)
tt = tt[which(tt$type == "psycm"),]

ttHip = tt[tt$region == "HIPPO",]
ttDLPFC = tt[tt$region == "DLPFC",]

ttDLPFC$FDR = p.adjust(ttDLPFC$TWAS.P,"fdr")
ttHip$FDR = p.adjust(ttHip$TWAS.P,"fdr")

ttSig_hip = ttHip[which(ttHip$FDR < 0.05),]
ttSig_hip = ttSig_hip[order(ttSig_hip$TWAS.P),]
ttSig_dlp = ttDLPFC[which(ttDLPFC$FDR < 0.05),]
ttSig_dlp = ttSig_dlp[order(ttSig_dlp$TWAS.P),]

dim(ttSig_dlp)
dim(ttSig_hip)

length(unique(ttSig_dlp$geneid)) #  1519
length(unique(ttSig_hip$geneid)) #  1256

## Continue


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
