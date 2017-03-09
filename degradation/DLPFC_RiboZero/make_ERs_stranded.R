###

library(rtracklayer)
library(derfinder)

## load stranded mean bigwigs
meanBWs = paste0("Coverage/mean.",c("Forward","Reverse"), ".bw")
names(meanBWs) = c("Forward","Reverse")
meanList = lapply(meanBWs, import)

## make ERs
