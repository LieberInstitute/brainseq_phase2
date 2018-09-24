######

library(rtracklayer)
library(derfinder)
library(jaffelab)
library(SummarizedExperiment)
library(recount.bwtool)

## chromosome info
chrInfo <- read.table('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode',
    header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
chrInfo <- subset(chrInfo, chr %in% paste0('chr', c(1:22, 'X', 'Y','M')))

## load stranded mean bigwigs
degPaths = rep(paste0("/dcl01/lieber/ajaffe/lab/degradation_experiments/",
	c("Hippo", "DLPFC"), "_RiboZero/"),each=2)
meanBWs = paste0(degPaths, "Coverage/mean.",c("Forward","Reverse"), ".bw")
names(meanBWs) = paste0(rep(c("Hippo", "DLPFC"), each=2), "_", c("Forward","Reverse"))
all(file.exists(meanBWs))

## read in
meanList_byRegion = mclapply(meanBWs, import,mc.cores=4)

## take mean across regions
meanList = mclapply(list(Forward = c(1,3), Reverse = c(2,4)), function(ii) {
	first = coverage(meanList_byRegion[[ii[1]]])
	second = coverage(meanList_byRegion[[ii[2]]])
}, mc.cores=2)
strand(meanList$Forward) = Rle("+")
strand(meanList$Reverse) = Rle("-")

## make ERs
meanListFilter = GRangesList(lapply(meanList, function(x) x[abs(x$score) >= 5]))
reduceList = endoapply(meanListFilter, reduce,min.gapwidth=2)

## at least 12 BP and on main chromosomes
erList = endoapply(reduceList, function(x) 
	x[width(x) >= 12 & seqnames(x) %in% chrInfo$chr])

## write out
dir.create("bed")
export(erList$Forward, con = "bed/HIPPO_RiboZero_ERs_cut5_Forward.bed")
export(erList$Reverse, con = "bed/HIPPO_RiboZero_ERs_cut5_Reverse.bed")
