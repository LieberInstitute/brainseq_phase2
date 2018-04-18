## Based on the following script by Emily Burke:
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/wgcna/get_ind_exons_jxns.R

library('jaffelab')
library('SummarizedExperiment')
library('sva')
library('edgeR')
library('limma')
library('recount')
library('WGCNA')
library('parallel')
library('devtools')
allowWGCNAThreads()

cores <- 4
dir.create('rda', showWarnings = FALSE)

## load expression data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata", verbose = TRUE)

###############################################################
################## get qSVA model terms
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs.Rdata", verbose = TRUE)

min(rowMeans(assays(rse_gene)$rpkm)) 					## already cutoff to 0.25
min(rowMeans(assays(rse_exon)$rpkm)) 					## already cutoff to 0.30
min(rowMeans(assays(rse_jxn)$rp10m))  					## already cutoff to 0.46


stopifnot(identical(rse_gene$SAMPLE_ID, rse_exon$SAMPLE_ID))
stopifnot(identical(rse_gene$SAMPLE_ID, rse_jxn$SAMPLE_ID))
stopifnot(identical(colnames(rse_gene), rownames(modQsva)))

##################
## filter for age
keepIndex = which(rse_gene$Age > 17)
length(keepIndex)
rse_gene = rse_gene[,keepIndex]
rse_exon = rse_exon[,keepIndex]
rse_jxn = rse_jxn[,keepIndex]

mod <- mod[keepIndex, ]
modQsva <- modQsva[keepIndex, ]

## WGCNA options
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads= cores)

exonRpkm = assays(rse_exon)$rpkm
exonMap = rowRanges(rse_exon)
exonExprs = log2(exonRpkm+1) # transform

## normalize
exonExprsAdj = cleaningY(exonExprs, modQsva, P=3)

## drop based on correlation
exonList = split(as.data.frame(exonExprsAdj), exonMap$gencodeID)

## keep mono exonic
exonListOne = exonList[sapply(exonList,nrow) == 1]
e1 = sapply(exonListOne, rownames)

## cluster and cut exons for multiexonic genes
exonListTwo = exonList[sapply(exonList,nrow) > 1]
exonDistList = mclapply(exonListTwo, function(x) as.dist(1 - cor(t(x))), mc.cores = cores)
exonHclustList = mclapply(exonDistList, hclust, mc.cores= cores)
exonCutList = mclapply(exonHclustList, cutree, h=0.2, mc.cores= cores)
e2 = unlist(mclapply(exonCutList, function(x) names(x[!duplicated(x)]), mc.cores = cores))

# extreact
exonExprsAdj2 = exonExprsAdj[c(e1,e2),]
exonMap2 = exonMap[c(e1,e2),]

save(exonExprsAdj2, exonMap2, file="rda/independent_exons.rda")





###############################################################
############# filter features - junctions
############# (same for both regions)


## WGCNA options
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = cores)

jRpkm = assays(rse_jxn)$rp10m
jMap = rowRanges(rse_jxn)
jExprs = log2(jRpkm+1) # transform

## normalize
jExprsAdj = cleaningY(jExprs, modQsva, P=3)

## drop based on correlation
jList = split(as.data.frame(jExprsAdj), jMap$newGeneID)

## keep mono jic
jListOne = jList[sapply(jList,nrow) == 1]
j1 = sapply(jListOne, rownames)

## cluster and cut js for multijic genes
jListTwo = jList[sapply(jList,nrow) > 1]
jDistList = mclapply(jListTwo, function(x) as.dist(1 - cor(t(x))), mc.cores = cores)
jHclustList = mclapply(jDistList, hclust, mc.cores = cores)
jCutList = mclapply(jHclustList, cutree, h=0.2, mc.cores = cores)
j2 = unlist(mclapply(jCutList, function(x) names(x[!duplicated(x)]), mc.cores = cores))

# extract
jExprsAdj2 = jExprsAdj[c(j1,j2),]
jMap2 = jMap[c(j1,j2),]

save(jExprsAdj2, jMap2, file="rda/independent_jxns.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
