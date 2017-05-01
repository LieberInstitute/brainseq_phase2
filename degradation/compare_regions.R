###

library(jaffelab)
library(limma)

## load pheno
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/overall_degradation_pheno.rda")

## load DLPFC
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/rpkmCounts_DLPFC_riboZero_degradation_hg38_n20.rda")
pdD = pd[rownames(metrics),]
pdD = cbind(pdD, metrics)
geneRpkmD = geneRpkm

### load hippocampus
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/rpkmCounts_HIPPO_riboZero_degradation_hg38_n12.rda")
pdH = pd[rownames(metrics),]
pdH = cbind(pdH, metrics)
geneRpkmH = geneRpkm

## filter expression
m = rowMeans(cbind(geneRpkmD, geneRpkmH))
geneRpkmD = geneRpkmD[m > 0.05,]
geneRpkmH = geneRpkmH[m > 0.05,]
geneMap = geneMap[m > 0.05,]

## fit model
modD = model.matrix(~DegradationTime + factor(BrNum),data=pdD)
fitD = lmFit(log2(geneRpkmD+1), modD)
ebD = ebayes(fitD)

modH = model.matrix(~DegradationTime + factor(BrNum),data=pdH)
fitH = lmFit(log2(geneRpkmH+1), modH)
ebH = ebayes(fitH) 

## compare
plot(fitD$coef[,2], fitH$coef[,2])
plot(ebD$t[,2], ebH$t[,2])