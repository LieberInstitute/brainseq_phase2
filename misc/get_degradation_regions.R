###

library(jaffelab)
library(rtracklayer)
library(recount)
library(recount.bwtool)
library(BiocParallel)
library(SummarizedExperiment)

#######################
##### HIPPO ###########
#######################

## read manifest
load("count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda")
pdHippo = colData(rse_gene)
	
## import degradation regions
bed = import("/dcl01/lieber/ajaffe/lab/degradation_experiments/Hippo_RiboZero/bed/HIPPO_RiboZero_degradation_regions_top1000.bed")

## designate bigwigs
forwardBw = paste0("preprocessed_data/Hippo_RiboZero/Coverage/",
	sapply(pdHippo$SAMPLE_ID,"[",1),".Forward.bw")
reverseBw = paste0("preprocessed_data/Hippo_RiboZero/Coverage/",
	sapply(pdHippo$SAMPLE_ID,"[",1), ".Reverse.bw")
all(file.exists(c(forwardBw,reverseBw))) # TRUE
names(forwardBw) = names(reverseBw) = sapply(pdHippo$SAMPLE_ID,"[",1)

## try coverage tool
covForward = coverage_bwtool(forwardBw, bed, strand = "+", 
	sumsdir = "degradation", bpparam = MulticoreParam(8))
covForward$bigwig_path = NULL
covForward$bigwig_file = NULL

covReverse = coverage_bwtool(reverseBw, bed, strand = "-", 
	sumsdir = "degradation", bpparam = MulticoreParam(8))
covReverse$bigwig_path = NULL
covReverse$bigwig_file = NULL

## combine
cov_rse = rbind(covForward, covReverse)	
rownames(cov_rse) = rowData(cov_rse)$name
cov_rse = cov_rse[bed$name,]

## divide by number of reads
assays(cov_rse)$counts = assays(cov_rse)$counts/100 # divide by read length

## make positive
assays(cov_rse)$counts = abs(assays(cov_rse)$counts) 

## add which ones are bonf
bedBonf= import("/dcl01/lieber/ajaffe/lab/degradation_experiments/Hippo_RiboZero/bed/HIPPO_RiboZero_degradation_regions_bonf.bed")
rowData(cov_rse)$bonfSig = FALSE
rowData(cov_rse)$bonfSig[rownames(cov_rse) %in% bedBonf$name] = TRUE

## save to final people
cov_rse_hippo = cov_rse
save(cov_rse_hippo, file = "count_data/degradation_rse_phase2_hippo.rda")

#####################################
### DLPFC ########################
############################

## import degradation regions
bedDlpfc = import("/dcl01/lieber/ajaffe/lab/degradation_experiments/DLPFC_RiboZero/bed/DLPFC_RiboZero_degradation_regions_bonf.bed")

## designate bigwigs
forwardBw = paste0("preprocessed_data/Hippo_RiboZero/Coverage/",
	sapply(pdHippo$SAMPLE_ID,"[",1),".Forward.bw")
reverseBw = paste0("preprocessed_data/Hippo_RiboZero/Coverage/",
	sapply(pdHippo$SAMPLE_ID,"[",1), ".Reverse.bw")
all(file.exists(c(forwardBw,reverseBw))) # TRUE
names(forwardBw) = names(reverseBw) = sapply(pdHippo$SAMPLE_ID,"[",1)

## try coverage tool
covForward = coverage_bwtool(forwardBw, bedDlpfc, strand = "+", 
	sumsdir = "degradation", bpparam = MulticoreParam(8))
covForward$bigwig_path = NULL
covForward$bigwig_file = NULL

covReverse = coverage_bwtool(reverseBw, bedDlpfc, strand = "-", 
	sumsdir = "degradation", bpparam = MulticoreParam(8))
covReverse$bigwig_path = NULL
covReverse$bigwig_file = NULL

## combine
cov_rse = rbind(covForward, covReverse)	
rownames(cov_rse) = rowData(cov_rse)$name
cov_rse = cov_rse[bedDlpfc$name,]

## divide by number of reads
assays(cov_rse)$counts = assays(cov_rse)$counts/100 # divide by read length

## make positive
assays(cov_rse)$counts = abs(assays(cov_rse)$counts) 

## save to final people
cov_rse_hippo = cov_rse
save(cov_rse_hippo, file = "count_data/degradation_rse_phase2_hippo.rda")


##############################################
## HIPPO DATA w/ DLPFC DEGRADATION REGIONS ###
##############################################

## designate bigwigs from hippo
forwardBw = paste0("preprocessed_data/Hippo_RiboZero/Coverage/",
	sapply(pdHippo$SAMPLE_ID,"[",1),".Forward.bw")
reverseBw = paste0("preprocessed_data/Hippo_RiboZero/Coverage/",
	sapply(pdHippo$SAMPLE_ID,"[",1), ".Reverse.bw")
all(file.exists(c(forwardBw,reverseBw))) # TRUE
names(forwardBw) = names(reverseBw) = sapply(pdHippo$SAMPLE_ID,"[",1)

## try coverage tool w/ dlpfc regions
covForward = coverage_bwtool(forwardBw, bedDlpfc, strand = "+", 
	sumsdir = "degradation", bpparam = MulticoreParam(8))
covForward$bigwig_path = NULL
covForward$bigwig_file = NULL

covReverse = coverage_bwtool(reverseBw, bedDlpfc, strand = "-", 
	sumsdir = "degradation", bpparam = MulticoreParam(8))
covReverse$bigwig_path = NULL
covReverse$bigwig_file = NULL

## combine
cov_rse = rbind(covForward, covReverse)	
rownames(cov_rse) = rowData(cov_rse)$name
cov_rse = cov_rse[bedDlpfc$name,]

## divide by number of reads
assays(cov_rse)$counts = assays(cov_rse)$counts/100 # divide by read length

## make positive
assays(cov_rse)$counts = abs(assays(cov_rse)$counts) 

## save to final people
cov_rse_hippo_dlpfc_qSV = cov_rse
save(cov_rse_hippo_dlpfc_qSV, file = "count_data/degradation_rse_phase2_hippo_usingDlpfcDegrade.rda")


#################
#### explore ####
#################

## get gene RPKMs
load("count_data/astellas_dg_hg38_rseGene_n263.rda")
getRPKM = recount::getRPKM
geneRpkm = getRPKM(rse_gene, length="Length")
pd = colData(rse_gene)
pd$Dx = factor(pd$Dx)
pd$Dx = relevel(pd$Dx, "Control")

# do pca on genes
pca = prcomp(t(log2(geneRpkm+1)))
pcaVars = getPcaVars(pca)

## check output
plot(pca$x[,1] ~ pd$totalAssignedGene)
plot(pca$x[,2] ~ pd$RIN)
plot(pca$x[,3] ~ pd$mitoRate)
plot(pca$x[,3] ~ factor(pd$BatchLab))


## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse[rowData(cov_rse)$bonfSig,])$counts+1)))
getPcaVars(qsvBonf)[1:5]

pdf("qcChecks/qSVs_bonfRegions.pdf")
mypar(2,2)
plot(qsvBonf$x[,1] ~ pd$totalAssignedGene,
	xlab = "Gene Assignment Rate",pch=19,bg="grey",
	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,2] ~ factor(pd$BatchLab),
		xlab = "Batch",	ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"))
plot(qsvBonf$x[,3] ~ factor(pd$BatchLab),
		xlab = "Batch",	ylab=paste0("qSV3: ",getPcaVars(qsvBonf)[3],"% Var Expl"))
plot(qsvBonf$x[,4] ~ factor(pd$BatchLab),
		xlab = "Batch",	ylab=paste0("qSV4: ",getPcaVars(qsvBonf)[4],"% Var Expl"))
dev.off()

summary(lm(qsvBonf$x[,1] ~ pd$totalAssignedGene))
summary(lm(qsvBonf$x[,1] ~ pd$Dx))
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$totalAssignedGene))


## get qSVs for top 1000
qsvTop = prcomp(t(log2(assays(cov_rse)$counts+1)))
getPcaVars(qsvTop)[1:5]
plot(qsvTop$x[,1] ~ pd$totalAssignedGene)
plot(qsvTop$x[,1] ~ factor(pd$Dx))
plot(qsvTop$x[,2] ~ factor(pd$Batch))
plot(qsvTop$x[,3] ~ factor(pd$Batch))
plot(qsvTop$x[,4] ~ factor(pd$Batch))
plot(qsvTop$x[,4] ~ pd$overallMapRate)

cc = cor(qsvBonf$x[,1:10], qsvTop$x[,1:10])
signif(cc,2)