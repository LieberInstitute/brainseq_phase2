library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library('devtools')
library('readxl')
library('sva')
library('gplots')



#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')

#load qsvBonf, qSVs, mod, modQsva object from brainseq data
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs.Rdata", verbose = TRUE)

filt_amy <- rse_gene$Age>17 & rse_gene$Region == "HIPPO"
filt_kit <- rse_gene$Age>17 & rse_gene$Region == "HIPPO" & as.integer(gsub('R', '', colnames(rse_gene))) %in% hipxl$RNum[hipxl$Protocol == 'RiboZeroHMR']

addmargins(table(filt_amy, filt_kit))

#filter for age and hippocampus; no age in rse_gene
keepIndex = which(filt_kit)
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_hippo_qsva_onlyHMR.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('hippo_voom_noqsva_onlyHMR.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('hippo_voom_noadj_onlyHMR.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_hippo_filtered_qSVA_geneLevel_onlyHMR.rda")






### Re-make qSVs without prenatal samples
rm(list = ls())

#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)

## Drop prenatal
keepIndex <- which(rse_gene$Age > 0)
rse_gene <- rse_gene[, keepIndex]

## model
mod = model.matrix(~Dx + Age + Sex + mitoRate + Region + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(rse_gene))

#load cov_rse object from brainseq data 
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)

cov_rse <- cov_rse[, keepIndex]

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
#  [1] 67.800  5.580  3.090  2.840  1.370  1.220  1.020  0.956  0.823  0.753
# [11]  0.680  0.540  0.515  0.450  0.342  0.326  0.322  0.311  0.289  0.258
# [21]  0.246
modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, file = 'brainseq_phase2_qsvs_noPrenatal.Rdata')
pdf('qsvs_var_explained_noPrenatal.pdf')
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()


### Re-make qSVs with age > 17
rm(list = ls())

#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)

## Drop age<=17
keepIndex <- which(rse_gene$Age > 17)
rse_gene <- rse_gene[, keepIndex]

## model
mod = model.matrix(~Dx + Age + Sex + mitoRate + Region + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(rse_gene))

#load cov_rse object from brainseq data 
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)

cov_rse <- cov_rse[, keepIndex]

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
#  [1] 68.500  5.590  3.110  2.910  1.390  1.260  0.966  0.835  0.778  0.716
# [11]  0.543  0.508  0.442  0.361  0.333  0.328  0.313  0.283  0.258  0.244
# [21]  0.236  0.207
modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, file = 'brainseq_phase2_qsvs_age17.Rdata')
pdf('qsvs_var_explained_age17.pdf')
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()



### Re-make qSVs without age>17 and no HIPPO gold samples
rm(list = ls())

#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)

hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')

## Drop age<=17 and HIPPO Gold samples
hipHMR <- as.integer(gsub('R', '', colnames(rse_gene))) %in% hipxl$RNum[hipxl$Protocol == 'RiboZeroHMR']
hipHMR[rse_gene$Region == 'DLPFC'] <- TRUE

keepIndex <- which(rse_gene$Age > 17 & hipHMR)
rse_gene <- rse_gene[, keepIndex]

## model
mod = model.matrix(~Dx + Age + Sex + mitoRate + Region + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(rse_gene))

#load cov_rse object from brainseq data 
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)

cov_rse <- cov_rse[, keepIndex]

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
#  [1] 70.200  5.540  3.070  2.600  1.180  1.070  0.954  0.770  0.646  0.577
# [11]  0.512  0.422  0.400  0.334  0.316  0.315  0.278  0.276  0.243  0.235
# [21]  0.208  0.204
modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, keepIndex, file = 'brainseq_phase2_qsvs_age17_noHGold.Rdata')
pdf('qsvs_var_explained_age17_noHGold.pdf')
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()





### Re-make qSVs without age>17 and no HIPPO gold samples
## Only HIPPO
rm(list = ls())

#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)

hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')

## Drop age<=17 and HIPPO Gold samples
hipHMR <- as.integer(gsub('R', '', colnames(rse_gene))) %in% hipxl$RNum[hipxl$Protocol == 'RiboZeroHMR']
hipHMR[rse_gene$Region == 'DLPFC'] <- TRUE

keepIndex <- which(rse_gene$Age > 17 & hipHMR & rse_gene$Region == 'HIPPO')
rse_gene <- rse_gene[, keepIndex]

## model
mod = model.matrix(~Dx + Age + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(rse_gene))

#load cov_rse object from brainseq data 
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)

cov_rse <- cov_rse[, keepIndex]

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
#  [1] 60.200  7.850  4.630  2.460  2.040  1.640  1.240  1.060  0.974  0.770
# [11]  0.665  0.611  0.580  0.489  0.458  0.432
modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, keepIndex, file = 'brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata')
pdf('qsvs_var_explained_age17_noHGold_HIPPO.pdf')
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()




### Re-make qSVs without age>17 and no HIPPO gold samples
## Only DLPFC
rm(list = ls())

#load rse_gene object from brainseq data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)

hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')

## Drop age<=17 and HIPPO Gold samples
hipHMR <- as.integer(gsub('R', '', colnames(rse_gene))) %in% hipxl$RNum[hipxl$Protocol == 'RiboZeroHMR']
hipHMR[rse_gene$Region == 'DLPFC'] <- TRUE

keepIndex <- which(rse_gene$Age > 17 & hipHMR & rse_gene$Region == 'DLPFC')
rse_gene <- rse_gene[, keepIndex]

## model
mod = model.matrix(~Dx + Age + Sex + mitoRate + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5,
	data = colData(rse_gene))

#load cov_rse object from brainseq data 
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)

cov_rse <- cov_rse[, keepIndex]

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
#  [1] 69.600  5.620  2.730  1.900  1.420  1.130  0.875  0.869  0.845  0.633
# [11]  0.576  0.498  0.469  0.424  0.401
modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, keepIndex, file = 'brainseq_phase2_qsvs_age17_noHGold_DLPFC.Rdata')
pdf('qsvs_var_explained_age17_noHGold_DLPFC.pdf')
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()






### Run with Gold samples using the qSVs made without the pre-natal samples
rm(list = ls())
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl_HIPPO_checks/brainseq_phase2_qsvs_noPrenatal.Rdata', verbose = TRUE)

## Drop prenatal
keepIndex <- which(rse_gene$Age > 0)
rse_gene <- rse_gene[, keepIndex]

keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "HIPPO")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_hippo_qsva_noPrenatalQSV.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('hippo_voom_noqsva_noPrenatalQSV.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('hippo_voom_noadj_noPrenatalQSV.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_hippo_filtered_qSVA_geneLevel_noPrenatalQSV.rda")



### Run with Gold samples using the qSVs made without the pre-natal samples
## And adjust for kit
rm(list = ls())
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl_HIPPO_checks/brainseq_phase2_qsvs_noPrenatal.Rdata', verbose = TRUE)
hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')

## Drop prenatal
keepIndex <- which(rse_gene$Age > 0)
rse_gene <- rse_gene[, keepIndex]

keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "HIPPO")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

kit <- hipxl$Protocol[match(as.integer(gsub('R', '', colnames(rse_gene))), hipxl$RNum)]

mod <- cbind(mod, 'kitRiboZeroHMR' = model.matrix(~ kit)[, 2])
modQsva <- cbind(modQsva, 'kitRiboZeroHMR' = model.matrix(~ kit)[, 2])


##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_hippo_qsva_noPrenatalQSV_adjKit.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('hippo_voom_noqsva_noPrenatalQSV_adjKit.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('hippo_voom_noadj_noPrenatalQSV_adjKit.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_hippo_filtered_qSVA_geneLevel_noPrenatalQSV_adjKit.rda")





### Run without Gold samples using the qSVs made without the pre-natal samples
rm(list = ls())
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl_HIPPO_checks/brainseq_phase2_qsvs_noPrenatal.Rdata', verbose = TRUE)

## Drop prenatal
keepIndex <- which(rse_gene$Age > 0)
rse_gene <- rse_gene[, keepIndex]

hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')

filt_amy <- rse_gene$Age>17 & rse_gene$Region == "HIPPO"
filt_kit <- rse_gene$Age>17 & rse_gene$Region == "HIPPO" & as.integer(gsub('R', '', colnames(rse_gene))) %in% hipxl$RNum[hipxl$Protocol == 'RiboZeroHMR']

addmargins(table(filt_amy, filt_kit))

#filter for age and hippocampus; no age in rse_gene
keepIndex = which(filt_kit)
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_hippo_qsva_noPrenatalQSV_onlyHMR.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('hippo_voom_noqsva_noPrenatalQSV_onlyHMR.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('hippo_voom_noadj_noPrenatalQSV_onlyHMR.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_hippo_filtered_qSVA_geneLevel_noPrenatalQSV_onlyHMR.rda")





### Run with Gold samples using the qSVs made without the age>17 samples
rm(list = ls())
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl_HIPPO_checks/brainseq_phase2_qsvs_age17.Rdata', verbose = TRUE)

## Drop prenatal
keepIndex <- which(rse_gene$Age > 17)
rse_gene <- rse_gene[, keepIndex]

keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "HIPPO")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_hippo_qsva_age17QSV.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('hippo_voom_noqsva_age17QSV.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('hippo_voom_noadj_age17QSV.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_hippo_filtered_qSVA_geneLevel_age17QSV.rda")



### Run with Gold samples using the qSVs made without the age>17 samples
## And adjust for kit
rm(list = ls())
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl_HIPPO_checks/brainseq_phase2_qsvs_age17.Rdata', verbose = TRUE)

hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')

## Drop prenatal
keepIndex <- which(rse_gene$Age > 17)
rse_gene <- rse_gene[, keepIndex]

keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "HIPPO")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

kit <- hipxl$Protocol[match(as.integer(gsub('R', '', colnames(rse_gene))), hipxl$RNum)]

mod <- cbind(mod, 'kitRiboZeroHMR' = model.matrix(~ kit)[, 2])
modQsva <- cbind(modQsva, 'kitRiboZeroHMR' = model.matrix(~ kit)[, 2])

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_hippo_qsva_age17QSV_adjKit.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('hippo_voom_noqsva_age17QSV_adjKit.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('hippo_voom_noadj_age17QSV_adjKit.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_hippo_filtered_qSVA_geneLevel_age17QSV_adjKit.rda")




### Run without Gold samples using the qSVs made without the age>17 samples
rm(list = ls())
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl_HIPPO_checks/brainseq_phase2_qsvs_age17.Rdata', verbose = TRUE)

## Drop prenatal
keepIndex <- which(rse_gene$Age > 17)
rse_gene <- rse_gene[, keepIndex]

hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')

filt_amy <- rse_gene$Age>17 & rse_gene$Region == "HIPPO"
filt_kit <- rse_gene$Age>17 & rse_gene$Region == "HIPPO" & as.integer(gsub('R', '', colnames(rse_gene))) %in% hipxl$RNum[hipxl$Protocol == 'RiboZeroHMR']

addmargins(table(filt_amy, filt_kit))

#filter for age and hippocampus; no age in rse_gene
keepIndex = which(filt_kit)
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_hippo_qsva_age17QSV_onlyHMR.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('hippo_voom_noqsva_age17QSV_onlyHMR.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('hippo_voom_noadj_age17QSV_onlyHMR.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_hippo_filtered_qSVA_geneLevel_age17QSV_onlyHMR.rda")





### Run with Gold samples using the qSVs made without the age>17 samples
## and without the HIPPO Gold samples
rm(list = ls())
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl_HIPPO_checks/brainseq_phase2_qsvs_age17_noHGold.Rdata', verbose = TRUE)

## Drop prenatal
rse_gene <- rse_gene[, keepIndex]

keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "HIPPO")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_hippo_qsva_noHGoldQSV.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('hippo_voom_noqsva_noHGoldQSV.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('hippo_voom_noadj_noHGoldQSV.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV.rda")




### Run with Gold samples using the qSVs made without the age>17 samples
## and without the HIPPO Gold samples
rm(list = ls())
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl_HIPPO_checks/brainseq_phase2_qsvs_age17_noHGold.Rdata', verbose = TRUE)

## Drop prenatal
rse_gene <- rse_gene[, keepIndex]

keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "DLPFC")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_dlpfc_qsva_noHGoldQSV.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('dlpfc_voom_noqsva_noHGoldQSV.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('dlpfc_voom_noadj_noHGoldQSV.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV.rda")





### Run with Gold samples using the qSVs made without the age>17 samples
## and without the HIPPO Gold samples (qsv are HIPPO specific)
rm(list = ls())
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl_HIPPO_checks/brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata', verbose = TRUE)

## Drop prenatal
rse_gene <- rse_gene[, keepIndex]

keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "HIPPO")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, ]
modQsva <- modQsva[keepIndex, ]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_hippo_qsva_noHGoldQSV_HIPPO.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('hippo_voom_noqsva_noHGoldQSV_HIPPO.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('hippo_voom_noadj_noHGoldQSV_HIPPO.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV_HIPPO.rda")




### Run with Gold samples using the qSVs made without the age>17 samples
## and without the HIPPO Gold samples (qsv are DLPFC specific)
rm(list = ls())
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl_HIPPO_checks/brainseq_phase2_qsvs_age17_noHGold_DLPFC.Rdata', verbose = TRUE)

## Drop prenatal
rse_gene <- rse_gene[, keepIndex]

keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "DLPFC")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, ]
modQsva <- modQsva[keepIndex, ]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_dlpfc_qsva_noHGoldQSV_DLPFC.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('dlpfc_voom_noqsva_noHGoldQSV_DLPFC.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('dlpfc_voom_noadj_noHGoldQSV_DLPFC.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV_DLPFC.rda")







### Run with Gold samples using the qSVs made without the pre-natal samples
rm(list = ls())
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl_HIPPO_checks/brainseq_phase2_qsvs_noPrenatal.Rdata', verbose = TRUE)

## Drop prenatal
keepIndex <- which(rse_gene$Age > 0)
rse_gene <- rse_gene[, keepIndex]

keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "DLPFC")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_dlpfc_qsva_noPrenatalQSV.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('dlpfc_voom_noqsva_noPrenatalQSV.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('dlpfc_voom_noadj_noPrenatalQSV.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_dlpfc_filtered_qSVA_geneLevel_noPrenatalQSV.rda")



### Run with Gold samples using the qSVs made without the age>17 samples
rm(list = ls())
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl_HIPPO_checks/brainseq_phase2_qsvs_age17.Rdata', verbose = TRUE)

## Drop prenatal
keepIndex <- which(rse_gene$Age > 17)
rse_gene <- rse_gene[, keepIndex]

keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "DLPFC")
rse_gene <- rse_gene[, keepIndex]
mod <- mod[keepIndex, -which(colnames(mod) == 'RegionHIPPO')]
modQsva <- modQsva[keepIndex, -which(colnames(modQsva) == 'RegionHIPPO')]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('voom_dlpfc_qsva_age17QSV.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('dlpfc_voom_noqsva_age17QSV.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('dlpfc_voom_noadj_age17QSV.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Dx)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]


stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "dxStats_dlpfc_filtered_qSVA_geneLevel_age17QSV.rda")




### Check all
rm(list = ls())
files <- c( '/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl/dxStats_DLPFC_filtered_qSVA_geneLevel.rda',
    '/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl/dxStats_hippo_filtered_qSVA_geneLevel.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_geneLevel.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_geneLevel.rda',
    
    'dxStats_hippo_filtered_qSVA_geneLevel_onlyHMR.rda',
    'dxStats_hippo_filtered_qSVA_geneLevel_noPrenatalQSV.rda',
    'dxStats_hippo_filtered_qSVA_geneLevel_noPrenatalQSV_onlyHMR.rda',
    'dxStats_hippo_filtered_qSVA_geneLevel_age17QSV.rda',
    'dxStats_hippo_filtered_qSVA_geneLevel_age17QSV_onlyHMR.rda',
    'dxStats_dlpfc_filtered_qSVA_geneLevel_noPrenatalQSV.rda',
    'dxStats_dlpfc_filtered_qSVA_geneLevel_age17QSV.rda',
    'dxStats_hippo_filtered_qSVA_geneLevel_noPrenatalQSV_adjKit.rda',
    'dxStats_hippo_filtered_qSVA_geneLevel_age17QSV_adjKit.rda',
    'dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV.rda',
    'dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV.rda',
    'dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV_HIPPO.rda',
    'dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV_DLPFC.rda'

)
outGene <- lapply(files, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(outGene)
})
names(outGene) <- c('DLPFC_andrew', 'HIPPO_andrew', 'DLPFC_amy', 'HIPPO_amy', 'HIPPO_onlyHMR', 'HIPPO_noPrenatalQSV', 'HIPPO_noPrenatalQSV_onlyHMR', 'HIPPO_age17QSV', 'HIPPO_age17QSV_onlyHMR', 'DLPFC_noPrenatalQSV', 'DLPFC_age17QSV', 'HIPPO_noPrenatalQSV_adjKit', 'HIPPO_age17QSV_adjKit', 'HIPPO_noHGoldQSV', 'DLPFC_noHGoldQSV', 'HIPPO_matchQSV', 'DLPFC_matchQSV')


load('/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/caseControl/rdas/expressed_de_features.rda', verbose = TRUE)

prev <- list(
    'DLPFC_P1' = data.frame(
        ensemblID = names(outStatsExprs$Gene),
        adj.P.Val = outStatsExprs$Gene$fdr_qsva
        ),
    'CMC' = data.frame(
        ensemblID = names(outStatsExprs$Gene),
        adj.P.Val = p.adjust(outStatsExprs$Gene$CMC_pval_qsva, method = 'fdr')
        )
)

outGene <- c(outGene, prev)


options(width = 120)
res <- do.call(rbind, lapply(c(0.05, 0.1, 0.15, 0.2), function(cut) {
    xx <- sapply(outGene, function(x) {
        table(factor(x$adj.P.Val < cut, levels = c('FALSE', 'TRUE')))
    })
    cbind(xx, cutoff = cut)
}))
res

data.frame(colnames(res))

make_venn <- function(i, txt = 'HIPPO_', cut = 0.1) {
    vinfo <- lapply(outGene[i], function(x) {
        x$ensemblID[x$adj.P.Val < cut]
    })
    names(vinfo) <- gsub(txt, '', names(vinfo))
    venn(vinfo) + title(paste('FDR cutoff:', cut))
}

# rmote_off()

pdf('venn_checks.pdf', useDingbats = FALSE)
make_venn(c(2, 8, 9))
make_venn(c(2, 8, 9), cut = 0.05)

make_venn(c(2, 8, 9, 13, 14), txt = 'HIPPO_|QSV|adj')
make_venn(c(2, 8, 9, 13, 14), txt = 'HIPPO_|QSV|adj', cut = 0.05)

make_venn(c(8, 9, 13, 14), txt = 'HIPPO_|QSV|adj')
make_venn(c(8, 9, 13, 14), txt = 'HIPPO_|QSV|adj', cut = 0.05)

make_venn(c(8, 9, 11), txt = 'QSV|IPPO|LPFC')
make_venn(c(8, 9, 11), txt = 'QSV|IPPO|LPFC', cut = 0.05)

make_venn(c(14, 15, 18, 19), txt = 'QSV|IPPO|LPFC')
make_venn(c(14, 15, 18, 19), txt = 'QSV|IPPO|LPFC', cut = 0.05)

make_venn(c(9, 11, 18, 19), txt = 'QSV|IPPO|LPFC')
make_venn(c(9, 11, 18, 19), txt = 'QSV|IPPO|LPFC', cut = 0.05)

make_venn(c(9, 11, 14, 15), txt = 'QSV|IPPO|LPFC')
make_venn(c(9, 11, 14, 15), txt = 'QSV|IPPO|LPFC', cut = 0.05)

make_venn(c(11, 15, 18, 19), txt = 'QSV|IPPO|LPFC')
make_venn(c(11, 15, 18, 19), txt = 'QSV|IPPO|LPFC', cut = 0.05)


make_venn(c(2, 9, 14, 16), txt = 'HIPPO_|QSV|adj')
make_venn(c(2, 9, 14, 16), txt = 'HIPPO_|QSV|adj', cut = 0.05)


make_venn(c(15, 17, 18, 19), txt = 'QSV|IPPO|LPFC')
make_venn(c(15, 17, 18, 19), txt = 'QSV|IPPO|LPFC', cut = 0.05)

make_venn(c(16, 17, 18, 19), txt = 'QSV|IPPO|LPFC')
make_venn(c(16, 17, 18, 19), txt = 'QSV|IPPO|LPFC', cut = 0.05)


make_venn(c(14, 15, 16, 17), txt = 'QSV|IPPO|LPFC')
make_venn(c(14, 15, 16, 17), txt = 'QSV|IPPO|LPFC', cut = 0.05)


dev.off()

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

plot_dequal <- function(out_input, degrade_input, var = 't', xlabtxt = 'case-control', main, ylabtxt = '') {
	both <- intersect(rownames(out_input), rownames(degrade_input))
	degrade <- degrade_input[both, ]
	interest <- out_input[both, ]
	
	stopifnot(identical(rownames(degrade), rownames(interest)))
	corr = signif(cor( degrade[, var],  interest[, var]), 3)
	plot(y = degrade[, var], x = interest[, var], xlab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), xlabtxt), ylab = paste(ylabtxt, var, 'degradation'), main = main, col = add.alpha('black', 1/10), pch = 16)
	legend('topleft', legend = paste('r =', corr))
}

load("/dcl01/ajaffe/data/lab/qsva_brain/ERs/rdas/DLPFC_Plus_HIPPO_RiboZero_geneLevel_degradationStats_forDEqual_hg38.rda", verbose = TRUE)
load("/dcl01/ajaffe/data/lab/qsva_brain/ERs/rdas/DLPFC_HIPPO_degradationStats_hg38.rda", verbose = TRUE)



plot_dequal2 <- function(i) {
    par(mfcol = c(3, 2))
    v = 't'
    plot_dequal(outGene[[i]], degradeStats, var = v, xlabtxt = 'case-control (with qSVs)', main = names(outGene)[i], ylabtxt = 'Combined')
    plot_dequal(outGene[[i]], degradeStats_HIPPO, var = v, xlabtxt = 'case-control (with qSVs)', main = names(outGene)[i], ylabtxt = 'HIPPO')
    plot_dequal(outGene[[i]], degradeStats_DLPFC, var = v, xlabtxt = 'case-control (with qSVs)', main = names(outGene)[i], ylabtxt = 'DLPFC')
    
    v = 'logFC'
    plot_dequal(outGene[[i]], degradeStats, var = v, xlabtxt = 'case-control (with qSVs)', main = names(outGene)[i], ylabtxt = 'Combined')
    plot_dequal(outGene[[i]], degradeStats_HIPPO, var = v, xlabtxt = 'case-control (with qSVs)', main = names(outGene)[i], ylabtxt = 'HIPPO')
    plot_dequal(outGene[[i]], degradeStats_DLPFC, var = v, xlabtxt = 'case-control (with qSVs)', main = names(outGene)[i], ylabtxt = 'DLPFC')
}


pdf('dequal_checks.pdf', useDingbats = FALSE, width = 10, height = 15)
for(i in 1:17) plot_dequal2(i)
dev.off()


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

