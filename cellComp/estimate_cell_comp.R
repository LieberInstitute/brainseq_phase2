##

library(jaffelab)
library(SummarizedExperiment)
library(recount)
library(minfi)
library(RColorBrewer)

## load RDAs
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda")
rse_gene_hippo = rse_gene
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseGene_merged_n453.rda")
rse_gene_dlpfc = rse_gene
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata")
getRPKM=recount::getRPKM # overwrite in case

## load coefficients
load("singleCell_quake_coefEsts_calibration_Zscale_hg38.rda")
table(rownames(coefEsts) %in% rownames(rse_gene))

## expression filter
yExprs_dlpfc = log2(getRPKM(rse_gene_dlpfc,"Length")+1)
yExprs_hippo = log2(getRPKM(rse_gene_hippo,"Length")+1)

yExprs_dlpfc_Scaled = scale(yExprs_dlpfc[rownames(coefEsts),])
yExprs_hippo_Scaled = scale(yExprs_hippo[rownames(coefEsts),])

dlpfc_PropEsts = minfi:::projectCellType(yExprs_dlpfc_Scaled, coefEsts)
hippo_PropEsts = minfi:::projectCellType(yExprs_hippo_Scaled, coefEsts)
propEsts = rbind(dlpfc_PropEsts, hippo_PropEsts)
propEsts = propEsts[colnames(rse_gene),]
propEsts = as.data.frame(propEsts)
save(propEsts, file = "RNA_cell_proportions_brainSeq_phase2.rda")

## reload filtered data
load("RNA_cell_proportions_brainSeq_phase2.rda")
load("../expr_cutoff/rse_gene.Rdata")

propEsts = propEsts[colnames(rse_gene),]

palette(c("darkgoldenrod2", "steelblue1")) 

## ADULT, change to age bins
plot(propEsts$Neurons ~ log2(rse_gene$Age+1), 
	pch = 21, bg = factor(rse_gene$Region))
boxplot(propEsts$Neurons ~ rse_gene$Region, 
	subset = rse_gene$Age > 17)
plot(propEsts$Astrocytes ~ log2(rse_gene$Age+1), 
	pch = 21, bg = factor(rse_gene$Region))
boxplot(propEsts$Astrocytes ~ rse_gene$Region, 
	subset = rse_gene$Age > 17)	
plot(propEsts$Oligodendrocytes ~ log2(rse_gene$Age+1), 
	pch = 21, bg = factor(rse_gene$Region))
boxplot(propEsts$Oligodendrocytes ~ rse_gene$Region, 
	subset = rse_gene$Age > 17)	
plot(propEsts$Microglia ~ log2(rse_gene$Age+1), 
	pch = 21, bg = factor(rse_gene$Region))
boxplot(propEsts$Microglia ~ rse_gene$Region, 
	subset = rse_gene$Age > 17)	
plot(propEsts$Endothelial ~ log2(rse_gene$Age+1), 
	pch = 21, bg = factor(rse_gene$Region))
boxplot(propEsts$Endothelial ~ rse_gene$Region, 
	subset = rse_gene$Age > 17)	
	
### FETAL
plot(propEsts$Fetal_quiescent ~ log2(rse_gene$Age+1), 
	pch = 21, bg = factor(rse_gene$Region))
plot(propEsts$Fetal_quiescent ~ rse_gene$Age, 
	pch = 21, bg = factor(rse_gene$Region), 
	subset = rse_gene$Age < 0)
boxplot(propEsts$Fetal_quiescent ~ rse_gene$Region, 
	subset = rse_gene$Age < 0 )	
	
plot(propEsts$Fetal_replicating ~ log2(rse_gene$Age+1), 
	pch = 21, bg = factor(rse_gene$Region))
plot(propEsts$Fetal_replicating ~ rse_gene$Age, 
	pch = 21, bg = factor(rse_gene$Region), 
	subset = rse_gene$Age < 0)
boxplot(propEsts$Fetal_replicating ~ rse_gene$Region, 
	subset = rse_gene$Age < 0 )	

	
## add p-values using lmer for the various comparisons
plot(propEsts$Neurons ~ log2(rse_gene$Age+1), 
	pch = 21, bg = factor(rse_gene$Region))