###

library(jaffelab)
library(GenomicRanges)
library(limma)

## load results
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eqtl_tables/mergedEqtl_output_hippo_4features.rda")
hippoEqtl = allEqtl
hippoEqtl = hippoEqtl[hippoEqtl$FDR < 0.01,]

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda")
dlpfcEqtl = allEqtl
dlpfcEqtl = dlpfcEqtl[dlpfcEqtl$FDR < 0.01,]

rm(allEqtl)

##########################
## hippocampus only ######
##########################
hippoList = split(hippoEqtl, factor(hippoEqtl$Type,
	levels=c("Gene","Exon","Jxn", "Tx")))

sapply(hippoList, function(x) max(x$pvalue))

sapply(hippoList, function(x) length(unique(x$EnsemblGeneID)))
sapply(hippoList, function(x) length(unique(x$Symbol)))
sapply(hippoList, function(x) length(unique(x$snps)))

sapply(hippoList, function(x) signif(quantile(abs(x$beta), c(0.25,0.5,0.75)),3))

sapply(hippoList, function(x) table(x$Class[!duplicated(x$gene)]))$Jxn
sapply(hippoList, function(x) prop.table(table(x$Class)))$Jxn
