##

library(jaffelab)
library(SummarizedExperiment)
library(recount)
library(minfi)

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