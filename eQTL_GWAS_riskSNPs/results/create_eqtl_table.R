##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)


#####################
##### Subset of 881 SNPs from PGC
#####################

################
## load eQTLs


load("../eqtl_tables/mergedEqtl_output_hippo_raggr_4features.rda", verbose=TRUE)
sigEqtlHippo_sub = allEqtl[allEqtl$FDR < 0.01,]
load("../eqtl_tables/mergedEqtl_output_dlpfc_raggr_4features.rda", verbose=TRUE)
sigEqtlDlpfc_sub = allEqtl[allEqtl$FDR < 0.01,]

# some genes have extra number, e.g. ENSG###.2
sigEqtlHippo_sub$EnsemblGeneID = ss(sigEqtlHippo_sub$EnsemblGeneID, "\\.")
sigEqtlDlpfc_sub$EnsemblGeneID = ss(sigEqtlDlpfc_sub$EnsemblGeneID, "\\.")


################
## make csv

##### amygdala and sACC #####
dlp = sigEqtlDlpfc_sub
hippo = sigEqtlHippo_sub
dlp$EnsemblGeneID = ss(dlp$EnsemblGeneID, "\\.")
hippo$EnsemblGeneID = ss(hippo$EnsemblGeneID, "\\.")

## snpMap
load("../../genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")
snpMap$hg19POS = paste0(snpMap$CHR,":",snpMap$POS)
snpMap = snpMap[which(rownames(snpMap) %in% c(hippo$snps,dlp$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]


## featMap
load("../../expr_cutoff/rse_gene.Rdata")
load("../../expr_cutoff/rse_exon.Rdata")
load("../../expr_cutoff/rse_jxn.Rdata")
load("../../expr_cutoff/rse_tx.Rdata")
gMap = as.data.frame(rowRanges(rse_gene))[,c("seqnames","start","end","strand","Class")]
eMap = as.data.frame(rowRanges(rse_exon))[,c("seqnames","start","end","strand","Class")]
jMap = as.data.frame(rowRanges(rse_jxn))[,c("seqnames","start","end","strand","Class")]
txMap = as.data.frame(rowRanges(rse_tx))[,c("seqnames","start","end","strand","source")]
txMap$source = "InGen"
# rm(rse_gene, rse_exon, rse_jxn, rse_tx)
colnames(gMap) = colnames(eMap) = colnames(jMap) = colnames(txMap) = 
	c("feat_chr","feat_start","feat_end","strand","Class")
featMap = rbind(rbind(rbind(gMap, eMap),jMap),txMap)
featMap$Type = c(rep("Gene",nrow(gMap)),rep("Exon",nrow(eMap)),rep("Jxn",nrow(jMap)),rep("Tx",nrow(txMap)))

geneMap = as.data.frame(rowRanges(rse_gene))[,c("gencodeID","Symbol","ensemblID","gene_type")]

## put together
snpMap_temp = snpMap[hippo$snps,]
featMap_temp = featMap[hippo$gene,]
geneMap_temp = geneMap[match(hippo$EnsemblGeneID, geneMap$ensemblID),]
hippo2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),hippo)
hippo3 = hippo2[,c(1:4,16,10,5:9,22,21,14,17:20)]
write.csv(hippo3, "raggr_179_snps_hippo_eqtls_fdr01.csv")

snpMap_temp = snpMap[dlp$snps,]
featMap_temp = featMap[dlp$gene,]
geneMap_temp = geneMap[match(dlp$EnsemblGeneID, geneMap$ensemblID),]
dlp2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),dlp)
dlp3 = dlp2[,c(1:4,16,10,5:9,22,21,14,17:20)]
write.csv(dlp3, "raggr_179_snps_dlp_eqtls_fdr01.csv")





