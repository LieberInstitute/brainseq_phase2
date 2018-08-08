##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)


#####################
##### Subset of 881 SNPs from PGC
#####################

################
## load eQTLs

load("eqtl_tables/mergedEqtl_output_hippo_4features.rda", verbose=TRUE)
sigEqtlHippo_sub = allEqtl[allEqtl$FDR < 0.01,]
load("eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)
sigEqtlDlpfc_sub = allEqtl[allEqtl$FDR < 0.01,]
rm(allEqtl)

## Don't include tx eQTLs
sigEqtlHippo_sub = sigEqtlHippo_sub[which(sigEqtlHippo_sub$Type != "Tx"),]
sigEqtlDlpfc_sub = sigEqtlDlpfc_sub[which(sigEqtlDlpfc_sub$Type != "Tx"),]



################
## metrics

## total features
nrow(sigEqtlSacc_sub)  ## 38609
nrow(sigEqtlAmy_sub)   ## 22572
nrow(sigEqtlDlpfc_sub) ## 14669

## per feature
table(sigEqtlSacc_sub$Type)
# Exon 	Gene  Jxn   Tx
# 20014  4463  8641  5491
table(sigEqtlAmy_sub$Type)
# Exon 	 Gene  Jxn   Tx
# 10370  2638  6435  3129
table(sigEqtlDlpfc_sub$Type)
# Exon  Gene  Jxn   Tx
# 6019 1886 3927 2837

## unique ensemblIDs
tapply(sigEqtlSacc_sub$EnsemblGeneID, sigEqtlSacc_sub$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
#  234  175  177  160
tapply(sigEqtlAmy_sub$EnsemblGeneID, sigEqtlAmy_sub$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
#  145  105  130  105
tapply(sigEqtlDlpfc_sub$EnsemblGeneID, sigEqtlDlpfc_sub$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
#  127   89  106  121


################
## make csv

##### amygdala and sACC #####
dlp = sigEqtlDlpfc_sub
amyg = sigEqtlAmy_sub
sacc = sigEqtlSacc_sub
amyg$EnsemblGeneID = ss(amyg$EnsemblGeneID, "\\.")
sacc$EnsemblGeneID = ss(sacc$EnsemblGeneID, "\\.")

## snpMap
load("../../genotype_data/zandiHyde_bipolar_Genotypes_n511.rda")
snpMap$hg19POS = paste0(snpMap$CHR,":",snpMap$POS)
snpMap = snpMap[which(rownames(snpMap) %in% c(amyg$snps,sacc$snps,dlp$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

## featMap
load("../../data/zandiHypde_bipolar_rseTx_n511.rda")
load("../../data/zandiHypde_bipolar_rseJxn_n511.rda")
load("../../data/zandiHypde_bipolar_rseExon_n511.rda")
load("../../data/zandiHypde_bipolar_rseGene_n511.rda")
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
snpMap_temp = snpMap[amyg$snps,]
featMap_temp = featMap[amyg$gene,]
geneMap_temp = geneMap[match(amyg$EnsemblGeneID, geneMap$ensemblID),]
amyg2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),amyg)
amyg3 = amyg2[,c(1:4,16,10,5:9,12:14,17:20)]
write.csv(amyg3, "raggr_suggestive881_snps_amyg_eqtls_fdr01.csv")

snpMap_temp = snpMap[sacc$snps,]
featMap_temp = featMap[sacc$gene,]
geneMap_temp = geneMap[match(sacc$EnsemblGeneID, geneMap$ensemblID),]
sacc2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),sacc)
sacc3 = sacc2[,c(1:4,16,10,5:9,12:14,17:20)]
write.csv(sacc3, "raggr_suggestive881_snps_sacc_eqtls_fdr01.csv")

snpMap_amyg = snpMap
##### DLPFC #####
dlp = sigEqtlDlpfc_sub
dlp$EnsemblGeneID = ss(dlp$EnsemblGeneID, "\\.")

## load SNP data
load("/dcl01/ajaffe/data/lab/brainseq_phase1/genotype_data/brainseq_phase1_Genotypes_n732.rda", verbose=TRUE)
snpMap$hg19POS = paste0(snpMap$CHR,":",snpMap$POS)
snpMap2 = snpMap[which(rownames(snpMap) %in% c(amyg$snps,sacc$snps,dlp$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseTx_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseJxn_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseExon_merged_n732.rda")
load("/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda")
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


snpMap_temp = snpMap[dlp$snps,]
featMap_temp = featMap[dlp$gene,]
geneMap_temp = geneMap[match(dlp$EnsemblGeneID, geneMap$ensemblID),]
dlp2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),dlp)
dlp2 = dlp2[,-which(colnames(dlp2)=="gencodeTx")]

dlp3 = dlp2[,c(2,12,13,14,26,20,15:19,22:24,27:30)]

## fill in NA snps
missInd = which(is.na(dlp3$SNP))
missSnp = dlp$snps[missInd]
snpMap_miss = snpMap_amyg[missSnp,]
count = 1
for (i in missInd) {  
	rownames(dlp3)[i] = rownames(snpMap_miss)[count]
	dlp3[i,1:4] = c(snpMap_miss[count,1:4])
	count = count+1
}

write.csv(dlp3, "raggr_suggestive881_snps_dlpfc_eqtls_fdr01.csv")







