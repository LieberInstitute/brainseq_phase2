##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)
library(RColorBrewer)
library(VennDiagram)

#####################
##### Subset of 881 SNPs from PGC
#####################

################
## load eQTLs

load("eqtl_tables/mergedEqtl_output_hippo_raggr_4features.rda", verbose=TRUE)
sigEqtlHippo_sub = allEqtl[allEqtl$FDR < 0.01,]
load("eqtl_tables/mergedEqtl_output_dlpfc_raggr_4features.rda", verbose=TRUE)
sigEqtlDlpfc_sub = allEqtl[allEqtl$FDR < 0.01,]

# some genes have extra number, e.g. ENSG###.2
sigEqtlHippo_sub$EnsemblGeneID = ss(sigEqtlHippo_sub$EnsemblGeneID, "\\.")
sigEqtlDlpfc_sub$EnsemblGeneID = ss(sigEqtlDlpfc_sub$EnsemblGeneID, "\\.")

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("rAggr_results_179.csv", stringsAsFactors=FALSE)	# 10,981 snps
riskLoci_full = riskLoci
colnames(riskLoci) = colnames(riskLoci_full) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

################
## metrics

## unique SNPs and index SNPs
length(unique(sigEqtlHippo_sub$snps))
length(unique(sigEqtlHippo_sub$snps[which(sigEqtlHippo_sub$snps %in% riskLoci$SNP1_Name)]))
length(unique(sigEqtlDlpfc_sub$snps))
length(unique(sigEqtlDlpfc_sub$snps[which(sigEqtlDlpfc_sub$snps %in% riskLoci$SNP1_Name)]))

## unique features by type
tapply(sigEqtlHippo_sub$gene, sigEqtlHippo_sub$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
# 857  123  507  244
tapply(sigEqtlDlpfc_sub$gene, sigEqtlDlpfc_sub$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
#  1363  171  659  332

## total snp-feature pairs
## from xx genes
nrow(sigEqtlHippo_sub)   ## 66,923
length(unique(sigEqtlHippo_sub$EnsemblGeneID))
nrow(sigEqtlDlpfc_sub) ## 106,438
length(unique(sigEqtlDlpfc_sub$EnsemblGeneID))

## total features by type
table(sigEqtlHippo_sub$Type)
# Exon 	 Gene  Jxn   Tx
# 33122  4795 19377  9629
table(sigEqtlDlpfc_sub$Type)
# Exon  Gene  Jxn   Tx
# 57986  7358 27876 13218

## unique ensemblIDs
tapply(sigEqtlHippo_sub$EnsemblGeneID, sigEqtlHippo_sub$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
#  216  123  165  170
tapply(sigEqtlDlpfc_sub$EnsemblGeneID, sigEqtlDlpfc_sub$Type, function(x) length(unique(x)))
# Exon Gene  Jxn   Tx
#  270  171  204  216


pal = brewer.pal(8,"Set1")
venn.diagram(list(Hippo = unique(sigEqtlHippo_sub$snps), DLPFC = unique(sigEqtlDlpfc_sub$snps)), 
	fill = pal[1:2], main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_SNP.png")
venn.diagram(list(Hippo = unique(sigEqtlHippo_sub$snps[which(sigEqtlHippo_sub$snps %in% riskLoci$SNP1_Name)]), 
				DLPFC = unique(sigEqtlDlpfc_sub$snps[which(sigEqtlDlpfc_sub$snps %in% riskLoci$SNP1_Name)])), 
	fill = pal[1:2], main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_IndexSNP.png")

venn.diagram(list(Hippo = unique(sigEqtlHippo_sub$gene[which(sigEqtlHippo_sub$Type=="Gene")]), 
				DLPFC = unique(sigEqtlDlpfc_sub$gene[which(sigEqtlDlpfc_sub$Type=="Gene")])), 
	fill = pal[1:2], main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_feats_Gene.png")	
venn.diagram(list(Hippo = unique(sigEqtlHippo_sub$gene[which(sigEqtlHippo_sub$Type=="Exon")]), 
				DLPFC = unique(sigEqtlDlpfc_sub$gene[which(sigEqtlDlpfc_sub$Type=="Exon")])), 
	fill = pal[1:2], main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_feats_Exon.png")	
venn.diagram(list(Hippo = unique(sigEqtlHippo_sub$gene[which(sigEqtlHippo_sub$Type=="Jxn")]), 
				DLPFC = unique(sigEqtlDlpfc_sub$gene[which(sigEqtlDlpfc_sub$Type=="Jxn")])), 
	fill = pal[1:2], main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_feats_Jxn.png")	
venn.diagram(list(Hippo = unique(sigEqtlHippo_sub$gene[which(sigEqtlHippo_sub$Type=="Tx")]), 
				DLPFC = unique(sigEqtlDlpfc_sub$gene[which(sigEqtlDlpfc_sub$Type=="Tx")])), 
	fill = pal[1:2], main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_feats_Tx.png")	

	
## can you guys also tabulate how many of the SNPs only associate with 1 gene
h = sigEqtlHippo_sub[order(sigEqtlHippo_sub$snps),]
h_snps = data.frame(snps = unique(h$snps), ngenes = NA)
for (i in 1:nrow(h_snps)) {
	s = h_snps$snps[i]
	h_snps$ngenes[i] = length(unique(h$EnsemblGeneID[which(h$snps == s)]))
}
table(h_snps$ngenes)
   # 1    2    3    4    5    6    7    8    9   10   11   12	13
# 1982 1353  653  580  317   91  229   28  106   47   58   59    7

d = sigEqtlDlpfc_sub[order(sigEqtlDlpfc_sub$snps),]
d_snps = data.frame(snps = unique(d$snps), ngenes = NA)
for (i in 1:nrow(d_snps)) {
	s = d_snps$snps[i]
	d_snps$ngenes[i] = length(unique(d$EnsemblGeneID[which(d$snps == s)]))
}
table(d_snps$ngenes)
   # 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
# 2447 1790  679  686  329  292  107  136   50   49   70   73   36   28    7    1


################
## make csv

##### dlp and hippo #####
hippo = sigEqtlHippo_sub
dlp = sigEqtlDlpfc_sub
# dlp$EnsemblGeneID = ss(dlp$EnsemblGeneID, "\\.")
# hippo$EnsemblGeneID = ss(hippo$EnsemblGeneID, "\\.")

## snpMap
load("../genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")
snpMap$hg19POS = paste0(snpMap$CHR,":",snpMap$POS)
snpMap = snpMap[which(rownames(snpMap) %in% c(hippo$snps,dlp$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

## featMap
load("../expr_cutoff/rse_gene.Rdata")
load("../expr_cutoff/rse_exon.Rdata")
load("../expr_cutoff/rse_jxn.Rdata")
load("../expr_cutoff/rse_tx.Rdata")
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
hippo3 = hippo2[,c(1:4,16,10,5:9,12:14,17:20)]
write.csv(hippo3, "raggr179snps_hippo_eqtls_fdr01.csv")

snpMap_temp = snpMap[dlp$snps,]
featMap_temp = featMap[dlp$gene,]
geneMap_temp = geneMap[match(dlp$EnsemblGeneID, geneMap$ensemblID),]
dlp2 = cbind(cbind(cbind(snpMap_temp,featMap_temp),geneMap_temp),dlp)
dlp3 = dlp2[,c(1:4,16,10,5:9,12:14,17:20)]
write.csv(dlp3, "raggr179snps_dlp_eqtls_fdr01.csv")







