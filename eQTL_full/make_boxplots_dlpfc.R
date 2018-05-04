##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)
library(VennDiagram)
library(RColorBrewer)

## load
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_dlpfc.rda")
pd = colData(rse_gene)

## load SNP data ## same 10777 snps used in other regions
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/dlpfc_n167_snps10777_Genotypes.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

### make mds and snp dimensions equal to N
###(repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("../raggr/rAggr_results_881.csv", stringsAsFactors=FALSE)
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

snpMap$maf = rowSums(snp, na.rm=TRUE)/(2*rowSums(!is.na(snp))) 


################
## load table
amyg = read.csv("raggr_31_snps_amyg_eqtls_fdr01.csv", row.names=1)
sacc = read.csv("raggr_31_snps_sacc_eqtls_fdr01.csv", row.names=1)
dlp = read.csv("raggr_31_snps_dlpfc_eqtls_fdr01.csv", row.names=1)


##
pd$Dx = factor(pd$Dx,
	levels = c("Control", "Bipolar"))
mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]), data = pd)


################
## load expression
## dlpfc
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/eqtl_exprs_cutoffs/eQTL_expressed_rse_dlpfc.rda")
geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

## residualize expression		
gExprs = log2(geneRpkm+1)		
gExprs = cleaningY(gExprs, mod, P=2)

eExprs = log2(exonRpkm+1)		
eExprs = cleaningY(eExprs, mod, P=2)

jExprs = log2(jxnRp10m+1)		
jExprs = cleaningY(jExprs, mod, P=2)

tExprs = log2(txTpm+1)		
tExprs = cleaningY(tExprs, mod, P=2)


exprsAdj = rbind(gExprs,eExprs,jExprs,tExprs)
dlp$Symbol = as.character(dlp$Symbol)

pdf("dlpfc_top_eqtl_adj2.pdf", h=6, w=10)
par(mfrow=c(2,3), cex.main=1.2, cex.lab=1.2)
palette(brewer.pal(8,"Spectral"))			
## plot
for (i in 1:nrow(dlp)) {
	symi = dlp[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = dlp[i,"SNP"]
	feati = dlp[i,"gene"]
	p_i = signif(dlp[i,"pvalue"],3)
	typei = dlp[i,"Type.1"]

	boxplot(exprsAdj[feati,] ~ snp[snpi,],
			xlab=snpi, ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (",typei,")"), 
			ylim = c(range(exprsAdj[feati,])), outline=FALSE)
	points(exprsAdj[feati,] ~ jitter(snp[snpi,]+1),
			   pch=21, 
			   bg=as.numeric(snp[snpi,])+2,cex=1.5)
	legend("top",paste0("p=",p_i))
}
dev.off()









