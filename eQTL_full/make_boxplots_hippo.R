####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(RColorBrewer)

######################
### load data ####
######################

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata")


## keep adult samples & correct region
keepInd = which(colData(rse_gene)$Age > 13 & colData(rse_gene)$Region == "HIPPO")
rse_gene = rse_gene[,keepInd]
rse_exon = rse_exon[,keepInd]
rse_jxn = rse_jxn[,keepInd]
rse_tx = rse_tx[,keepInd]

## extract pd and rpkms
pd = colData(rse_gene)



######################
### snp data ####
######################

## load SNP data
load("../genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")

### make mds and snp dimensions equal to N
###(repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum

## drop SNPs not mapping to hg38
keepIndex = which(!is.na(snpMap$chr_hg38))
snpMap = snpMap[keepIndex,]
snp = snp[keepIndex,]



################
## load table
load("eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)
dlp = allEqtl

load("eqtl_tables/mergedEqtl_output_hippo_4features.rda", verbose=TRUE)
hippo = allEqtl

hippo = allEqtl[which(allEqtl$Type=="Gene"),]
hippo$Symbol = as.character(hippo$Symbol)

################
## load expression

pd$Dx = factor(pd$Dx,
	levels = c("Control", "Schizo"))

mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]),
	data = pd)
colnames(mod)[4:8] = colnames(mds)[1:5]

geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

load('eqtl_tables/rdas/pcs_hippo_4features_filtered_over13.rda', verbose = TRUE)

## residualize expression        
gExprs = log2(geneRpkm+1)        
gExprs = cleaningY(gExprs, cbind(mod, genePCs), P=1)

eExprs = log2(exonRpkm+1)        
eExprs = cleaningY(eExprs, cbind(mod, exonPCs), P=1)

jExprs = log2(jxnRp10m+1)        
jExprs = cleaningY(jExprs, cbind(mod, jxnPCs), P=1)

tExprs = log2(txTpm+1)        
tExprs = cleaningY(tExprs, cbind(mod, txPCs), P=1)


exprsAdj = rbind(gExprs,eExprs,jExprs,tExprs)



pdf("hippo_top_eqtl_gene.pdf", h=6, w=6)
par(mar=c(5,5,5,2), cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
palette(brewer.pal(8,"Spectral"))			
## plot
for (i in 1:25) {
	symi = hippo[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = hippo[i,"snps"]
	feati = hippo[i,"gene"]
	p_i = signif(hippo[i,"pvalue"],3)

	boxplot(gExprs[feati,] ~ snp[snpi,],
			xlab=snpi, ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (Gene)"), 
			ylim = c(range(gExprs[feati,])), outline=FALSE)
	points(gExprs[feati,] ~ jitter(snp[snpi,]+1),
			   pch=21, 
			   bg=as.numeric(snp[snpi,])+2,cex=1.5)
	legend("top",paste0("p=",p_i))
}
dev.off()



#################
### PGC index snps

pgc = read.csv("../eQTL_GWAS_riskSNPs/41588_2018_59_MOESM3_ESM.csv", stringsAsFactors=FALSE)
index = pgc$"Index.SNP..dbSNP.b141."

hippo2 = allEqtl[which(allEqtl$snps %in% index),]
hippo2 = hippo2[order(hippo2$FDR, decreasing=FALSE),]

dlp2 = dlp[which(dlp$snps %in% index),]
dlp2 = dlp2[order(dlp2$FDR, decreasing=FALSE),]



pdf("hippo_top_eqtl_PGC_indexSNPs.pdf", h=6, w=6, useDingbats=FALSE)
par(mar=c(5,5,5,2), cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
palette(brewer.pal(8,"Spectral"))			
## plot
for (i in 1:15) {
	symi = hippo2[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = hippo2[i,"snps"]
	feati = hippo2[i,"gene"]
	p_i = signif(hippo2[i,"pvalue"],3)
	typei = hippo2[i,"Type"]
	
	boxplot(exprsAdj[feati,] ~ snp[snpi,],
			xlab=snpi, ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (",typei,")"), 
			ylim = c(range(exprsAdj[feati,])), outline=FALSE)
	points(exprsAdj[feati,] ~ jitter(snp[snpi,]+1),
			   pch=21, 
			   bg=as.numeric(snp[snpi,])+2,cex=1.5)
	legend("topleft",paste0("p=",p_i))
}
dev.off()





