##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)
library(VennDiagram)
library(RColorBrewer)

## load
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

## load SNP data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

# ## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
# snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
# snpMap = snpMap[-snpInd,]
# snp = snp[-snpInd,]

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("rAggr_results_179.csv", stringsAsFactors=FALSE)	# 10,981 snps
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 
riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

# ## keep SNPs from list
# keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS2)	# keep 9735 snps from snpMap
# snpMap = snpMap[keepIndex,]
# snp = snp[keepIndex,]

# snpMap$maf = rowSums(snp, na.rm=TRUE)/(2*rowSums(!is.na(snp))) 

#####################
# filter brain region
# make mds and snp dimensions equal to N
# (repeat rows or columns for BrNum replicates)
pd = colData(rse_gene)

mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum



################
## load table

load("eqtl_tables/mergedEqtl_output_hippo_raggr_4features.rda", verbose=TRUE)

# hippo = read.csv("raggr_179_snps_hippo_eqtls_fdr01.csv", row.names=1)
hippo = read.csv("raggr_179_snps_hippo_eqtls_fdr01.csv", row.names=1)
hippo$Symbol = as.character(hippo$Symbol)
hippo$SNP = as.character(hippo$SNP)
hippo$gene = as.character(hippo$gene)
hippo$IndexSNP = as.character(hippo$IndexSNP)

load("eqtl_tables/rdas/pcs_4features_hippo.rda", verbose=TRUE)
##
modG = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]) + genePCs, data = pd)
modE = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]) + exonPCs, data = pd)
modJ = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]) + jxnPCs, data = pd)
modT = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]) + txPCs, data = pd)

################
## load expression
## hippo
geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

## residualize expression		
gExprs = log2(geneRpkm+1)		
gExprs = cleaningY(gExprs, modG, P=1)

eExprs = log2(exonRpkm+1)		
eExprs = cleaningY(eExprs, modE, P=1)

jExprs = log2(jxnRp10m+1)		
jExprs = cleaningY(jExprs, modJ, P=1)

tExprs = log2(txTpm+1)		
tExprs = cleaningY(tExprs, modT, P=1)

exprsAdj = rbind(gExprs,eExprs,jExprs,tExprs)

hippo2 = hippo[order(hippo$pvalue, decreasing=FALSE),]

pdf("hippo_top_eqtl_adj.pdf", h=6, w=10)
par(mfrow=c(2,3), cex.main=1.2, cex.lab=1.2)
palette(brewer.pal(6,"Spectral"))			
## plot
for (i in 1:100) {
	symi = hippo2[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = hippo2[i,"SNP"]
	feati = hippo2[i,"gene"]
	p_i = signif(hippo2[i,"pvalue"],3)
	typei = hippo2[i,"Type"]

	boxplot(exprsAdj[feati,] ~ snp[snpi,],
			xlab=snpi, ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (",typei,")"), 
			ylim = c(range(exprsAdj[feati,])), outline=FALSE)
	points(exprsAdj[feati,] ~ jitter(snp[snpi,]+1),
			   pch=21, 
			   bg=as.numeric(snp[snpi,])+4,cex=1.5)
	legend("top",paste0("p=",p_i))
}

dev.off()


#### plot proxy and index next to each other
pdf("hippo_top_eqtl_adj_index.pdf", h=8, w=8)
par(mfrow=c(2,2), cex.main=1.2, cex.lab=1.2)
palette(brewer.pal(6,"Spectral"))			
## plot
for (i in 1:10) {
	symi = hippo2[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = hippo2[i,"SNP"]
	snpicoord = hippo2[i,"hg19POS"]
	indexi = hippo2[i,"IndexSNP"]
	indexicoord = hippo2[i,"IndexSNP_hg19POS"]
	feati = hippo2[i,"gene"]
	pindex_i = signif(hippo2$pvalue[which(hippo2$SNP==indexi & hippo2$gene==feati)],3)
	p_i = signif(hippo2[i,"pvalue"],3)
	typei = hippo2[i,"Type"]

	boxplot(exprsAdj[feati,] ~ snp[snpi,],
			xlab=paste0(snpi,"\n",snpicoord), ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (",typei,")"), 
			ylim = c(range(exprsAdj[feati,])), outline=FALSE)
	points(exprsAdj[feati,] ~ jitter(snp[snpi,]+1),
			   pch=21, 
			   bg=as.numeric(snp[snpi,])+4,cex=1.5)
	legend("top",paste0("p=",p_i))
	
	boxplot(exprsAdj[feati,] ~ snp[indexi,],
			xlab=paste0(indexi,"\n",indexicoord), ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (",typei,")"), 
			ylim = c(range(exprsAdj[feati,])), outline=FALSE)
	points(exprsAdj[feati,] ~ jitter(snp[indexi,]+1),
			   pch=21, 
			   bg=as.numeric(snp[indexi,])+4,cex=1.5)
	legend("top", paste0("p=",pindex_i))
}
dev.off()



