###

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(edgeR)
library(limma)
library(recount)

## load expression data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata")

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene,"[",1)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,"[",1)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate,"[",1)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,"[",1)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,"[",1)
colData(rse_gene)$Kit = ifelse(colData(rse_gene)$mitoRate < 0.05, "Gold", "HMR")

##################
## filter for age and to hippocampus
keepIndex = which(rse_gene$Age > 17 & 
	rse_gene$Kit == "HMR" & rse_gene$Region == "HIPPO")
rse_gene = rse_gene[,keepIndex]
rse_exon = rse_exon[,keepIndex]
rse_jxn = rse_jxn[,keepIndex]
rse_tx = rse_tx[,keepIndex]

## load qSVs and line up
load("count_data/degradation_rse_phase2_hippo.rda")
cov_rse_hippo = cov_rse_hippo[,sapply(rse_gene$SAMPLE_ID, "[", 1)]

###################
## filter low expression
geneIndex = which(rowMeans(getRPKM(rse_gene,"Length")) > 0.2)
rse_gene = rse_gene[geneIndex,]

exonIndex = which(rowMeans(getRPKM(rse_exon,"Length")) > 0.5)
rse_exon = rse_exon[exonIndex,]

rowData(rse_jxn)$bp_length=100 # to trick getRPKM to be RP10M
jxnIndex = which(rowMeans(getRPKM(rse_jxn)) > 0.5 & 
	rowData(rse_jxn)$Class != "Novel")
rse_jxn = rse_jxn[jxnIndex,]

txIndex = which(rowMeans(assays(rse_tx)$tpm) > 0.2)
rse_tx = rse_tx[txIndex,]

## Ns
table(rse_gene$Dx)

## add mds data
mds = read.table("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6.mds",
	header=TRUE,as.is=TRUE, row.names=1)
mds = mds[colData(rse_gene)$BrNum,3:7]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

## model
mod = model.matrix(~Dx + Age + Sex + mitoRate + 
	rRNA_rate + totalAssignedGene + RIN + 
	snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5,
	data = colData(rse_gene))
	
## qSVA
pcaDeg = prcomp(t(log2(assays(cov_rse_hippo)$count + 1)))
k = num.sv(log2(assays(cov_rse_hippo)$count + 1), mod)
qSVs = pcaDeg$x[,1:k]
getPcaVars(pcaDeg)[1:k]
modQsva = cbind(mod, qSVs)

###################
## modeling #######

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts, 
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=TRUE)
fitGene = lmFit(vGene)
ebGene = ebayes(fitGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]
save(outGene, file = "caseControl/dxStats_hippo_filtered_qSVA.rda")

##### Exon ######
dee = DGEList(counts = assays(rse_exon)$counts, 
	genes = rowData(rse_exon))
dee = calcNormFactors(dee)
vExon = voom(dee,modQsva, plot=TRUE)
fitExon = lmFit(vExon)
ebExon = ebayes(fitExon)
eBExon = eBayes(fitExon)
sigExon = topTable(eBExon,coef=2,
	p.value = 1,number=nrow(rse_exon))
outExon = sigExon[rownames(rse_exon),]

##### Junction ######
dje = DGEList(counts = assays(rse_jxn)$counts, 
	genes = rowData(rse_jxn))
dje = calcNormFactors(dje)
vJxn = voom(dje,modQsva, plot=TRUE)
fitJxn = lmFit(vJxn)
ebJxn = ebayes(fitJxn)
eBJxn = eBayes(fitJxn)
sigJxn = topTable(eBJxn,coef=2,
	p.value = 1,number=nrow(rse_jxn))
outJxn = sigJxn[rownames(rse_jxn),]

##### Transcript ######
dte = DGEList(counts = assays(rse_tx)$tpm, 
	genes = rowData(rse_tx))
# dte = calcNormFactors(dte) # already TPM
vTx = voom(dte,modQsva, plot=TRUE)
fitTx = lmFit(vTx)
ebTx = ebayes(fitTx)
eBTx = eBayes(fitTx)
sigTx = topTable(eBTx,coef=2,
	p.value = 1,number=nrow(rse_tx))
outTx = sigTx[rownames(rse_tx),]
save(outGene, outExon, outJxn,outTx,
	file = "caseControl/dxStats_hippo_filtered_qSVA.rda")