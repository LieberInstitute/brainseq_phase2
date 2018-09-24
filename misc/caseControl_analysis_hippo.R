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
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)
colData(rse_gene)$Kit = ifelse(colData(rse_gene)$mitoRate < 0.05, "Gold", "HMR")

##################
## filter for age and to hippocampus
keepIndex = which(rse_gene$Age > 17 &
	rse_gene$Kit == "HMR" & rse_gene$Region == "HIPPO")
rse_gene = rse_gene[,keepIndex]
rse_exon = rse_exon[,keepIndex]
rse_jxn = rse_jxn[,keepIndex]
rse_tx = rse_tx[,keepIndex]

load("count_data/degradation_rse_phase2_hippo.rda")
cov_rse_hippo = cov_rse_hippo[,sapply(rse_gene$SAMPLE_ID, "[", 1)]

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

## no qSVA
vGene0 = voom(dge,mod, plot=TRUE)
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]
save(outGene, outGene0,file = "caseControl/dxStats_hippo_filtered_qSVA_geneLevel.rda")

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
fitTx = lmFit(assays(rse_tx)$tpm, modQsva)
ebTx = ebayes(fitTx)
eBTx = eBayes(fitTx)
sigTx = topTable(eBTx,coef=2,
	p.value = 1,number=nrow(rse_tx))
outTx = sigTx[rownames(rse_tx),]
save(outGene, outExon, outJxn,outTx,
	file = "caseControl/dxStats_hippo_filtered_qSVA.rda")
