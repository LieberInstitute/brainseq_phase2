###

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(edgeR)
library(limma)

## load gene
load("count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n442.rda")
load("count_data/hippo_brainseq_phase2_hg38_rseExon_merged_n442.rda")
load("count_data/hippo_brainseq_phase2_hg38_rseJxn_merged_n442.rda")
load("count_data/hippo_brainseq_phase2_hg38_rseTx_merged_n442.rda")
colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene,"[",1)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,"[",1)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate,"[",1)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,"[",1)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,"[",1)
colData(rse_gene)$Kit = ifelse(colData(rse_gene)$mitoRate < 0.05, "Gold", "HMR")
getRPKM = recount::getRPKM

## load qSVs
load("count_data/degradation_rse_phase2_hippo.rda")

## filter for age
keepIndex = which(rse_gene$Age > 17 & 
	rse_gene$Kit == "HMR")
rse_gene = rse_gene[,keepIndex]
rse_exon = rse_exon[,keepIndex]
rse_jxn = rse_jxn[,keepIndex]
rse_tx = rse_tx[,keepIndex]
cov_rse_hippo = cov_rse_hippo[,keepIndex]

## filter low expression
geneIndex = which(rowMeans(getRPKM(rse_gene,"Length")) > 0.2)
rse_gene = rse_gene[geneIndex,]

exonIndex = which(rowMeans(getRPKM(rse_exon,"Length")) > 0.5)
rse_exon = rse_exon[exonIndex,]

jxnIndex = which(rowMeans(getRPM(rse_jxn)) > 0.5 & rowData(jMap)$Class != "Novel")
rse_gene = rse_gene[geneIndex,]

## Ns
table(rse_gene$Dx)

## add mds data
mds = read.table("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n546_maf05_geno10_hwe1e6.mds",
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

## modeling
dge = DGEList(counts = assays(rse_gene)$counts, 
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)

## mean-variance
vGene = voom(dge,modQsva, plot=TRUE)

## do analysis
fitGene = lmFit(vGene)
ebGene = ebayes(fitGene)

## top table
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]
save(outGene, file = "caseControl/dxStats_hippo_filtered_qSVA.rda")