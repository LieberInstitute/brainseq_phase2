# Adapted from https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/casectrl_HIPPO_allFeatures.R
library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library('sessioninfo')

dir.create('pdf', showWarnings = FALSE)
dir.create('rdas', showWarnings = FALSE)

### Run with Gold samples using the qSVs made without the age>17 samples
## and without the HIPPO Gold samples (qsv are HIPPO specific)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata", verbose = TRUE)
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata', verbose = TRUE)

## Load cell type info
load('../methprop_pd.Rdata', verbose = TRUE)

## Plug the cell type info
colData(rse_gene) <- colData(rse_exon) <- colData(rse_jxn) <- colData(rse_tx) <- pd

## Dropped Fetal_replicating + Fetal_quiescent to make the model tractable
cellMod <- with(pd[keepIndex, ], model.matrix( ~ OPC + Neurons + Astrocytes + Oligodendrocytes + Microglia + Endothelial))
modQsva <- cbind(modQsva, cellMod[, -1])
stopifnot(is.fullrank(modQsva))
mod <- cbind(mod, cellMod[, -1])
stopifnot(is.fullrank(mod))

## Drop samples absent in mod and modQsVA
rse_gene <- rse_gene[, keepIndex]
rse_jxn <- rse_jxn[, keepIndex]
rse_exon <- rse_exon[, keepIndex]
rse_tx <- rse_tx[, keepIndex]

## Keep region-specific samples
keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "HIPPO")
rse_gene <- rse_gene[, keepIndex]
rse_jxn <- rse_jxn[, keepIndex]
rse_exon <- rse_exon[, keepIndex]
rse_tx <- rse_tx[, keepIndex]

mod <- mod[keepIndex, ]
modQsva <- modQsva[keepIndex, ]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)
pdf('pdf/hippo_voom_qsva_noHGoldQSV_matchHIPPO_gene.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]


##### Exon ######
dee = DGEList(counts = assays(rse_exon)$counts,
	genes = rowData(rse_exon))
dee = calcNormFactors(dee)
pdf('pdf/hippo_voom_qsva_noHGoldQSV_matchHIPPO_exon.pdf', useDingbats = FALSE)
vExon = voom(dee,modQsva, plot=TRUE)
dev.off()
fitExon = lmFit(vExon)
eBExon = eBayes(fitExon)
sigExon = topTable(eBExon,coef=2,
	p.value = 1,number=nrow(rse_exon))
outExon = sigExon[rownames(rse_exon),]

##### Junction ######
dje = DGEList(counts = assays(rse_jxn)$counts,
	genes = rowData(rse_jxn))
dje = calcNormFactors(dje)
pdf('pdf/hippo_voom_qsva_noHGoldQSV_matchHIPPO_jxn.pdf', useDingbats = FALSE)
vJxn = voom(dje,modQsva, plot=TRUE)
dev.off()
fitJxn = lmFit(vJxn)
eBJxn = eBayes(fitJxn)
sigJxn = topTable(eBJxn,coef=2,
	p.value = 1,number=nrow(rse_jxn))
outJxn = sigJxn[rownames(rse_jxn),]

##### Transcript ######
fitTx = lmFit(log2(assays(rse_tx)$tpm + 0.5), modQsva)
eBTx = eBayes(fitTx)
sigTx = topTable(eBTx,coef=2,
	p.value = 1,number=nrow(rse_tx))
outTx = sigTx[rownames(rse_tx),]
outTx <- cbind(outTx, rowData(rse_tx))


save(outGene, outExon, outJxn,outTx,
	file = "rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
