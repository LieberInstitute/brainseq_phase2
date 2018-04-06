###

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(edgeR)
library(limma)
library(recount)
library(WGCNA)
allowWGCNAThreads()

## load expression data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata")

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)
colData(rse_gene)$Kit = ifelse(colData(rse_gene)$mitoRate < 0.05, "Gold", "HMR")


###############################################################
################## get qSVA model terms

min(rowMeans(assays(rse_exon)$rpkm)) 					## already cutoff to 0.30
min(rowMeans(assays(rse_jxn)$rp10m))  					## already cutoff to 0.46	


##################
## filter for age and race
keepIndex = which(rse_gene$Age > 17 & rse_gene$Race %in% c("AA", "CAUC") & (
	(rse_gene$Kit=="Gold" & rse_gene$Region=="DLPFC") | (rse_gene$Kit=="HMR" & rse_gene$Region=="HIPPO")) )
rse_gene = rse_gene[,keepIndex]
rse_exon = rse_exon[,keepIndex]
rse_jxn = rse_jxn[,keepIndex]

# ##################
# ## load qSVs
# load("../count_data/degradation_rse_phase2_dlpfc.rda")
# cov_rse_dlpfc = cov_rse_dlpfc[,sapply(rse_gene$SAMPLE_ID, "[", 1)]

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
	
# ## qSVA
# pcaDeg = prcomp(t(log2(assays(cov_rse_dlpfc)$count + 1)))
# k = num.sv(log2(assays(cov_rse_dlpfc)$count + 1), mod)
# qSVs = pcaDeg$x[,1:k]
# getPcaVars(pcaDeg)[1:k]
# modQsva = cbind(mod, qSVs)


###############################################################
############# filter features - exons
############# (same for both regions)

library(parallel)

## WGCNA options
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=8)

exonRpkm = assays(rse_exon)$rpkm
exonMap = rowRanges(rse_exon)
exonExprs = log2(exonRpkm+1) # transform

## normalize
exonExprsAdj = cleaningY(exonExprs, mod, P=3)

## drop based on correlation
exonList = split(as.data.frame(exonExprsAdj), exonMap$gencodeID)

## keep mono exonic
exonListOne = exonList[sapply(exonList,nrow) == 1]
e1 = sapply(exonListOne, rownames)

## cluster and cut exons for multiexonic genes
exonListTwo = exonList[sapply(exonList,nrow) > 1]
exonDistList = mclapply(exonListTwo, function(x) as.dist(1 - cor(t(x))), mc.cores=6)
exonHclustList = mclapply(exonDistList, hclust, mc.cores=6)
exonCutList = mclapply(exonHclustList, cutree, h=0.2, mc.cores=6)
e2 = unlist(mclapply(exonCutList, function(x) names(x[!duplicated(x)]), mc.cores=6))

# extreact
exonExprsAdj2 = exonExprsAdj[c(e1,e2),]
exonMap2 = exonMap[c(e1,e2),]

save(exonExprsAdj2, exonMap2, file="independent_exons.rda")





###############################################################
############# filter features - junctions
############# (same for both regions)

library(parallel)

## WGCNA options
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=8)

jRpkm = assays(rse_jxn)$rp10m
jMap = rowRanges(rse_jxn)
jExprs = log2(jRpkm+1) # transform

## normalize
jExprsAdj = cleaningY(jExprs, mod, P=3)

## drop based on correlation
jList = split(as.data.frame(jExprsAdj), jMap$newGeneID)

## keep mono jic
jListOne = jList[sapply(jList,nrow) == 1]
j1 = sapply(jListOne, rownames)

## cluster and cut js for multijic genes
jListTwo = jList[sapply(jList,nrow) > 1]
jDistList = mclapply(jListTwo, function(x) as.dist(1 - cor(t(x))), mc.cores=6)
jHclustList = mclapply(jDistList, hclust, mc.cores=6)
jCutList = mclapply(jHclustList, cutree, h=0.2, mc.cores=6)
j2 = unlist(mclapply(jCutList, function(x) names(x[!duplicated(x)]), mc.cores=6))

# extreact
jExprsAdj2 = jExprsAdj[c(j1,j2),]
jMap2 = jMap[c(j1,j2),]

save(jExprsAdj2, jMap2, file="independent_jxns.rda")






