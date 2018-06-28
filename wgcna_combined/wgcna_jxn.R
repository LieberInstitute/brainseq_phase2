## Based on the following script by Emily Burke:
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/wgcna/get_ind_exons_jxns.R

library('jaffelab')
library('SummarizedExperiment')
library('sva')
library('edgeR')
library('limma')
library('recount')
library('RColorBrewer')
library('WGCNA')
library('clusterProfiler')
library('org.Hs.eg.db')
library('devtools')

allowWGCNAThreads()
cores <- 4
enableWGCNAThreads(nThreads = cores)

dir.create('rda', showWarnings = FALSE)
dir.create('pdf', showWarnings = FALSE)
dir.create('top200', showWarnings = FALSE)

## load expression data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata")
# load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata")
# load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata")

## load qSVs
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold.Rdata', verbose = TRUE)

## Drop samples absent in mod and modQsVA
rse_gene <- rse_gene[, keepIndex]
rse_jxn <- rse_jxn[, keepIndex]
# rse_exon <- rse_exon[, keepIndex]
#rse_tx <- rse_tx[, keepIndex]

###############################################################
############# filter features
############# (same for both regions)

# min(rowMeans(assays(rse_exon)$rpkm))                  ## already cutoff to 0.30
min(rowMeans(assays(rse_jxn)$rp10m))  				## already cutoff to 0.46

### Filter to independent features
# load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/wgcna_combined/rda/independent_exons.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/wgcna_combined/rda/independent_jxns.rda",verbose=TRUE)
# eInd = which(rownames(rse_exon) %in% names(exonMap2))
jInd = which(rownames(rse_jxn) %in% names(jMap2))

##################
## filter for age and to hippocampus
keepIndex = which(rse_gene$Age > 17)						# 712
rse_gene = rse_gene[,keepIndex]
# rse_exon = rse_exon[eInd,keepIndex]            # dim: 211455 712
rse_jxn = rse_jxn[jInd,keepIndex]			# dim: 238975 712
mod <- mod[keepIndex, ]
modQsva <- modQsva[keepIndex, ]

dim(modQsva)


###############################################################
############# run WGCNA -- jxn ##################

jRpkm = assays(rse_jxn)$rp10m
jMap = rowRanges(rse_jxn)
jExprs = log2(jRpkm+1)
## normalize
jExprsQsva = cleaningY(jExprs, modQsva, P=3)  ## keep intercept, Dx, Age

################
## thresholding
powers = c(c(1:10), seq(from = 12, to=20, by=2))
if(!file.exists('rda/wgcna_soft_threshold_jxn.rda')) {
    jsftQsva = pickSoftThreshold(t(jExprsQsva),
    	powerVector = powers, verbose = 5)
    save(jsftQsva, file="rda/wgcna_soft_threshold_jxn.rda")
} else {
    load("rda/wgcna_soft_threshold_jxn.rda", verbose = TRUE)
}

## threshold
jsftQsva$powerEstimate
## =

## clustering LIBD
if(!file.exists('rda/jxn_wgcnaModules.rda')) {
    jnetQsva = blockwiseModules(t(jExprsQsva), power = jsftQsva$powerEstimate,
    	TOMType = "unsigned", minModuleSize = 30,
    	reassignThreshold = 0, mergeCutHeight = 0.25,
    	numericLabels = TRUE, pamRespectsDendro = FALSE,
    	saveTOMs = TRUE, 	verbose = 3, maxBlockSize = 30000,
    	saveTOMFileBase = "rda/jxn_qSVA")
    save(jnetQsva, jMap, file = "rda/jxn_wgcnaModules.rda")
} else {
    load("rda/jxn_wgcnaModules.rda", verbose = TRUE)
}

head(table(jnetQsva$colors), 11)  ## 39 clusters. top 10:
#      0      1      2      3      4      5      6      7      8      9     10
# 187830   7526   6955   2840   2585   2455   2088   2014   2013   1686   1610





################
# plot pattern #
################

pd = colData(rse_gene)
pd$roundAge = round(pd$Age)

datME = moduleEigengenes(t(jExprsQsva), jnetQsva$colors)$eigengenes
datME = datME[,-1] ## remove unassigned

pdf("pdf/wgcna_eigengenes_jxn.pdf",h=12,w=12)
par(mfrow=c(2,1), mar=c(5,5,5,2),cex.axis=1.3,cex.lab=1.8,cex.main=2)
palette(brewer.pal(3,"Set1"))
for (i in seq_len(ncol(datME))) {
 plot(pd$Age, datME[,i], cex=1.5, ylim=c(-.37,.67), xaxt="n",
		  pch=16, col=as.factor(pd$Dx),
          ylab="Eigen jxn", xlab="Age", main=paste0("Cluster ",i,"  (", table(jnetQsva$colors)[i+1],")") )
 if (i==1) { legend("topleft", levels(as.factor(pd$Dx)), pch=15, cex=1.5, col=1:2) }
 axis(1, at=unique(pd$roundAge), labels = unique(pd$roundAge))
 # abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
}
dev.off()




################
####   pca  ####
################

jMap = as.data.frame(jMap)
#rownames(jMap) = jMap$gencodeID

moduleInds = split(seq_len(nrow(jExprsQsva)), jnetQsva$colors)
names(moduleInds) = paste0("Cluster_", names(moduleInds))
moduleInds$"Cluster_0" = NULL

jxnExprsQsvaModules = lapply(moduleInds, function(x) exonExprsQsva[x,])
pcaModules = lapply(jxnExprsQsvaModules, function(x) prcomp(t(x)) )

#
# Plot 200 with largest rotation for each cluster
for (k in 1:5) {

pca_k = abs(pcaModules[[k]]$rotation[,1:2])
sigOrderMat = data.frame(PC1 = order(pca_k[,1],decreasing=TRUE)[1:200])
sigOrderMat$gene = rownames(pca_k)[sigOrderMat$PC1]
sigOrderMat$Symbol = jMap[sigOrderMat$gene,"newGeneSymbol"]

maxi = min(100, table(jnetQsva$colors)[k+1] )

## Clusters
pdf(paste0("top200/jxn_cluster",k,"_topRotation.pdf"),h=6,w=12)
par(mar=c(5,6,5,2),cex.axis=1.5,cex.lab=2,cex.main=2)
palette(brewer.pal(3,"Set1"))
for (i in 1:maxi) {
gene = sigOrderMat$gene[i]
symbol = sigOrderMat$Symbol[i]
plot(jExprsQsva[gene,] ~ jitter(pd$Age,.5), xaxt="n",
		pch = 21, bg=as.factor(pd$Dx),
        cex=2,xlab="Age",
        ylab="Residualized log2(Exprs+1)",
		main=paste0(symbol,"\n",gene) )
  axis(1, at=unique(pd$roundAge), labels = unique(pd$roundAge))
  abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
}
dev.off()

}




################
# associations #
################

# split jxns into modules, dropping grey
moduleJxnList_adj = split(gsub('\\..*', '', jMap$newGeneID), jnetQsva$colors)
moduleJxnList_adj = lapply(moduleJxnList_adj, function(x) x[!is.na(x)])
moduleJxnList_adj = moduleJxnList_adj[-1]

## set universe of expressed genes
jxnUniverse = gsub('\\..*', '', as.character(jMap$newGeneID))
jxnUniverse = jxnUniverse[!is.na(jxnUniverse)]

if(!file.exists('rda/wgcna_compareCluster_signed_jxn.rda')) {
    goBP_Adj <- compareCluster(moduleJxnList_adj[1:8], fun = "enrichGO",
        universe = jxnUniverse, OrgDb = org.Hs.eg.db,
        ont = "BP", pAdjustMethod = "BH",
        pvalueCutoff  = .1, qvalueCutoff  = .1,
		readable= TRUE, keyType = 'ENSEMBL')
    goMF_Adj <- compareCluster(moduleJxnList_adj[1:8], fun = "enrichGO",
        universe = jxnUniverse, OrgDb = org.Hs.eg.db,
        ont = "MF", pAdjustMethod = "BH",
        pvalueCutoff  = .1, qvalueCutoff  = .1,
		readable= TRUE, keyType = 'ENSEMBL')
    goCC_Adj <- compareCluster(moduleJxnList_adj[1:8], fun = "enrichGO",
        universe = jxnUniverse, OrgDb = org.Hs.eg.db,
        ont = "CC", pAdjustMethod = "BH",
        pvalueCutoff  = .1, qvalueCutoff  = .1,
		readable= TRUE, keyType = 'ENSEMBL')
    kegg_Adj <- compareCluster(moduleJxnList_adj[1:8], fun = "enrichKEGG",
        universe = jxnUniverse,  pAdjustMethod = "BH",
        pvalueCutoff  = .1, qvalueCutoff  = .1, keyType = 'ENSEMBL')
    save(goBP_Adj, goMF_Adj, goCC_Adj, kegg_Adj, file="rda/wgcna_compareCluster_signed_jxn.rda")
} else {
    load('rda/wgcna_compareCluster_signed_jxn.rda', verbose = TRUE)
}

pdf("pdf/wgcna_enrichments_jxn.pdf",h=6,w=10)
dotplot(goBP_Adj, includeAll="TRUE")
dotplot(goMF_Adj, includeAll="TRUE")
dotplot(goCC_Adj, includeAll="TRUE")
dotplot(kegg_Adj, includeAll="TRUE")
dev.off()



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
