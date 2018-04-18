###

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(edgeR)
library(limma)
library(recount)
library(WGCNA)
library(RColorBrewer)

## load expression data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata")
# load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata")

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)
colData(rse_gene)$Kit = ifelse(colData(rse_gene)$mitoRate < 0.05, "Gold", "HMR")
stopifnot(identical(colnames(rse_gene), colnames(rse_exon)))
stopifnot(identical(colnames(rse_gene), colnames(rse_jxn)))

###############################################################
############# filter features
############# (same for both regions)

min(rowMeans(assays(rse_exon)$rpkm)) 	 			## already cutoff to 0.30	
min(rowMeans(assays(rse_jxn)$rp10m))  				## already cutoff to 0.46		

### Filter to independent features
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/wgcna/independent_exons.rda",verbose=TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/wgcna/independent_jxns.rda",verbose=TRUE)
eInd = which(rownames(rse_exon) %in% names(exonMap2))
jInd = which(rownames(rse_jxn) %in% names(jMap2))

##################
## filter for age and to hippocampus
keepIndex = which(rse_gene$Age > 17 & rse_gene$Race %in% c("AA", "CAUC") & 
	rse_gene$Kit == "Gold" & rse_gene$Region == "DLPFC")						# 357 
rse_gene = rse_gene[,keepIndex]
rse_exon = rse_exon[eInd,keepIndex]			# dim: 211455 357
rse_jxn = rse_jxn[jInd,keepIndex]			# dim: 238975 357


##################
## load qSVs
load("../count_data/degradation_rse_phase2_dlpfc.rda")
cov_rse_dlpfc = cov_rse_dlpfc[,sapply(rse_gene$SAMPLE_ID, "[", 1)]

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
pcaDeg = prcomp(t(log2(assays(cov_rse_dlpfc)$count + 1)))
k = num.sv(log2(assays(cov_rse_dlpfc)$count + 1), mod) 		# k = 16
qSVs = pcaDeg$x[,1:k]
getPcaVars(pcaDeg)[1:k]
modQsva = cbind(mod, qSVs)



###############################################################
############# run WGCNA -- exon ##################

exonRpkm = assays(rse_exon)$rpkm
exonMap = rowRanges(rse_exon)
exonExprs = log2(exonRpkm+1)
## normalize
exonExprsQsva = cleaningY(exonExprs, modQsva, P=3)  ## keep intercept, Dx, Age

################
## thresholding
powers = c(c(1:10), seq(from = 12, to=20, by=2))
esftQsva = pickSoftThreshold(t(exonExprsQsva), 
	powerVector = powers, verbose = 5)
save(esftQsva, file="rdas/wgcna_soft_threshold_DLPFC_exon.rda")
# load("rdas/wgcna_soft_threshold_DLPFC_exon.rda")

## threshold
esftQsva$powerEstimate
## = 

## clustering LIBD
enetQsva = blockwiseModules(t(exonExprsQsva), power = esftQsva$powerEstimate,
	TOMType = "unsigned", minModuleSize = 30,
	reassignThreshold = 0, mergeCutHeight = 0.25,
	numericLabels = TRUE, pamRespectsDendro = FALSE,
	saveTOMs = TRUE, 	verbose = 3, maxBlockSize = 30000,
	saveTOMFileBase = "rdas/DLPFC_exon_qSVA")
save(enetQsva, exonMap, file = "rdas/DLPFC_exon_wgcnaModules.rda")

####

# load("rdas/DLPFC_exon_wgcnaModules.rda", verbose=TRUE)

table(enetQsva$colors)  ## 39 clusters. top 10:
#    0     1     2     3     4     5     6     7     8     9    10
# 14928   913   775   642   632   612   605   550   406   397   359



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
jsftQsva = pickSoftThreshold(t(jExprsQsva), 
	powerVector = powers, verbose = 5)
save(jsftQsva, file="rdas/wgcna_soft_threshold_DLPFC_jxn.rda")
# load("rdas/wgcna_soft_threshold_DLPFC_jxn.rda")

## threshold
jsftQsva$powerEstimate
## = 

## clustering LIBD
jnetQsva = blockwiseModules(t(jExprsQsva), power = jsftQsva$powerEstimate,
	TOMType = "unsigned", minModuleSize = 30,
	reassignThreshold = 0, mergeCutHeight = 0.25,
	numericLabels = TRUE, pamRespectsDendro = FALSE,
	saveTOMs = TRUE, 	verbose = 3, maxBlockSize = 30000,
	saveTOMFileBase = "rdas/DLPFC_jxn_qSVA")
save(jnetQsva, jMap, file = "rdas/DLPFC_jxn_wgcnaModules.rda")

####

# load("rdas/DLPFC_jxn_wgcnaModules.rda", verbose=TRUE)

table(jnetQsva$colors)  ## 39 clusters. top 10:
#    0     1     2     3     4     5     6     7     8     9    10
# 14928   913   775   642   632   612   605   550   406   397   359














# ################
# # plot pattern #
# ################

# pd = colData(rse_gene)[,c(1,27:50)]
# pd$roundAge = round(pd$Age)

# datME = moduleEigengenes(t(geneExprsQsva), netQsva$colors)$eigengenes
# datME = datME[,-1] ## remove unassigned

# pdf("wgcna_eigenegenes_hippo_gene.pdf",h=12,w=12)
# par(mfrow=c(2,1), mar=c(5,5,5,2),cex.axis=1.3,cex.lab=1.8,cex.main=2)
# palette(brewer.pal(3,"Set1"))
# for (i in 1:max(netQsva$colors)) {
 # plot(pd$Age, datME[,i], cex=1.5, ylim=c(-.28,.72), xaxt="n",
		  # pch=16, col=as.factor(pd$Dx),
          # ylab="Eigengene", xlab="Age", main=paste0("Cluster ",i,"  (", table(netQsva$colors)[i+1],")") )
 # if (i==1) { legend("topleft", levels(as.factor(pd$Dx)), pch=15, cex=1.5, col=1:2) }
 # axis(1, at=unique(pd$roundAge), labels = unique(pd$roundAge))
 # # abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
# }
# dev.off()


# ################
# ####   pca  ####
# ################

# geneMap = as.data.frame(geneMap)
# rownames(geneMap) = geneMap$gencodeID

# moduleInds = split(1:nrow(geneExprsQsva), netQsva$colors)
# names(moduleInds) = paste0("Cluster_", names(moduleInds))
# moduleInds$"Cluster_0" = NULL

# geneExprsQsvaModules = lapply(moduleInds, function(x) geneExprsQsva[x,])
# pcaModules = lapply(geneExprsQsvaModules, function(x) prcomp(t(x)) )

# #
# # Plot 200 with largest rotation for each cluster
# for (k in 1:5) {

# pca_k = abs(pcaModules[[k]]$rotation[,1:2])
# sigOrderMat = data.frame(PC1 = order(pca_k[,1],decreasing=TRUE)[1:200])
# sigOrderMat$gene = rownames(pca_k)[sigOrderMat$PC1]
# sigOrderMat$Symbol = geneMap[sigOrderMat$gene,"Symbol"]

# maxi = min(100, table(netQsva$colors)[k+1] )

# ## Clusters
# pdf(paste0("top200/hippo_gene_cluster",k,"_topRotation.pdf"),h=6,w=12)
# par(mar=c(5,6,5,2),cex.axis=1.5,cex.lab=2,cex.main=2)
# palette(brewer.pal(3,"Set1"))
# for (i in 1:maxi) {
# gene = sigOrderMat$gene[i]
# symbol = sigOrderMat$Symbol[i]
# plot(geneExprsQsva[gene,] ~ jitter(pd$Age,.5), xaxt="n",
		# pch = 21, bg=as.factor(pd$Dx),
        # cex=2,xlab="Day",
        # ylab="Residualized log2(Exprs+1)",
		# main=paste0(symbol,"\n",gene) )
  # axis(1, at=unique(pd$roundAge), labels = unique(pd$roundAge))
  # abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
# }
# dev.off()

# }





# ################
# # associations #
# ################

# # gene set associations
# library(clusterProfiler)
# library(org.Hs.eg.db)

# # split genes into modules, dropping grey
# moduleGeneList_adj = split(geneMap$EntrezID, netQsva$colors)
# moduleGeneList_adj = lapply(moduleGeneList_adj, function(x) x[!is.na(x)])
# moduleGeneList_adj = moduleGeneList_adj[-1]

# ## set universe of expressed genes
# geneUniverse = as.character(geneMap$EntrezID)
# geneUniverse = geneUniverse[!is.na(geneUniverse)]

# ############################## 
# ## run enrichment analysis ###
# ##############################

# ## and adjusted
# goBP_Adj <- compareCluster(moduleGeneList_adj[1:8], fun = "enrichGO",
                # universe = geneUniverse, OrgDb = org.Hs.eg.db,
                # ont = "BP", pAdjustMethod = "BH",
                # pvalueCutoff  = .1, qvalueCutoff  = .1,
				# readable= TRUE)
# goMF_Adj <- compareCluster(moduleGeneList_adj[1:8], fun = "enrichGO",
                # universe = geneUniverse, OrgDb = org.Hs.eg.db,
                # ont = "MF", pAdjustMethod = "BH",
                # pvalueCutoff  = .1, qvalueCutoff  = .1,
				# readable= TRUE)
# goCC_Adj <- compareCluster(moduleGeneList_adj[1:8], fun = "enrichGO",
                # universe = geneUniverse, OrgDb = org.Hs.eg.db,
                # ont = "CC", pAdjustMethod = "BH",
                # pvalueCutoff  = .1, qvalueCutoff  = .1,
				# readable= TRUE)
# kegg_Adj <- compareCluster(moduleGeneList_adj[1:8], fun = "enrichKEGG",
                # universe = geneUniverse,  pAdjustMethod = "BH",
                # pvalueCutoff  = .1, qvalueCutoff  = .1)
# # save(goBP_Adj, goMF_Adj, goCC_Adj, kegg_Adj, file="wgcna_compareCluster_signed.rda")
				
# pdf("wgcna_enrichments_hippo_gene.pdf",h=6,w=10)
# dotplot(goBP_Adj, includeAll="TRUE")
# dotplot(goMF_Adj, includeAll="TRUE")
# dotplot(goCC_Adj, includeAll="TRUE")
# dotplot(kegg_Adj, includeAll="TRUE")
# dev.off()				
		

















