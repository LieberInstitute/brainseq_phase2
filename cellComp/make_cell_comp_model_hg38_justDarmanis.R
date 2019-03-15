###

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(genefilter)
library(RColorBrewer)

# load gene-level darmanis
load("/dcl01/ajaffe/data/lab/singleCell/Darmanis/rna-seq-pipeline_run2/rse_gene_Darmanis_scRNASeq_Darmanis_n466.Rdata")
rse_geneQuake = rse_gene
rm(rse_gene)

# phenotype
load("/dcl01/lieber/ajaffe/PublicData/Brain/Darmanis_Quake/rdas/phenotype_quake_processed.rda")
colData(rse_geneQuake) = cbind(colData(rse_geneQuake), 
	pd[match(colnames(rse_geneQuake), pd$RunName),2:14])

rse_geneQuake$Cell_type[rse_geneQuake$Cell_type == "hybrids"] = "hybrid"
## exclude hybrid
rse_geneQuake = rse_geneQuake[,rse_geneQuake$Cell_type != "hybrid"]

# get expression
geneExprsQuake = log2(getRPKM(rse_geneQuake, "Length")+1)
	
## make factor
rse_geneQuake$Cell_type = factor(rse_geneQuake$Cell_type,
	levels =c("Fetal_replicating", "Fetal_quiescent", 
		"OPC", "Neurons", "Astrocytes",
		"Oligodendrocytes", "Microglia", "Endothelial"))

## #split by age cat and region					
group = rse_geneQuake$Cell_type
tIndexes <- splitit(group)

tstatList <- lapply(tIndexes, function(i) {
	x <- rep(0, ncol(geneExprsQuake))
	x[i] <- 1
	return(genefilter::rowttests(geneExprsQuake, factor(x)))
})

numProbes=20
probeList <- lapply(tstatList, function(x) {
	y <- x[which(x[, "p.value"] < 1e-15), ]
	yUp <- y[order(y[, "dm"], decreasing = FALSE),] # signs are swapped
	rownames(yUp)[1:numProbes]
})

## filter
trainingProbes <- unique(unlist(probeList))
mergeMarkerExprs <- geneExprsQuake[trainingProbes, ]
mergeMarkerMeanExprs <- colMeans(mergeMarkerExprs)

form <- as.formula(sprintf("y ~ %s - 1", paste(levels(group),collapse = "+")))
phenoDF <- as.data.frame(model.matrix(~group - 1))
colnames(phenoDF) <- sub("^group", "", colnames(phenoDF))

## try z-score
mergeMarkerExprsZ = scale(mergeMarkerExprs)
mergeMarkerMeanExprsZ = colMeans(mergeMarkerExprsZ)

## do calibration
coefEsts <- minfi:::validationCellType(Y = mergeMarkerExprsZ, 
	pheno = phenoDF, modelFix = form)$coefEsts

save(coefEsts, mergeMarkerMeanExprsZ, 
	file = "singleCell_quake_coefEsts_calibration_Zscale_hg38.rda")
write.csv(coefEsts, file = "singleCell_quake_coefEsts_calibration_Zscale_hg38.csv")

###################
## plots ##########
###################

## make plot
library(lattice)
theSeq = seq(-3.5,3.5,by=0.2)
mat =  coefEsts
colnames(mat) = c("Fetal:Repl", "Fetal:Quies", "Adult:OPC",
	"Adult:Neuron", "Adult:Astro", "Adult:Oligo", "Adult:Microglia", "Adult:Endothelial")
rownames(mat) = rowData(rse_geneQuake[rownames(coefEsts),])$Symbol
	my.col <- colorRampPalette(c("blue","white","red"))(length(theSeq))
pdf("singleCellGroup_exprsMatZ.pdf",w=24)
print(levelplot(mat, aspect = "fill", at = theSeq,pretty=TRUE,
	panel = panel.levelplot.raster, col.regions = my.col,
		scales=list(x=list(rot=90)),	ylab = "Cell Type", xlab = ""))
dev.off()

library(pheatmap)
eMat =geneExprsQuake[rownames(coefEsts),]
colnames(eMat) = colnames(mat)[match(group, colnames(coefEsts))]
rownames(eMat) = rowData(rse_geneQuake[rownames(coefEsts),])$Symbol
pdf("singleCell_exprsMat.pdf", h=30,w=70)
print(pheatmap(eMat,
	color = colorRampPalette(brewer.pal(n = 7, name = "PuRd"))(100)))
dev.off()

### mean profile:
plot(mergeMarkerMeanExprsZ ~ group)