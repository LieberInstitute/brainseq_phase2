###
library(GenomicRanges)

## load in data
load("dxStats_hippo_filtered_qSVA_geneLevel.rda")
outGene_HIPPO = outGene
outGene0_HIPPO = outGene0
load("dxStats_DLPFC_filtered_qSVA_geneLevel.rda")
outGene_DLPFC = outGene
outGene0_DLPFC = outGene0
rm(outGene, outGene0)
identical(rownames(outGene_DLPFC), rownames(outGene_HIPPO))


### compare
sum(outGene_HIPPO$adj.P.Val < 0.05)
sum(outGene_DLPFC$adj.P.Val < 0.05)
table(outGene_HIPPO$adj.P.Val < 0.05, 
	outGene_DLPFC$adj.P.Val < 0.05, dnn = c("HIPPO", "DLPFC"))
table(outGene_HIPPO$adj.P.Val < 0.1, 
	outGene_DLPFC$adj.P.Val < 0.1, dnn = c("HIPPO", "DLPFC"))
table(outGene0_HIPPO$adj.P.Val < 0.05, 
	outGene0_DLPFC$adj.P.Val < 0.05, dnn = c("HIPPO", "DLPFC"))

#########
## plots

## cross region, qSVA
pdf("HIPPO_vs_DLPFC_scatterplots_geneLevel.pdf")
par(mar=c(5,6,2,2),cex.lab=2,cex.axis=2, cex.main=1.8)
cor(outGene_HIPPO$t, outGene_DLPFC$t)
plot(outGene_HIPPO$t, outGene_DLPFC$t,
	pch = 21, bg="grey", main = "SZ T-stats (qSVA)",
	xlab="HIPPO", ylab="DLPFC")
# cor(outGene_HIPPO$logFC, outGene_DLPFC$logFC)
# plot(outGene_HIPPO$logFC, outGene_DLPFC$logFC,
	# pch = 21, bg="grey")

## cross region, no qSVA
cor(outGene0_HIPPO$t, outGene0_DLPFC$t)
plot(outGene0_HIPPO$t, outGene0_DLPFC$t,
	pch = 21, bg="grey", main = "SZ T-stats (no qSVA)",
	xlab="HIPPO", ylab="DLPFC")
dev.off()

## compare to polyA+ data
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/caseControl/rdas/expressed_de_features.rda")
geneStats_DLPFC_Old = outStatsExprs$Gene
mmEns = match(outGene_DLPFC$ensemblID, names(geneStats_DLPFC_Old))

plot(outGene_DLPFC$t, geneStats_DLPFC_Old$tstat_qsva[mmEns],
	xlab = "DLPFC RiboZero", ylab = "DLPFC polyA+",pch=21,bg="grey")

cor(outGene_DLPFC$t, geneStats_DLPFC_Old$tstat_qsva[mmEns],use="comp")
	
plot(outGene_DLPFC$t, geneStats_DLPFC_Old$CMC_tstat_qsva[mmEns])
plot(outGene0_DLPFC$t, geneStats_DLPFC_Old$tstat_adj[mmEns])
plot(outGene0_DLPFC$t, geneStats_DLPFC_Old$CMC_tstat_adj[mmEns])

## gene ontology
library(clusterProfiler)
univ = as.character(outGene_HIPPO$Entrez[!is.na(outGene_HIPPO$Entrez)])

sigGene_HIPPO = outGene_HIPPO[outGene_HIPPO$adj.P.Val < 0.05,]
table(outGene_DLPFC$P.Value[outGene_HIPPO$adj.P.Val < 0.05] < 0.05,
		sign(outGene_DLPFC$t[outGene_HIPPO$adj.P.Val < 0.05]) == 
		sign(sigGene_HIPPO$t) )
which(outGene_DLPFC$P.Value[outGene_HIPPO$adj.P.Val < 0.05] < 0.05 &
		sign(outGene_DLPFC$t[outGene_HIPPO$adj.P.Val < 0.05]) !=  
			sign(sigGene_HIPPO$t) )


table(outGene_DLPFC$adj.P.Val[outGene_HIPPO$adj.P.Val < 0.05] < 0.05,
	sign(outGene_DLPFC$t[outGene_HIPPO$adj.P.Val < 0.05]) == 
		sign(sigGene_HIPPO$t) )

plot(sigGene_HIPPO$t, outGene_DLPFC$t[outGene_HIPPO$adj.P.Val < 0.05])

geneList = split(sigGene_HIPPO$Entrez, sign(sigGene_HIPPO$logFC))
names(geneList) = c("SZ<Cont", "SZ>Cont")
geneList = lapply(geneList, function(x) as.character(x[!is.na(x)]))

goBP = compareCluster(geneList, fun = "enrichGO", 
	OrgDb = "org.Hs.eg.db", ont = "BP",universe = univ)
goMF = compareCluster(geneList, fun = "enrichGO", 
	OrgDb = "org.Hs.eg.db", ont = "MF",universe = univ)
kegg = compareCluster(geneList, fun = "enrichKEGG", 
	universe = univ)
	
## features?
load("dxStats_hippo_filtered_qSVA.rda")
length(unique(outExon$ensemblID[outExon$adj.P.Val < 0.05]))

sigExon = outExon[outExon$adj.P.Val < 0.05,]

geneList = split(sigExon$Entrez, sign(sigExon$logFC))
names(geneList) = c("SZ<Cont", "SZ>Cont")
geneList = lapply(geneList, function(x) as.character(x[!is.na(x)]))
univ = unique(as.character(outExon$Entrez[!is.na(outExon$Entrez)]))

goBP = compareCluster(geneList, fun = "enrichGO", 
	OrgDb = "org.Hs.eg.db", ont = "BP",universe = univ)
goMF = compareCluster(geneList, fun = "enrichGO", 
	OrgDb = "org.Hs.eg.db", ont = "MF",universe = univ)
kegg = compareCluster(geneList, fun = "enrichKEGG", 
	universe = univ)
	