###
library(limma)
library(GenomicRanges)
library(jaffelab)
library(RColorBrewer)

## load data
load("expressedRegions_DLPFC_RiboZero_degradation_cut5_hg38_n20.rda")

## filter
keepIndex = which(width(regions) > 50 & 
	regions$annoClass %in% c("strictExonic","exonIntron"))
regionMat = regionMat[keepIndex,]
regions = regions[keepIndex]

## model degradation time
yCov = log2(regionMat + 1)
mod = model.matrix(~DegradationTime + factor(BrNum), data=pd)
fit = lmFit(yCov, mod)
eb = ebayes(fit)

## extract
stats = regions
stats$log2FC_min = fit$coef[,2]
stats$tstat = eb$t[,2]
stats$pvalue = eb$p[,2]
stats$fdr = p.adjust(stats$pvalue, "fdr")
stats$bonf = p.adjust(stats$pvalue, "bonf")
stats$regionNum = seq(along=stats)

## filter to significance
sig = stats[order(stats$pvalue)[1:500]]

## plots
yCovSig = yCov[names(sig),]
bIndexes = splitit(pd$BrNum)

pdf("degradation_region_plots.pdf")
palette(brewer.pal(6,"Dark2"))
par(mar=c(5,6,3,2),	cex.axis=2, cex.lab=2, cex.main=2)
for( i in 1:500) {
	plot(yCovSig[i,] ~ pd$DegradationTime,cex=3,
		pch = 21, bg= as.numeric(factor(pd$BrNum)),
		xlab="Degradation Time (min)",		ylab="log2 Expression", 
		main = paste(sig$nearestSymbol[i],"-", sig$annoClass[i]))
	for(j in seq(along=bIndexes)) {
		jj= bIndexes[[j]]
		lines(yCovSig[i,]  ~ pd$DegradationTime, subset=jj, col=j,lwd=8)
	}
}
dev.off()