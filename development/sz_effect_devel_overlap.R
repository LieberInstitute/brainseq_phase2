###
library(GenomicRanges)

## load dx stats and combine
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda", verbose=TRUE)
szStats = list(DLPFC_Gene = outGene, DLPFC_Exon = outExon,
	DLPFC_Jxn = outJxn, DLPFC_Tx = outTx)
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda", verbose=TRUE)
szStats = c(szStats, list(HIPPO_Gene = outGene, HIPPO_Exon = outExon,
	HIPPO_Jxn = outJxn, HIPPO_Tx = outTx))

## check lengths	
sapply(szStats, nrow)

## load in devel stats
devObjs = list.files("rda", pattern = "_only",ful=TRUE)
devObjs = devObjs[c(3,1,5,7,4,2,6,8)]
names(devObjs) = names(szStats)
devStats = lapply(devObjs, function(x) {
	cat(".")
	load(x)
	return(top)
})

## add bonf
devStats = lapply(devStats, function(x) {
	x$p_bonf = p.adjust(x$P.Value, method = "bonf")
	return(x)
})
sapply(devStats, function(x) mean(x$p_bonf < 0.05))

## line up
szStats = mapply(function(s, d) s[rownames(d),], 
	szStats, devStats, SIMPLIFY=FALSE)

## add devel stats
szStats = mapply(function(s, d) {
	s$devReg = d$adj.P.Val < 0.05
	s$negCorr = d$adj.P.Val < 0.05 & d$ageCorr < 0
	s$posCorr = d$adj.P.Val < 0.05 & d$ageCorr > 0
	return(s)
}, szStats, devStats, SIMPLIFY=FALSE)

### find overlaps by features
oStats = lapply(szStats, function(x) {
	o = rbind(summary(lm(x$t ~ x$negCorr))$coef[2,c(1,4)],
		summary(lm(x$t ~ x$posCorr))$coef[2,c(1,4)])
	rownames(o) = c("fetal", "postnatal")
	return(o)
})

out = t(sapply(oStats, function(x) c(x[,1], x[,2])))
colnames(out) = paste0(colnames(out), "_", rep(c("shift", "pval"),each=2))
out = out[c(1,5,2,6,3,7,4,8),] # put features together

dir.create("table")
write.csv(out, file = "table/szEffect_develEnrich.csv")

## all plots
pdf("pdf/densityPlots_dxEffects_devRegByDir.pdf",w=10,h=4)
par(mar=c(5,6,2,2), mfrow = c(1,2), cex.axis=2,cex.lab=2)
for(i in seq(along=szStats)) {
	x = szStats[[i]]
		plot(density(x$t[x$posCorr==0]),
		col="grey",lwd=2,main=names(szStats)[i],
		xlab="SZ vs Control",xlim=c(-7,7), ylim = c(0, 0.4))
	lines(density(x$t[x$posCorr==1]),
		col="darkorange",lwd=3)
	legend("topright", paste0("p=", signif(oStats[[i]][2,2],3)),
		bty="n", cex=1.5)
		
	plot(density(x$t[x$negCorr==0]),
		col="grey",lwd=2,main=names(szStats)[i],
		xlab="SZ vs Control",xlim=c(-7,7), 	ylim = c(0, 0.4))
	lines(density(x$t[x$negCorr==1]),
		col="darkblue",lwd=3)
	legend("topright", paste0("p=", signif(oStats[[i]][1,2],3)),
		bty="n", cex=1.5)
		

}
dev.off()
