###
library(GenomicRanges)
library(jaffelab)

## load results
load("rdas/PGC_SZ_hits_allEqtl_subset_fdr.rda")
load("rdas/gwas_hits_allEqtl_subset_fdr.rda")

###################
### PGC first ####
###################
pgcEqtls$Symbol[is.na(pgcEqtls$EnsemblID)] = NA # some ERs

## check coverage of regions
length(unique(theSnps$finalHitIndex))
length(unique(theSnpsDropped$finalHitIndex))

## more metrics
pgcEqtlsFinal = pgcEqtls[!is.na(pgcEqtls$pgcFinalRank),]
nrow(pgcEqtlsFinal)
table(pgcEqtlsFinal$Type)

length(unique(pgcEqtlsFinal$EnsemblID)) - 1
length(unique(pgcEqtlsFinal$Symbol)) - 1

# number genes
table(lengths(sapply(split(pgcEqtlsFinal, pgcEqtlsFinal$pgcFinalRank), function(x) 
	unique(x$EnsemblID[!is.na(x$EnsemblID)]))))

# num tx
t(sapply(split(pgcEqtlsFinal, pgcEqtlsFinal$pgcFinalRank), function(x) 
	c(length(unique(x$EnsemblID[!is.na(x$EnsemblID)])), 
		length(unique(unlist(x$WhichTx[!is.na(x$WhichTx)]))))))
	
## three tiers of eQTLs
pgcEqtlList = list(FDR=pgcEqtls, 
	FDR_rep = pgcEqtls[pgcEqtls$CMC_MetaPval < 1e-8,],
	Bonf = pgcEqtls[pgcEqtls$bonf < 0.05,])
	
pgcTable = sapply(pgcEqtlList, function(x) {
	rownames(x) = NULL
	## only keep final
	x = x[!is.na(x$pgcFinalRank),]
	
	x$Class = as.character(x$Class)
	x$WhichTx[x$Type == "Gene"] = NA
	ttNovel = table(x$pgcFinalRank, x$Class)	## novelty
	ttType = table(x$pgcFinalRank, x$Type)	## type

	## summary stats
	data.frame(numGwas = length(unique(x$pgcFinalRank)),
		numGwasNoGene = sum(ttType[,"Gene"] == 0),
		numGwasNoTxOrGene = sum(rowSums(ttType[,c("Gene", "Transcript")]) == 0),
		numGwasUn = sum(rowSums(ttNovel[,-3]) > 0),
		numGwasOnlyUn = sum(ttNovel[,"InEns"] == 0),
		numGwasTx = sum(lengths(sapply(split(x$WhichTx, x$pgcFinalRank), 
			function(y) unlist(unique(y)))) <= 1))
})
pgcTable


## examples of loci
examples = pgcEqtlList$FDR_rep
examples$WhichTx[examples$Type == "Gene"] = NA
examples = examples[!is.na(examples$pgcFinalRank),]
rownames(examples) = NULL

exampleList = split(examples, examples$pgcFinalRank)

## number of genes
checkInd1 = which(lengths(sapply(exampleList, function(x) 
	unique(x$EnsemblID[!is.na(x$EnsemblID)]))) == 1)
exampleList[names(checkInd1)]
	
sapply(exampleList, function(x) unique(unlist(x$WhichTx)))

# novelty
ttNovel = table(examples$pgcFinalRank, examples$Class)	## novelty

# tx specifity 
txSpec = lengths(sapply(split(examples$WhichTx, examples$pgcFinalRank), 
			function(y) unlist(unique(y)))) <= 1

examples = examples[which(examples$pgcFinalRank %in% 
	unique(c(names(which(ttNovel[,"InEns"] == 0)), names(which(txSpec))))),] 
exampleList = split(examples,examples$pgcFinalRank)

### bring in case-control
load("/users/ajaffe/Lieber/Projects/RNAseq/SzControl_DE_paper/rdas/all_de_features.rda")
outStatsSub = unlist(endoapply(outStats, function(x) x[names(x) %in% pgcEqtls$Feature]))
outStatsSub$Type = ss(names(outStatsSub), "\\.",1)
outStatsSub$Feature = ss(names(outStatsSub), "\\.",2)
outStatsMatch = outStatsSub[match(pgcEqtls$Feature, outStatsSub$Feature)]

pgcEqtls$SzDir_qSVA = ifelse(outStatsMatch$log2FC_qsva > 0, "UP", "DOWN")
pgcEqtls$SzPval_qSVA = outStatsMatch$pval_qsva
pgcEqtls$SzDir_Adj = ifelse(outStatsMatch$log2FC_adj > 0, "UP", "DOWN")
pgcEqtls$SzPval_Adj = outStatsMatch$pval_adj

pgcEqtlsFinal = pgcEqtls[!is.na(pgcEqtls$pgcFinalRank),]

## tables
chisq.test(table(pgcEqtlsFinal$SzDir_qSVA, pgcEqtlsFinal$riskDir))

chisq.test(table(pgcEqtlsFinal$SzDir_qSVA[pgcEqtlsFinal$GTEx_MetaPval < 1e-8], 
	pgcEqtlsFinal$riskDir[pgcEqtlsFinal$GTEx_MetaPval < 1e-8]))

chisq.test(table(pgcEqtlsFinal$SzDir_qSVA[pgcEqtlsFinal$bonf < 0.05], 
	pgcEqtlsFinal$riskDir[pgcEqtlsFinal$bonf < 0.05]))

chisq.test(table(pgcEqtlsFinal$SzDir_qSVA[pgcEqtlsFinal$SzPval_qSVA < 0.05], 
	pgcEqtlsFinal$riskDir[pgcEqtlsFinal$SzPval_qSVA < 0.05]))	
	
chisq.test(table(pgcEqtlsFinal$SzDir_qSVA[pgcEqtlsFinal$SzPval_qSVA < 0.05 & pgcEqtlsFinal$bonf < 0.05], 
	pgcEqtlsFinal$riskDir[pgcEqtlsFinal$SzPval_qSVA < 0.05& pgcEqtlsFinal$bonf < 0.05]))

## by locus
sapply(split(pgcEqtlsFinal, pgcEqtlsFinal$pgcFinalRank), function(x) {
	data.frame(fdr = mean(x$riskDir == x$SzDir_qSVA),
		meta = mean(x$riskDir[x$CMC_MetaPval < 1e-8] == x$SzDir_qSVA[x$CMC_MetaPval < 1e-8]),
		bonf = mean(x$riskDir[x$bonf < 1e-8] == x$SzDir_qSVA[x$bonf < 1e-8]))
})

## by type
sapply(split(pgcEqtlsFinal, pgcEqtlsFinal$Type), function(x) {
	data.frame(fdr = mean(x$riskDir == x$SzDir_qSVA),
		meta = mean(x$riskDir[x$CMC_MetaPval < 1e-8] == x$SzDir_qSVA[x$CMC_MetaPval < 1e-8]),
		bonf = mean(x$riskDir[x$bonf < 1e-8] == x$SzDir_qSVA[x$bonf < 1e-8]))
})

pgcEqtlsFinal[pgcEqtlsFinal$SzPval_qSVA < 0.05,]
pgcEqtlsFinal[pgcEqtlsFinal$SzPval_qSVA < 0.05,]
	
## by locus

###################
### gwas second ####
###################

## three tiers of eQTLs
gwasEqtlList = list(FDR=allEqtlGwas, 
	FDR_rep = allEqtlGwas[allEqtlGwas$CMC_MetaPval < 1e-8,],
	Bonf = allEqtlGwas[allEqtlGwas$bonf < 0.05,])
	
gwasTable = sapply(gwasEqtlList, function(x) {
	x$Class = as.character(x$Class)
	x$txSpecific = x$NumTxEqtl <= 1 & x$NumTxGene > 1
	x$txSpecific[x$Type == "ER"] = FALSE
 	
	ttNovel = table(x$gwasIndex, x$Class)	## novelty
	ttSpecific = table(x$gwasIndex, x$txSpecific)	## tx specific
	ttType = table(x$gwasIndex, x$Type)	## type

	## summary stats
	data.frame(numGwas = length(unique(x$gwasIndex)),
		numGwasNoGene = sum(ttType[,"Gene"] == 0),
		numGwasNoTxOrGene = sum(rowSums(ttType[,c("Gene", "Transcript")]) == 0),
		numGwasUn = sum(rowSums(ttNovel[,-3]) > 0),
		numGwasOnlyUn = sum(ttNovel[,"InEns"] == 0),
		numGwasTx = sum(lengths(sapply(split(x$WhichTx, x$gwasIndex), 
			function(y) unlist(unique(y)))) <= 1))
})
gwasTable
