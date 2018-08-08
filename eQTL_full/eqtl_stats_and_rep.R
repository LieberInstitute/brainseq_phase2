###	
library(GenomicRanges)
library(jaffelab)

## load eQTLs w/ CMC data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/allEqtls_withGtex.rda")

## change back to character
allEqtls$Type = as.character(allEqtls$Type)
allEqtls$Type = factor(allEqtls$Type, levels = c("Gene", "Exon", "Transcript","Junction","ER"))

### LD independence ###
pruned = read.table("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/plink/LIBD_Brain_DLPFC_szControls_imputed_indep.prune.in")
allEqtls$ldIndep = allEqtls$SNP %in% pruned$V1

## snp by gene
allEqtls$snpByGene = paste0(allEqtls$SNP, ".", allEqtls$EnsemblID)

## snp chrpos
snpMap_hg19= readr::read_delim("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/plink/LIBD_Brain_DLPFC_szControls_imputed.bim",
				delim = "\t", col_names=FALSE)
snpMap_hg19$chrpos = paste0("chr", snpMap_hg19$X1, ":", snpMap_hg19$X4)
allEqtls$snp_chrpos = snpMap_hg19$chrpos[match(allEqtls$SNP, snpMap_hg19$X2)]

## add biallelic
spTab = table(snpMap_hg19$chrpos)
allEqtls$biallelic = ! (allEqtls$snp_chrpos %in% names(spTab[spTab>1]))

# and hg38 for gwas catalog
load("/dcl01/lieber/ajaffe/Brain/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10_updatedMap.rda")
snpMap$chrpos38 = paste0(snpMap$chr_hg38, ":", snpMap$pos_hg38)
allEqtls$snp_chrpos38 = snpMap$chrpos38[match(allEqtls$SNP, snpMap$SNP)]
snpMap = snpMap[match(snpMap_hg19$X2, snpMap$SNP),]
snpMap[is.na(snpMap$SNP),1:6] = snpMap_hg19[is.na(snpMap$SNP),1:6]
snpMap$chrpos = paste0("chr", snpMap$CHR, ":", snpMap$POS)

## replication stuff up front
gtexTstat = allEqtls[,grep("GTEx", colnames(allEqtls))]
for(i in 1:ncol(gtexTstat)) gtexTstat[is.na(gtexTstat[,i]),i] = 0
allEqtls$GTEx_MetaZ = rowSums(as.data.frame(gtexTstat))/sqrt(ncol(gtexTstat))
allEqtls$GTEx_MetaPval = 2*pnorm(-abs(allEqtls$GTEx_MetaZ))

CMC_stat0 = allEqtls$CMC_statistic
CMC_stat0[is.na(CMC_stat0)] = 0
allEqtls$CMC_MetaZ = (allEqtls$statistic +  CMC_stat0)/sqrt(2)
allEqtls$CMC_MetaPval = 2*pnorm(-abs(allEqtls$CMC_MetaZ))

#######################
# bonf significance ###
#######################

bonfEqtls = allEqtls[allEqtls$bonf < 0.05,]
bonfEqtlsList = split(bonfEqtls, bonfEqtls$Type)

# stats
bonfTable = data.frame(snpFeaturePairs = sapply(bonfEqtlsList,nrow),
	numSnps = sapply(bonfEqtlsList, function(x) length(unique(x$SNP))),
	numFeatures = sapply(bonfEqtlsList, function(x) length(unique(x$Feature))),
	pvalCutoff = sapply(bonfEqtlsList, function(x) max(x$pvalue)),
	numEnsGene = sapply(bonfEqtlsList, function(x) length(unique(x$EnsemblID))),
	numSymGene = sapply(bonfEqtlsList, function(x) length(unique(x$Symbol))),
	medianEffect = sapply(bonfEqtlsList, function(x) signif(median(abs(x$beta)),2)),
	iqrEffect = sapply(bonfEqtlsList, function(x) 
		paste(signif(quantile(abs(x$beta),c(0.25,0.75)),2),collapse="-")),
	propNovel = sapply(bonfEqtlsList, function(x) mean(x$Class[!duplicated(x$Feature)] != "InEns")))
bonfTable

# transcript specificity
numTxBonfUnique = sapply(bonfEqtlsList, function(x) {
	x = x[!duplicated(x$Feature),]
	c(nrow(x), sum(x$NumTxGene > 1), sum(x$NumTxEqtl <= 1 & x$NumTxGene > 1))
})
numTxBonfUnique[2,]/numTxBonfUnique[1,]
txUniqueBonf = numTxBonfUnique[3,]/numTxBonfUnique[2,]

# by gene
table(unique(bonfEqtlsList$Gene$EnsemblID) %in% unique(bonfEqtlsList$Exon$EnsemblID))
table(unique(bonfEqtlsList$Exon$EnsemblID) %in% unique(bonfEqtlsList$Gene$EnsemblID))

# novel annotation
sapply(bonfEqtlsList, function(x) table(x$Class[!duplicated(x$Feature)]))

## load regions for annotation
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda")
rm(regionMat, pd)
table(regions$annoClass[names(regions) %in% bonfEqtlsList$ER$Feature])

## how many novel when no annotated
ttNovelty = table( bonfEqtls$EnsemblID, bonfEqtls$Class,bonfEqtls$Type)
ttNovelty = cbind(ttNovelty[,3,1], ttNovelty[,3,2], ttNovelty[,-4,4], ttNovelty[,-4,5])
colnames(ttNovelty) = c("Gene_InEns", "Exon_InEns", "Jxn_Alt", "Jxn_Skip",
	"Jxn_InEns", "ER_Alt", "ER_Skip", "ER_InEns")
table(rowSums(ttNovelty[,grepl("InEns", colnames(ttNovelty))] > 0), 
	rowSums(ttNovelty[,!grepl("InEns", colnames(ttNovelty))] > 0),
	dnn = c("Annotated", "Unannotated"))
	
#################
# replication ###
#################

## overall
bonfRepTable = data.frame(snpFeaturePairs = sapply(bonfEqtlsList,nrow),
	numRepPairs = sapply(bonfEqtlsList, function(x) sum(!is.na(x$CMC_beta))),
	numSnps = sapply(bonfEqtlsList, function(x) length(unique(x$SNP))),
	numRepSnps = sapply(bonfEqtlsList, function(x) length(unique(x$SNP[!is.na(x$CMC_beta)]))),
	numFeatures = sapply(bonfEqtlsList, function(x) length(unique(x$Feature))),
	numRepFeatures = sapply(bonfEqtlsList, function(x) length(unique(x$Feature[!is.na(x$CMC_beta)]))),
	propDir = sapply(bonfEqtlsList, function(x) mean(sign(x$CMC_beta) == sign(x$beta),na.rm=TRUE)),
	propRep05 = sapply(bonfEqtlsList, function(x) mean((sign(x$CMC_beta) == sign(x$beta) & 
		x$CMC_pvalue < 0.05),na.rm=TRUE)),
	propRep01 = sapply(bonfEqtlsList, function(x) mean((sign(x$CMC_beta) == sign(x$beta) & 
		x$CMC_pvalue < 0.01),na.rm=TRUE)),
	propRep001 = sapply(bonfEqtlsList, function(x) mean((sign(x$CMC_beta) == sign(x$beta) & 
		x$CMC_pvalue < 0.001),na.rm=TRUE)),
	propRep1e5 = sapply(bonfEqtlsList, function(x) mean((sign(x$CMC_beta) == sign(x$beta) & 
		x$CMC_pvalue < 0.00001),na.rm=TRUE)))
bonfRepTable		

bonfRepTable[,2] / bonfRepTable[,1] # pairs
bonfRepTable[,6] / bonfRepTable[,5] # pairs

## meta analysis
mean(bonfEqtls$CMC_MetaPval < 1e-5)
mean(bonfEqtls$CMC_MetaPval < 1e-9)

plot(bonfEqtls$statistic[!duplicated(bonfEqtls$Feature) & bonfEqtls$biallelic], 
	bonfEqtls$CMC_statistic[!duplicated(bonfEqtls$Feature) & bonfEqtls$biallelic],
		xlim = c(-60,60), ylim = c(-60,60),
		xlab = "LIBD", ylab = "CMC")

pdf("plots/suppFigure6_CMC_replication.pdf")
par(mar=c(5,6,2,2), cex.axis=1.8,cex.lab=1.8)
plot(bonfEqtls$CMC_statistic[!duplicated(bonfEqtls$Feature)  & bonfEqtls$biallelic], 
	bonfEqtls$statistic[!duplicated(bonfEqtls$Feature) & bonfEqtls$biallelic],
		xlim = c(-60,60), ylim = c(-60,60),pch=21,bg="grey",
		xlab = "CMC DLPFC", ylab = "LIBD DLPFC")
dev.off()

########### 
## GTEx ###
###########

mean(bonfEqtls$GTEx_MetaPval < 1e-5)
bonfEqtls$GTEx_Frontal_Cortex.pvalue = 2*pnorm(-abs(bonfEqtls$GTEx_Frontal_Cortex.statistic))

pdf("plots/gtex_crossRegion_plots.pdf")
par(mar=c(5,6,2,2), cex.axis=1.8,cex.lab=1.8)
plot(bonfEqtls$GTEx_Frontal_Cortex.statistic[!duplicated(bonfEqtls$Feature) & bonfEqtls$biallelic], 
	bonfEqtls$GTEx_MetaZ[!duplicated(bonfEqtls$Feature) & bonfEqtls$biallelic],
		xlim = c(-60,60), ylim = c(-60,60), pch=21,bg="grey",
		xlab = "GTEx Frontal Cortex", ylab = "GTEx Meta-Analysis")

plot(bonfEqtls$GTEx_Frontal_Cortex.statistic[!duplicated(bonfEqtls$Feature) & bonfEqtls$biallelic], 
	bonfEqtls$statistic[!duplicated(bonfEqtls$Feature) & bonfEqtls$biallelic],
		xlim = c(-60,60), ylim = c(-60,60),pch=21,bg="grey",
		xlab = "GTEx Frontal Cortex", ylab = "LIBD DLPFC")

plot(bonfEqtls$GTEx_MetaZ[!duplicated(bonfEqtls$Feature)  & bonfEqtls$biallelic], 
	bonfEqtls$statistic[!duplicated(bonfEqtls$Feature) & bonfEqtls$biallelic],
		xlim = c(-60,60), ylim = c(-60,60),pch=21,bg="grey",
		xlab = "GTEx Meta-Analysis", ylab = "LIBD DLPFC")
dev.off()		
		
table(bonfEqtls$GTEx_MetaPval < 1e-5 & bonfEqtls$GTEx_Frontal_Cortex.pvalue > 0.05)
table(allEqtls$GTEx_MetaZ > 2 & allEqtls$GTEx_Frontal_Cortex.statistic > 2)
## FDR significant

##########################
## GWAS catalog ##########
##########################

gwasStudy =read.delim("gwas_catalog_v1.0.1-studies_r2017-07-24.tsv",as.is=TRUE)
gwas =DataFrame(read.delim("gwas_catalog_v1.0.1-associations_e89_r2017-07-31.tsv",as.is=TRUE))

gwas$SNPS = CharacterList(strsplit(gwas$SNPS, "; "))
gwas$MAPPED_GENE = CharacterList(strsplit(gwas$MAPPED_GENE, "; "))
gwas$CHR_ID = CharacterList(strsplit(gwas$CHR_ID, ";"))
gwas$CHR_POS = CharacterList(strsplit(gwas$CHR_POS, ";"))

## add chr pos
gwas$chrpos = CharacterList(vector("list", nrow(gwas)))
gwas$chrpos[lengths(gwas$CHR_ID) == 1] = paste0("chr",gwas$CHR_ID[lengths(gwas$CHR_ID) == 1],
	":", gwas$CHR_POS[lengths(gwas$CHR_ID) == 1])
for(i in which(lengths(gwas$CHR_ID) > 1)) {
	gwas$chrpos[[i]] = paste0("chr",gwas$CHR_ID[[i]],":", gwas$CHR_POS[[i]])
}

# subset
allEqtlGwas = allEqtls[which(allEqtls$snpRsNum %in% unlist(gwas$SNPS) | 
	allEqtls$SNP %in% unlist(gwas$SNPS) | 
	allEqtls$snp_chrpos38 %in% unlist(gwas$chrpos)),]

## line back up to Dx
check = allEqtlGwas[!duplicated(allEqtlGwas$SNP), c("SNP","snpRsNum", "snp_chrpos38")]
gwasDx1 = data.frame(rsNum = unlist(gwas$SNPS), dx = rep(gwas$DISEASE.TRAIT,
	lengths(gwas$SNPS)), index = rep(1:nrow(gwas), lengths(gwas$SNPS)))
gwasDx2 = data.frame(chrPos = unlist(gwas$chrpos), dx = rep(gwas$DISEASE.TRAIT,
	lengths(gwas$chrpos)), index = rep(1:nrow(gwas), lengths(gwas$chrpos)))
check$gwasTrait = gwasDx1$dx[match(check$snpRsNum, gwasDx1$rsNum)]
check$gwasIndex = gwasDx1$index[match(check$snpRsNum, gwasDx1$rsNum)]
check$gwasTrait[is.na(check$gwasTrait)] = gwasDx2$dx[match(check$snp_chrpos38[is.na(check$gwasTrait)], gwasDx2$chrPos)]
check$gwasIndex[is.na(check$gwasIndex)] = gwasDx2$index[match(check$snp_chrpos38[is.na(check$gwasIndex)], gwasDx2$chrPos)]
check$gwasTrait[is.na(check$gwasTrait)] = gwasDx1$dx[match(check$SNP[is.na(check$gwasTrait)], gwasDx1$rsNum)]
check$gwasIndex[is.na(check$gwasIndex)] = gwasDx1$index[match(check$SNP[is.na(check$gwasIndex)], gwasDx1$rsNum)]
allEqtlGwas$gwasTrait = check$gwasTrait[match(allEqtlGwas$snpRsNum, check$snpRsNum)]
allEqtlGwas$gwasIndex = check$gwasIndex[match(allEqtlGwas$snpRsNum, check$snpRsNum)]
save(allEqtlGwas, gwas, file = "rdas/gwas_hits_allEqtl_subset_fdr.rda")

## enrichments
snpMap$isEqtl_fdr = snpMap$SNP %in% allEqtls$SNP
snpMap$isEqtl_bonf = snpMap$SNP %in% bonfEqtls$SNP
snpMap$isGWAS = FALSE
snpMap$isGWAS[snpMap$rsNumGuess %in% gwasDx1$rsNum] = TRUE
snpMap$isGWAS[snpMap$chrpos38 %in% gwasDx2$chrPos] = TRUE
snpMap$isGWAS[snpMap$SNP %in% gwasDx1$rsNum] = TRUE

tt_bonf = table(snpMap$isGWAS, snpMap$isEqtl_bonf, dnn = c("GWAS","eQTL"))
getOR(tt_bonf)
tt_fdr = table(snpMap$isGWAS, snpMap$isEqtl_fdr, dnn = c("GWAS","eQTL"))
getOR(tt_fdr)

###################
## sz risk SNPs ###

### load SZ GWAS clumps ###
pgc =read.delim("/users/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/pgc/daner_PGC_SCZ52_0513a.gz.p4.clump.areator.sorted.1mhc.txt",
	as.is=TRUE, header=TRUE)
tmp = unlist(strsplit(pgc$LD.friends.0.1..p0.001,","))
tmp2 = strsplit(pgc$LD.friends.0.1..p0.001,",")

theSnps = data.frame(name = c(pgc$SNP, ss(tmp, "\\(")),
	R2 = as.numeric(c(rep(1,nrow(pgc)), ss(ss(tmp, "\\(",2),"/"))),
	Dist = as.numeric(c(rep(0,nrow(pgc)), 
		gsub(")", "", ss(ss(tmp, "\\(",2),"/",2),fixed=TRUE))))

# map back
theSnps$hitIndex=c(1:nrow(pgc), rep(seq(along=tmp2), times=sapply(tmp2,length)))
theSnps$Genes = pgc$genes.6.50kb.dist2index.[theSnps$hitIndex]
theSnps$pvalue = pgc$P[theSnps$hitIndex]

# drop those that don't map, and add coordinate info
library(GenomicRanges)
load("/users/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/rdas/granges_pgc2.rda")
pgcFullGR$riskAllele = ifelse(pgcFullGR$or > 1, 
	pgcFullGR$A1, pgcFullGR$A2) # risk allele

## match up
mm = match(theSnps$name, names(pgcFullGR))
theSnps = theSnps[!is.na(mm),]
pgcStats = pgcFullGR[mm[!is.na(mm)]]
theSnps$chrpos = paste0(seqnames(pgcStats), ":", start(pgcStats))
theSnps$riskAllele = pgcStats$riskAllele

## drop those not with R^2 > 0.6
theSnps = theSnps[which(theSnps$R2 > 0.6),]

## just signif
pgcSig = read.delim("/users/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/pgc/pgc2_128loci.txt")
pgcSig = pgcSig[order(pgcSig$Rank),]
mm = match(pgcSig$Index_SNP, names(pgcFullGR))
pgcStatsSig = pgcFullGR[mm]
pgcSig$chrpos = paste0(seqnames(pgcStatsSig), 
	":", start(pgcStatsSig))
theSnps$finalHitIndex = match(theSnps$chrpos, pgcSig$chrpos)

## in our data
theSnps$inLibd = theSnps$chrpos %in% snpMap$chrpos  | 
	theSnps$name %in% snpMap$name
table(theSnps$inLibd)
theSnpsDropped = theSnps[!theSnps$inLibd,]
theSnps = theSnps[theSnps$inLibd,] # filter

#### filter our eQTLs
pgcEqtls = allEqtls[(allEqtls$snp_chrpos %in% theSnps$chrpos |
		allEqtls$snpRsNum %in% theSnps$name),]
		
## match up alleles
pgcEqtls$snpCountedAllele = snpMap$COUNTED[match(pgcEqtls$SNP, snpMap$SNP)]
pgcEqtls$snpRefAllele = snpMap$ALT[match(pgcEqtls$SNP, snpMap$SNP)]

## add stats
pgcEqtls$pgcRow = match(pgcEqtls$snp_chrpos, theSnps$chrpos)
pgcEqtls$pgcDiscRank = theSnps$hitIndex[pgcEqtls$pgcRow]
pgcEqtls$pgcFinalRank = theSnps$finalHitIndex[pgcEqtls$pgcRow]
pgcEqtls$riskAllele = theSnps$riskAllele[pgcEqtls$pgcRow]


#################
## flip slopes  #
pgcEqtls$flipSlopes = 0
pgcEqtls$flipSlopes[pgcEqtls$snpCountedAllele == pgcEqtls$riskAllele] = 1
pgcEqtls$flipSlopes[pgcEqtls$snpRefAllele == pgcEqtls$riskAllele] = -1

## CNVs
pgcEqtls$flipSlopes[nchar(pgcEqtls$snpCountedAllele) > nchar(pgcEqtls$snpRefAllele) &
	grepl("I", pgcEqtls$riskAllele)] = 1
pgcEqtls$flipSlopes[nchar(pgcEqtls$snpCountedAllele) < nchar(pgcEqtls$snpRefAllele) &
	grepl("I", pgcEqtls$riskAllele)] = -1
pgcEqtls$flipSlopes[nchar(pgcEqtls$snpCountedAllele) < nchar(pgcEqtls$snpRefAllele) &
	grepl("D", pgcEqtls$riskAllele)] = 1
pgcEqtls$flipSlopes[nchar(pgcEqtls$snpCountedAllele) > nchar(pgcEqtls$snpRefAllele) &
	grepl("D", pgcEqtls$riskAllele)] = -1
table(pgcEqtls$flipSlopes)

## look manually, 21 SNPs
as.data.frame(pgcEqtls[!duplicated(pgcEqtls$SNP) & pgcEqtls$flipSlopes==0,])
pgcEqtls = pgcEqtls[pgcEqtls$flipSlopes != 0,] # and drop

## flip
pgcEqtls$riskDir = ifelse(sign(pgcEqtls$beta*
	pgcEqtls$flipSlopes)==1,"UP","DOWN")
	
save(pgcEqtls, theSnps, theSnpsDropped, file = "rdas/PGC_SZ_hits_allEqtl_subset_fdr.rda")

## list by feature type
allEqtlsList = split(allEqtls, allEqtls$Type)

fdrTable = data.frame(snpFeaturePairs = lengths(allEqtlsList),
	numSnps = sapply(allEqtlsList, function(x) length(unique(x$SNP))),
	numFeatures = sapply(allEqtlsList, function(x) length(unique(x$Feature))),
	pvalCutoff = sapply(allEqtlsList, function(x) max(x$pvalue)),
	numEnsGene = sapply(allEqtlsList, function(x) length(unique(x$EnsemblID))),
	numSymGene = sapply(allEqtlsList, function(x) length(unique(x$Symbol))),
	medianEffect = sapply(allEqtlsList, function(x) signif(median(abs(x$beta)),2)),
	iqrEffect = sapply(allEqtlsList, function(x) 
		paste(signif(quantile(abs(x$beta),c(0.25,0.75)),2),collapse="-")),
	propNovel = sapply(allEqtlsList, function(x) mean(x$Class[!duplicated(x$Feature)] != "InEns")))
fdrTable

## gene specificity
fdrSnpByGene = table(allEqtls$snpByGene, allEqtls$Type)	
table(rowSums(fdrSnpByGene > 0))
table(rowSums(fdrSnpByGene[,c(2,4)] > 0))


