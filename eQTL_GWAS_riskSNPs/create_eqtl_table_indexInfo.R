##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)
library(VennDiagram)

################
## load SNP data
load("../genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)
## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap = snpMap[-snpInd,]

## risk loci from PGC paper
indexLoci = read.csv("pgc_riskLoci.csv", stringsAsFactors=FALSE) ## 179
indexLoci$hg19POS = paste0(indexLoci$Chromosome, ":", indexLoci$snp_pos_hg19)
indexIndex = which(snpMap$pos_hg19 %in% indexLoci$hg19POS)	# keep 135

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("rAggr_results_179.csv", stringsAsFactors=FALSE)	# 10,981 snps
riskLoci_full = riskLoci
colnames(riskLoci) = colnames(riskLoci_full) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## keep SNPs from list
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS2)	# keep 9735 snps from snpMap
snpMap = snpMap[keepIndex,]
keepIndex = which(riskLoci$hg19POS2 %in% snpMap$pos_hg19)	# keep 9698 snps from riskLoci
riskLoci = riskLoci[keepIndex,]

snpMap$Status = ifelse(snpMap$pos_hg19 %in% indexLoci$hg19POS, "Index","Proxy")
riskLoci$Status1 = ifelse(riskLoci$hg19POS1 %in% indexLoci$hg19POS, "Index","Proxy")
riskLoci$Status2 = ifelse(riskLoci$hg19POS2 %in% indexLoci$hg19POS, "Index","Proxy")

## Also keep track of full list before dropping
riskLoci_full$hg19POS1 = paste0(riskLoci_full$SNP1_Chr, ":", riskLoci_full$SNP1_Pos) 
riskLoci_full$hg19POS2 = paste0(riskLoci_full$SNP2_Chr, ":", riskLoci_full$SNP2_Pos) 
riskLoci_full$SNP2_missing = "missing"
riskLoci_full$SNP2_missing[keepIndex] = "analyzed"
riskLoci_full$Status1 = ifelse(riskLoci_full$hg19POS1 %in% indexLoci$hg19POS, "Index","Proxy")
riskLoci_full$Status2 = ifelse(riskLoci_full$hg19POS2 %in% indexLoci$hg19POS, "Index","Proxy")


### numbers:
# 881 indexSNPs
# 456 indexSNPs in our SNPs used in eQTLS				#  length(unique(riskLoci$hg19POS1)) 490?
# 13,592 rAggr riskLoci (including index SNPs)  		#  nrow(riskLoci_full)
# 10,777 rAggr riskLoci in our SNPs used in eQTLS 		#  nrow(snpMap)


################
## load table
dlp = read.csv("raggr179snps_dlp_eqtls_fdr01.csv", row.names=1)
hippo = read.csv("raggr179snps_hippo_eqtls_fdr01.csv", row.names=1)

## unique SNPs
length(unique(sacc$hg19POS))
length(unique(amyg$hg19POS))
length(unique(dlp$hg19POS))

venn.diagram(list(Amygdala = amyg$hg19POS, sACC = sacc$hg19POS, DLPFC = dlp$hg19POS), 
	fill = c("slateblue", "skyblue3", "lightcyan2"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_SNP.png")

## unique index SNPs
length(unique(sacc$hg19POS[sacc$hg19POS %in% riskLoci$hg19POS1]))
length(unique(amyg$hg19POS[amyg$hg19POS %in% riskLoci$hg19POS1]))
length(unique(dlp$hg19POS[dlp$hg19POS %in% riskLoci$hg19POS1]))

venn.diagram(list(Amygdala = unique(amyg$hg19POS[amyg$hg19POS %in% riskLoci$hg19POS1]), 
				sACC = unique(sacc$hg19POS[sacc$hg19POS %in% riskLoci$hg19POS1]), 
				DLPFC = unique(dlp$hg19POS[dlp$hg19POS %in% riskLoci$hg19POS1])), 
	fill = c("slateblue", "skyblue3", "lightcyan2"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_index.png")

	
## unique features
tapply(sacc$gene, sacc$Type, function(x) length(unique(x)))
tapply(amyg$gene, amyg$Type, function(x) length(unique(x)))
tapply(dlp$gene, dlp$Type, function(x) length(unique(x)))

venn.diagram(list(Amygdala = amyg$gene, sACC = sacc$gene, DLPFC = dlp$gene), 
	fill = c("slateblue", "skyblue3", "lightcyan2"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_features.png")

## SNP-feature pairs
nrow(sacc)  ## 38609
nrow(amyg)   ## 22572
nrow(dlp)   ## 14669
table(sacc$Type)
table(amyg$Type)
table(dlp$Type)

venn.diagram(list(Amygdala = paste0(amyg$SNP, amyg$gene), 
				sACC = paste0(sacc$SNP, sacc$gene), 
				DLPFC = paste0(dlp$SNP, dlp$gene) ), 
	fill = c("slateblue", "skyblue3", "lightcyan2"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_snpfeatpairs.png")



## Unique symbols in SNP-feature pairs
length(unique(sacc$Symbol))
length(unique(amyg$Symbol))
length(unique(dlp$Symbol))
tapply(sacc$Symbol, sacc$Type, function(x) length(unique(x)))
tapply(amyg$Symbol, amyg$Type, function(x) length(unique(x)))
tapply(dlp$Symbol, dlp$Type, function(x) length(unique(x)))

venn.diagram(list(Amygdala = amyg$Symbol, sACC = sacc$Symbol, DLPFC = dlp$Symbol), 
	fill = c("slateblue", "skyblue3", "lightcyan2"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_symbol.png")

	
	
	
	
################################################################
###### Index SNP info ##########
################################################################

## Sort riskLoci by R2 so highest linked are chosen (i.e. index matches with itself)
riskLoci = riskLoci[order(riskLoci$Status2, decreasing=FALSE),]
riskLoci = riskLoci[order(riskLoci$R_squared, decreasing=TRUE),]

region = hippo

## note which proxy snps have a significant result
riskLoci$proxy_FDRsig = "na"
for (i in 1:nrow(riskLoci)) {
	pos = riskLoci$hg19POS2[i]
	sig = which(region$hg19POS == pos)
	riskLoci$proxy_FDRsig[i] = ifelse(length(sig) > 0, TRUE, FALSE)
}	

## Is SNP index snp or proxy
region$Status = ifelse(region$hg19POS %in% riskLoci$hg19POS1, "Index", "Proxy")

## What is the index snp for each row
proxInd = match(region$hg19POS, riskLoci$hg19POS2)
region$IndexSNP = riskLoci$SNP1_Name[proxInd]
region$IndexSNP_hg19POS = riskLoci$hg19POS1[proxInd]

## was index snp checked in eqtl analysis at all
indexInd = match(region$IndexSNP_hg19POS, riskLoci_full$hg19POS2) ## row of proxy
region$IndexSNP_indata = riskLoci_full$SNP2_missing[indexInd]

## does index snp have any significant eqtl result
region$IndexSNP_fdrSig = "na"
for (i in 1:nrow(region)) {
	if (region$IndexSNP_indata[i] == "analyzed") {
		pos = region$IndexSNP_hg19POS[i]
		sig = which(region$hg19POS == pos)
		region$IndexSNP_fdrSig[i] = ifelse(length(sig) > 0, TRUE, FALSE)
	}
}
## what is the most significant eqtl result for the index snp
region$IndexSNP_mostSigFeat = NA
region$IndexSNP_mostSigFeat_gene = NA
for (i in 1:nrow(region)) {
	if (region$IndexSNP_fdrSig[i] == "TRUE") {
		pos = region$IndexSNP_hg19POS[i]
		tmp = region[which(region$hg19POS == pos),]
		tmp = tmp[order(tmp$pvalue, decreasing=FALSE),]
		region$IndexSNP_mostSigFeat[i] = as.character(tmp$gene[1])
		region$IndexSNP_mostSigFeat_gene[i] = as.character(tmp$Symbol[1])
	}
}

## what is the lead variant for each SNP -- index if sig, else highest LD proxy
region$leadVariant = NA
for (i in 1:nrow(region)) {
	if (region$IndexSNP_fdrSig[i] == "TRUE") {  ## index snp is fdr significant
		region$leadVariant[i] = region$IndexSNP[i]
		
	} else {									## index snp is NOT fdr significant
		## find highest LD proxy snp that is sig
		pos = region$IndexSNP_hg19POS[i]
		tmp = riskLoci[which(riskLoci$hg19POS1 == pos),]
		t_ind = which(tmp$Distance==0 | tmp$proxy_FDRsig==FALSE)
		if (length(t_ind>0)) { tmp = tmp[-t_ind,] }
		tmp = tmp[order(tmp$R_squared, rev(tmp$Distance), decreasing=TRUE),]
		region$leadVariant[i] = as.character(tmp$SNP2_Name[1])
		}
}
## Is the SNP the lead variant? (Check by position)
leadVarInd = match(region$leadVariant, riskLoci$SNP2_Name)
leadVarPos = riskLoci$hg19POS1[leadVarInd]
region$leadVariant_indicator = (region$hg19POS == leadVarPos)


hippo = region


write.csv(dlp, file="raggr_179_snps_dlp_eqtls_fdr01.csv")
write.csv(hippo, file="raggr_179_snps_hippo_eqtls_fdr01.csv")











