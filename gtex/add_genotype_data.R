###
library(jaffelab)
library(readr)
library(stringr)

# load data and filter
load("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/gtexPd.Rdata")
gtexPd$sra_accession = str_trim(gtexPd$sra_accession)
gtexPd = gtexPd[!is.na(gtexPd$sra_accession),]

## drop those w/o SRR
gtexPd = gtexPd[which(gtexPd$SAMPLE_USE == 
	"Seq_RNA_WTSS; Seq_RNA_Expression"),]

## filter to 2 brain regions
pd = gtexPd[which(gtexPd$SMTSD %in% c("Brain - Frontal Cortex (BA9)",	"Brain - Hippocampus")),]

## get snp data
bfile = "/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Genotypes/Merged/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6"
snpMapGtex = read_delim(paste0(bfile, ".bim"),delim = "\t", col_names=FALSE)
snpMapGtex = as.data.frame(snpMapGtex)
snpMapGtex$chrpos = paste0("chr", snpMapGtex$X1, ":", snpMapGtex$X4)

# match to libd
snpMap = read_delim("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6.bim",
	delim = "\t", col_names = FALSE)
snpMap = as.data.frame(snpMap)
snpMap$chrpos = paste0("chr", snpMap$X1, ":", snpMap$X4)

mm = match(snpMap$chrpos, snpMapGtex$chrpos)
snpMap = snpMap[!is.na(mm),]
snpMapGtex = snpMapGtex[mm[!is.na(mm)],]
cat(snpMapGtex$X2, file="snps_to_extract.txt", sep="\n")

## pull snps
thecall = paste0("plink --bfile ", bfile, 
	" --extract snps_to_extract.txt --recode A-transpose --out GTEX_SNPs_LIBDphase2")
system(thecall)

##### read in SNPs
genotypes = read_delim("GTEX_SNPs_LIBDphase2.traw",delim="\t")
snpMapGtex = as.data.frame(genotypes[,1:6])
snpMapGtex$CHR[snpMapGtex$CHR == "23"] = "X"
snpMapGtex$chrpos = paste0("chr", snpMapGtex$CHR, ":", snpMapGtex$POS)
snpGtex = as.data.frame(genotypes[,-(1:6)])
## reformat
colnames(snpGtex) = ss(colnames(snpGtex), "_",1)
rownames(snpMapGtex) = rownames(snpGtex) = snpMapGtex$SNP

### add MDS
mds = read.table("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Genotypes/Merged/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf05_geno10_hwe1e6.mds",
	header=TRUE, as.is=TRUE)
	mds = mds[match(ss(colnames(snpGtex), "-", 2), ss(mds$IID, "-", 2)),]
mds= mds[match(colnames(snpGtex), mds$FID),]

## filter to people
pd$SUBJID = ss(pd$SAMPID, "-",2)
mm2 = match(pd$SUBJID, ss(colnames(snpGtex),"-",2))

pdGtex = pd[!is.na(mm2),]
snpGtex = snpGtex[,mm2[!is.na(mm2)]]
mds = mds[mm2[!is.na(mm2)],]

## combine
pdGtex = cbind(pdGtex, mds[,4:13])
colnames(pdGtex)[257:266] = paste0("snpPC", 1:10)
save(snpGtex, snpMapGtex, pdGtex, 	compress=TRUE,
	file = "genotypeData_GTEx_hippoPlusDlpfc.rda")
    
## Simplify some things for later (on the eQTL code)
load('genotypeData_GTEx_hippoPlusDlpfc.rda', verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda', verbose = TRUE)

## Match by chr position
snpMap$chrpos <- paste0('chr', snpMap$CHR, ':', snpMap$POS)
m <- match(snpMapGtex$chrpos, snpMap$chrpos)
table(is.na(m))
#   FALSE
# 6827646
identical(m, 1:length(m))
# FALSE

## Load sig results
library('GenomicRanges')

## DLPFC full first
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda', verbose = TRUE)
keepSnp <- unique(allEqtl$snps[allEqtl$FDR < 0.05])
m2 <- match(keepSnp, snpMap$SNP)
stopifnot(all(!is.na(m2)))
snpMap$fdr_dlpfc <- FALSE
snpMap$fdr_dlpfc[m2] <- TRUE
table(snpMap$fdr_dlpfc)
#   FALSE    TRUE
# 4289920 2733940

## HIPPO full second
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_hippo_4features.rda', verbose = TRUE)
keepSnp <- unique(allEqtl$snps[allEqtl$FDR < 0.05])
m2 <- match(keepSnp, snpMap$SNP)
stopifnot(all(!is.na(m2)))
snpMap$fdr_hippo <- FALSE
snpMap$fdr_hippo[m2] <- TRUE
table(snpMap$fdr_hippo)
#   FALSE    TRUE
# 4733895 2289965

## Interaction next
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/matrixEqtl_output_interaction_4features.rda', verbose = TRUE)
keepSnp <- unique(c(
    as.character(meGene$cis$eqtls$snps[meGene$cis$eqtls$FDR < 0.05]),
    as.character(meExon$cis$eqtls$snps[meGene$cis$eqtls$FDR < 0.05]),
    as.character(meJxn$cis$eqtls$snps[meGene$cis$eqtls$FDR < 0.05]),
    as.character(meTx$cis$eqtls$snps[meGene$cis$eqtls$FDR < 0.05])
))
m2 <- match(keepSnp, snpMap$SNP)
stopifnot(!all(is.na(m2)))
snpMap$fdr_interaction <- FALSE
snpMap$fdr_interaction[m2] <- TRUE
table(snpMap$fdr_interaction)
#   FALSE    TRUE
# 6543153  480707
rm(meGene, meExon, meJxn, meTx)

## PGC2 rAggr DLPFC subset
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/eqtl_tables/mergedEqtl_output_dlpfc_raggr_4features.rda', verbose = TRUE)
keepSnp <- unique(allEqtl$snps[allEqtl$FDR < 0.05])
m2 <- match(keepSnp, snpMap$SNP)
stopifnot(!all(is.na(m2)))
snpMap$fdr_dlpfc_gwas <- FALSE
snpMap$fdr_dlpfc_gwas[m2] <- TRUE
table(snpMap$fdr_dlpfc_gwas)
#   FALSE    TRUE
# 7016173    7687

## finally PGC2 rAggr HIPPO subset
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/eqtl_tables/mergedEqtl_output_hippo_raggr_4features.rda', verbose = TRUE)
keepSnp <- unique(allEqtl$snps[allEqtl$FDR < 0.05])
m2 <- match(keepSnp, snpMap$SNP)
stopifnot(!all(is.na(m2)))
snpMap$fdr_hippo_gwas <- FALSE
snpMap$fdr_hippo_gwas[m2] <- TRUE
table(snpMap$fdr_hippo_gwas)
#   FALSE    TRUE
# 7016883    6977

## Across all
snpMap$fdr_any <- with(snpMap, fdr_dlpfc | fdr_hippo | fdr_interaction | fdr_dlpfc_gwas | fdr_hippo_gwas)
table(snpMap$fdr_any)
#   FALSE    TRUE
# 3755771 3268089

## Add missing info
missinginfo <- snpMap[m, colnames(snpMap)[!colnames(snpMap) %in% colnames(snpMapGtex)]]
table(is.na(missinginfo$pos_hg38))
#   FALSE    TRUE
# 6827197     449
snpMapGtex <- cbind(snpMapGtex, missinginfo)

## Dims at this point
dim(snpGtex)
dim(snpMapGtex)
dim(pdGtex)
# > dim(snpGtex)
# [1] 6827646     190
# > dim(snpMapGtex)
# [1] 6827646      21
# > dim(pdGtex)
# [1] 190 266

## Save without subsetting to FDR sig res
save(snpGtex, snpMapGtex, pdGtex, 	compress=TRUE,
	file = "genotypeData_GTEx_hippoPlusDlpfc_withFDRinfo.rda")

## drop SNPs not mapping to hg38: actually none of these are in the fdr results
table(is.na(snpMapGtex$chr_hg38), snpMapGtex$fdr_any)
keepIndex = which(snpMapGtex$fdr_any)
length(keepIndex)
# [1] 3194576
snpMapGtex = snpMapGtex[keepIndex,]
snpGtex = snpGtex[keepIndex,]

## Final dims
dim(snpGtex)
dim(snpMapGtex)
dim(pdGtex)
# > dim(snpGtex)
# [1] 3194576     190
# > dim(snpMapGtex)
# [1] 3194576      21
# > dim(pdGtex)
# [1] 190 266

## Save again
save(snpGtex, snpMapGtex, pdGtex, 	compress=TRUE,
	file = "genotypeData_GTEx_hippoPlusDlpfc_simplified.rda")


