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

## keep those in LIBD by coordinates
snpMapGtex = snpMapGtex[snpMapGtex$chrpos %in% snpMap$chrpos,]
dim(snpMapGtex)

## stop and check
stopifnot(all(snpMapGtex$chrpos %in% snpMap$chrpos))

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

snpMap$X1[snpMap$X1 == "23"] = "X"
snpMap$chrpos = paste0("chr", snpMap$X1, ":", snpMap$X4)

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

## add LIBD SNP id
ttPos = table(snpMapGtex$chrpos)
table(snpMapGtex$chrpos %in% snpMap$chrpos)

snpMapGtex$chrpos_ref_count = paste(snpMapGtex$chrpos, snpMapGtex$ALT, snpMapGtex$COUNTED,sep="_") 
snpMap$chrpos_ref_count = paste(snpMap$chrpos, snpMap$X6, snpMap$X5, sep="_")
snpMap$chrpos_count_ref = paste(snpMap$chrpos, snpMap$X5, snpMap$X6, sep="_")
snpMapGtex$numPos=ttPos[snpMapGtex$chrpos]

snpMapGtex$mmGtexToLibd = match(snpMapGtex$chrpos_ref_count, snpMap$chrpos_ref_count)
snpMapGtex$mmGtexToLibdFlip = match(snpMapGtex$chrpos_ref_count, snpMap$chrpos_count_ref)

## check things out
table(is.na(snpMapGtex$mmGtexToLibd), snpMapGtex$numPos)
table(is.na(snpMapGtex$mmGtexToLibdFlip), snpMapGtex$numPos)
table(is.na(snpMapGtex$mmGtexToLibdFlip) & is.na(snpMapGtex$mmGtexToLibd), snpMapGtex$numPos)

## only keep those that also match on alleles (plus chrpos)
snpGtex = snpGtex[!(is.na(snpMapGtex$mmGtexToLibdFlip) & is.na(snpMapGtex$mmGtexToLibd)),]
snpMapGtex = snpMapGtex[!(is.na(snpMapGtex$mmGtexToLibdFlip) & is.na(snpMapGtex$mmGtexToLibd)),]
nrow(snpGtex)
identical(nrow(snpGtex), nrow(snpMapGtex))

## flip some alleles
flipIndex = which(!is.na(snpMapGtex$mmGtexToLibdFlip))
snpGtex[flipIndex,] = 2 - snpGtex[flipIndex,]
snpMapGtex$ALT[flipIndex] = snpMap$X5[snpMapGtex$mmGtexToLibdFlip[flipIndex]]
snpMapGtex$COUNTED[flipIndex] = snpMap$X6[snpMapGtex$mmGtexToLibdFlip[flipIndex]]
snpMapGtex$mmGtexToLibd[is.na(snpMapGtex$mmGtexToLibd)] = snpMapGtex$mmGtexToLibdFlip[is.na(snpMapGtex$mmGtexToLibd)]
table(is.na(snpMapGtex$mmGtexToLibd))
snpMapGtex$SNP = snpMap$X2[snpMapGtex$mmGtexToLibd]
snpMapGtex$OLDSNP = rownames(snpMapGtex)

save(snpGtex, snpMapGtex, pdGtex, 	compress=TRUE,
	file = "genotypeData_GTEx_hippoPlusDlpfc_LIBDmatched.rda")
    
## Simplify some things for later (on the eQTL code)
load('genotypeData_GTEx_hippoPlusDlpfc_LIBDmatched.rda', verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda', verbose = TRUE)

## Match by name
m <- match(snpMapGtex$SNP, snpMap$SNP)
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
    as.character(meExon$cis$eqtls$snps[meExon$cis$eqtls$FDR < 0.05]),
    as.character(meJxn$cis$eqtls$snps[meJxn$cis$eqtls$FDR < 0.05]),
    as.character(meTx$cis$eqtls$snps[meTx$cis$eqtls$FDR < 0.05])
))
m2 <- match(keepSnp, snpMap$SNP)
stopifnot(!all(is.na(m2)))
snpMap$fdr_interaction <- FALSE
snpMap$fdr_interaction[m2] <- TRUE
table(snpMap$fdr_interaction)
#   FALSE    TRUE
# 6889507  134353
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
# 3826525 3197335

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
# [1] 3125283     190
# > dim(snpMapGtex)
# [1] 3125283      21
# > dim(pdGtex)
# [1] 190 266
# > dim(snpMapGtex) # aj 9/7 update
# [1] 3135500      26
## Save again
save(snpGtex, snpMapGtex, pdGtex, 	compress=TRUE,
	file = "genotypeData_GTEx_hippoPlusDlpfc_simplified.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
devtools::session_info()

# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.5.0 Patched (2018-04-30 r74679)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       US/Eastern
#  date     2018-08-30
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package          * version   date       source
#  assertthat         0.2.0     2017-04-11 CRAN (R 3.5.0)
#  base             * 3.5.0     2018-05-02 local
#  bindr              0.1.1     2018-03-13 CRAN (R 3.5.0)
#  bindrcpp           0.2.2     2018-03-29 CRAN (R 3.5.0)
#  BiocGenerics     * 0.26.0    2018-05-03 Bioconductor
#  bitops             1.0-6     2013-08-17 CRAN (R 3.5.0)
#  colorout         * 1.2-0     2018-05-02 Github (jalvesaq/colorout@c42088d)
#  colorspace         1.3-2     2016-12-14 CRAN (R 3.5.0)
#  compiler           3.5.0     2018-05-02 local
#  crayon             1.3.4     2017-09-16 CRAN (R 3.5.0)
#  datasets         * 3.5.0     2018-05-02 local
#  devtools           1.13.6    2018-06-27 CRAN (R 3.5.0)
#  digest             0.6.15    2018-01-28 CRAN (R 3.5.0)
#  dplyr              0.7.6     2018-06-29 CRAN (R 3.5.0)
#  GenomeInfoDb     * 1.16.0    2018-05-03 Bioconductor
#  GenomeInfoDbData   1.1.0     2018-04-17 Bioconductor
#  GenomicRanges    * 1.32.6    2018-07-20 Bioconductor
#  ggplot2            3.0.0     2018-07-03 CRAN (R 3.5.0)
#  glue               1.3.0     2018-07-17 CRAN (R 3.5.0)
#  graphics         * 3.5.0     2018-05-02 local
#  grDevices        * 3.5.0     2018-05-02 local
#  grid               3.5.0     2018-05-02 local
#  gtable             0.2.0     2016-02-26 CRAN (R 3.5.0)
#  htmltools          0.3.6     2017-04-28 CRAN (R 3.5.0)
#  htmlwidgets        1.2       2018-04-19 CRAN (R 3.5.0)
#  httpuv             1.4.5     2018-07-19 CRAN (R 3.5.0)
#  IRanges          * 2.14.10   2018-05-17 Bioconductor
#  later              0.7.3     2018-06-08 CRAN (R 3.5.0)
#  lattice            0.20-35   2017-03-25 CRAN (R 3.5.0)
#  lazyeval           0.2.1     2017-10-29 CRAN (R 3.5.0)
#  magrittr           1.5       2014-11-22 CRAN (R 3.5.0)
#  memoise            1.1.0     2017-04-21 CRAN (R 3.5.0)
#  methods          * 3.5.0     2018-05-02 local
#  munsell            0.5.0     2018-06-12 CRAN (R 3.5.0)
#  parallel         * 3.5.0     2018-05-02 local
#  pillar             1.3.0     2018-07-14 CRAN (R 3.5.0)
#  pkgconfig          2.0.1     2017-03-21 CRAN (R 3.5.0)
#  plyr               1.8.4     2016-06-08 CRAN (R 3.5.0)
#  png                0.1-7     2013-12-03 CRAN (R 3.5.0)
#  promises           1.0.1     2018-04-13 CRAN (R 3.5.0)
#  purrr              0.2.5     2018-05-29 CRAN (R 3.5.0)
#  R6                 2.2.2     2017-06-17 CRAN (R 3.5.0)
#  Rcpp               0.12.18   2018-07-23 CRAN (R 3.5.0)
#  RCurl              1.95-4.11 2018-07-15 CRAN (R 3.5.0)
#  rlang              0.2.1     2018-05-30 cran (@0.2.1)
#  rmote            * 0.3.4     2018-05-02 deltarho (R 3.5.0)
#  S4Vectors        * 0.18.3    2018-06-13 Bioconductor
#  scales             0.5.0     2017-08-24 CRAN (R 3.5.0)
#  servr              0.10      2018-05-30 CRAN (R 3.5.0)
#  stats            * 3.5.0     2018-05-02 local
#  stats4           * 3.5.0     2018-05-02 local
#  tibble             1.4.2     2018-01-22 CRAN (R 3.5.0)
#  tidyselect         0.2.4     2018-02-26 CRAN (R 3.5.0)
#  tools              3.5.0     2018-05-02 local
#  utils            * 3.5.0     2018-05-02 local
#  withr              2.1.2     2018-03-15 CRAN (R 3.5.0)
#  xfun               0.3       2018-07-06 CRAN (R 3.5.0)
#  XVector            0.20.0    2018-05-03 Bioconductor
#  zlibbioc           1.26.0    2018-05-02 Bioconductor
