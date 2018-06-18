##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)

################
## load PGC snps
pgc = read.csv("41588_2018_59_MOESM3_ESM.csv", stringsAsFactors=FALSE)
names(pgc) = gsub("\\.\\.","_", names(pgc))
names(pgc) = gsub("\\.","_", names(pgc))


################
## load SNP data
load("../genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n546.rda", verbose=TRUE)
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

matchInd = which(pgc$Index_SNP_dbSNP_b141_ %in% snpMap$name)

pgc$snp_chr = snpMap$CHR[match(pgc$Index_SNP_dbSNP_b141_, snpMap$name)]
pgc$snp_pos_hg19 = snpMap$POS[match(pgc$Index_SNP_dbSNP_b141_, snpMap$name)]
pgc$snp_pos_hg38 = snpMap$pos_hg38[match(pgc$Index_SNP_dbSNP_b141_, snpMap$name)]

table(is.na(pgc$snp_chr))
## 45 still missing. Try loading another large dataset to get more?
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/genotype_data/zandiHyde_bipolar_Genotypes_n511.rda", verbose=TRUE)
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

pgc$snp_chr = snpMap$CHR[match(pgc$Index_SNP_dbSNP_b141_, snpMap$name)]
pgc$snp_pos_hg19 = snpMap$POS[match(pgc$Index_SNP_dbSNP_b141_, snpMap$name)]
pgc$snp_pos_hg38 = snpMap$pos_hg38[match(pgc$Index_SNP_dbSNP_b141_, snpMap$name)]

## 31 still missing
table(is.na(pgc$snp_chr))

pgc$coordsource = "manual"
pgc$coordsource[!is.na(pgc$snp_chr)] = "genotype_data"
table(pgc$coordsource)

### note 31 are still missing. Filled these in manually using https://www.ncbi.nlm.nih.gov/projects/SNP/
pgc$Gene_s_tagged = gsub(",",";",pgc$Gene_s_tagged) # remove commas to write csv
write.csv(pgc, file="pgc_riskLoci.csv", row.names=FALSE, quote=FALSE)
