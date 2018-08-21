##
library(jaffelab)
library(IRanges)
library(SummarizedExperiment)


#####################
##### Subset of 881 SNPs from PGC
#####################

################
## load eQTLs

load("eqtl_tables/mergedEqtl_output_hippo_4features.rda", verbose=TRUE)
hippo = allEqtl[allEqtl$FDR < 0.01,]
load("eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)
dlp = allEqtl[allEqtl$FDR < 0.01,]
load("eqtl_tables/mergedEqtl_output_interaction_4features.rda", verbose=TRUE)
inter = allEqtl[allEqtl$FDR < 0.05,]

rm(allEqtl)

## Don't include tx eQTLs
hippo = hippo[which(hippo$Type != "Tx"),]
dlp = dlp[which(dlp$Type != "Tx"),]
inter = inter[which(inter$Type != "Tx"),]

################
## load snp map

load("../genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")
snpMap$hg19POS = paste0(snpMap$CHR,":",snpMap$POS)
# snpMap2 = snpMap[which(rownames(snpMap) %in% c(hippo$snps,dlp$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

hippo$hg19POS = snpMap$hg19POS[match(hippo$snps, snpMap$SNP)]
dlp$hg19POS = snpMap$hg19POS[match(dlp$snps, snpMap$SNP)]
inter$hg19POS = snpMap$hg19POS[match(inter$snps, snpMap$SNP)]


################
## load risk snps

riskLoci = read.csv("../eQTL_GWAS_riskSNPs/rAggr_results_179.csv", stringsAsFactors=FALSE)	# 10,981 snps
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 



##  Among the HIPPO eQTLs, we found associations to 60 risk SNPs
indexsnps = riskLoci$hg19POS1[riskLoci$hg19POS1 %in% hippo$hg19POS]
length(unique(indexsnps))
# > length(unique(indexsnps))
# [1] 60

##  Among the DLPFC eQTLs, we found associations to 65 risk SNPs
indexsnps = riskLoci$hg19POS1[riskLoci$hg19POS1 %in% dlp$hg19POS]
length(unique(indexsnps))
# > length(unique(indexsnps))
# [1] 65


##  Among the interaction eQTLs, we found associations to xx risk SNPs
indexsnps = riskLoci$hg19POS1[riskLoci$hg19POS1 %in% inter$hg19POS]
length(unique(indexsnps))
# > length(unique(indexsnps))
# [1] 3
# > unique(riskLoci$SNP1_Name[which(riskLoci$hg19POS1 %in% unique(indexsnps))])
# [1] "rs4144797:233562197:T:C"  "rs12293670:124612932:A:G" "rs324015:57490100:T:C"











