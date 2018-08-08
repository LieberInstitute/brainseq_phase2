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