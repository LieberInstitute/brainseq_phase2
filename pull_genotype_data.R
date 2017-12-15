########################

## read packages
library(jaffelab)
library(readr)
library(SummarizedExperiment)
library(stringr)
library(GenomicRanges)

## load data
# load("count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseGene_merged_n449.rda")
# pdDlpfc = colData(rse_gene)
# load("count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n442.rda")
# pdHippo = colData(rse_gene)

##### now with a couple added samples that were originally dropped
## load data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseGene_merged_n453.rda")
pdDlpfc = colData(rse_gene)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda")
pdHippo = colData(rse_gene)
rm(rse_gene)

## get BrNums
BrNums = unique(c(pdDlpfc$BrNum, pdHippo$BrNum))

## update Br1060
BrNums[BrNums == "Br1061"] = "Br1060"

############# 
# fam file ##

### read in fam
fam = read.table("/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10.fam")
colnames(fam) = c("BrNum", "Platform", "MID","PID","Sex","Pheno")
table(BrNums %in% fam$BrNum) # keep samples w/ genotypes

fam = fam[!duplicated(fam$BrNum),] # remove 650 if they have 1M
famOut = fam[which(fam$BrNum %in% BrNums),]
write.table(famOut[,1:2], "samples_to_extract.txt",
	col.names=FALSE, row.names=FALSE, quote=FALSE)
	
#### overall extraction
bfile = "/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10"
newbfile = "/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6"

## extract
system(paste("/users/ajaffe/bin/plink --bfile", bfile, 
	"--keep samples_to_extract.txt --geno 0.1 --maf 0.05 --hwe 0.000001 --make-bed --out", 
	newbfile, " --memory 225000"))

# ## independent and cluster
system(paste("/users/ajaffe/bin/plink --bfile", newbfile, "--indep 100 10 1.25 --out", newbfile))

## MDS components	
system(paste0("/users/ajaffe/bin/plink --bfile ", newbfile, 
	" --cluster --mds-plot 10 --extract ",
	newbfile, ".prune.in --out ", newbfile))

# ## A transpose
system(paste("/users/ajaffe/bin/plink --bfile", newbfile,
	"--recode A-transpose --out", newbfile))
	
################
## read in #####

## read in genotypes
genotypes  = read_delim(paste0(newbfile, ".traw"), delim="\t")

snp = as.data.frame(genotypes[,-(1:6)])
colnames(snp) = ss(colnames(snp), "_")
snp = as.matrix(snp[,BrNums])

### update MAP
load("/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10_updatedMap.rda")
snpMap = snpMap[match(genotypes$SNP, snpMap$SNP),]
rownames(snp) = rownames(snpMap) =  snpMap$SNP

## check directionality
swapIndex = which(snpMap$ALT == genotypes$COUNTED)
snp[swapIndex,] = 2-snp[swapIndex,]

#### read in MDS
mds = read.table(paste0(newbfile, ".mds"), 
	header=TRUE,as.is=TRUE,row.names=1)[BrNums,-(1:2)]
colnames(mds) = paste0("snpPC",1:ncol(mds))

### fix brnum
ii = which(BrNums == "Br1060")
rownames(mds)[ii] = colnames(snp)[ii]= "Br1061"

### save
save(mds, snp, snpMap, compress=TRUE,
	file = "genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")
