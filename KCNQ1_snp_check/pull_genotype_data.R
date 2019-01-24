########################

## read packages
library(jaffelab)
library(readr)
library(SummarizedExperiment)
library(stringr)
library(GenomicRanges)


#### same 551 samples as main analysis, only 1 snp


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
fam = read.table("imputedplink/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_chr11.imputed.fam")
colnames(fam) = c("BrNum", "Platform", "MID","PID","Sex","Pheno")
table(BrNums %in% fam$BrNum) # keep samples w/ genotypes

fam = fam[!duplicated(fam$BrNum),] # remove 650 if they have 1M
famOut = fam[which(fam$BrNum %in% BrNums),]
write.table(famOut[,1:2], "samples_to_extract.txt",
	col.names=FALSE, row.names=FALSE, quote=FALSE)
	
#### overall extraction
bfile = "imputedplink/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_chr11.imputed"
newbfile = "BrainSeq_Phase2_chr11_n551_rs34097980"

## extract
system(paste("/users/ajaffe/bin/plink --bfile", bfile, 
	"--keep samples_to_extract.txt --maf 0.0005 --hwe 0.000001 --make-bed --out", 
	newbfile, " --memory 225000"))

# # ## independent and cluster
# system(paste("/users/ajaffe/bin/plink --bfile", newbfile, "--indep 100 10 1.25 --out", newbfile))

# ## MDS components	
# system(paste0("/users/ajaffe/bin/plink --bfile ", newbfile, 
	# " --cluster --mds-plot 10 --extract ",
	# newbfile, ".prune.in --out ", newbfile))

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


##
snpMap = as.data.frame(genotypes[,(1:6)])
snpMap$chr_hg38 = snpMap$CHR

snpInd = grep("rs34097980", genotypes$SNP)


### fix brnum
ii = which(BrNums == "Br1060")
colnames(snp)[ii]= "Br1061"


### save
save(snp, snpMap, snpInd, compress=TRUE,
	file = "BrainSeq_Phase2_chr11_n551.rda")


