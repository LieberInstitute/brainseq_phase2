## Based on /dcl01/lieber/ajaffe/lab/brainseq_phase2/get_degradation_regions.R
###

library('jaffelab')
library('rtracklayer')
library('recount.bwtool')
library('BiocParallel')
library('SummarizedExperiment')
library('devtools')

#####################################
### DLPFC ########################
############################

## read manifest
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseGene_merged_n449.rda")
pdDlpfc = colData(rse_gene)
	
## import degradation regions
bed = import("/dcl01/lieber/ajaffe/lab/degradation_experiments/DLPFC_RiboZero/bed/DLPFC_RiboZero_degradation_regions_bonf.bed")

## Drop 'score' otherwise recount.bwtool fails with some 'integer' error
mcols(bed) <- DataFrame('name' = bed$name)

## designate bigwigs
forwardBw = paste0("/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Coverage/",
	sapply(pdDlpfc$SAMPLE_ID,"[",1),".Forward.bw")
reverseBw = paste0("/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Coverage/",
	sapply(pdDlpfc$SAMPLE_ID,"[",1), ".Reverse.bw")
stopifnot(all(file.exists(c(forwardBw,reverseBw)))) # TRUE
names(forwardBw) = names(reverseBw) = sapply(pdDlpfc$SAMPLE_ID,"[",1)

## try coverage tool
covForward = coverage_bwtool(forwardBw, bed, strand = "+", 
	sumsdir = "degradation", bpparam = MulticoreParam(8))
covForward$bigwig_path = NULL
covForward$bigwig_file = NULL

covReverse = coverage_bwtool(reverseBw, bed, strand = "-", 
	sumsdir = "degradation", bpparam = MulticoreParam(8))
covReverse$bigwig_path = NULL
covReverse$bigwig_file = NULL

## combine
cov_rse = rbind(covForward, covReverse)	
rownames(cov_rse) = rowData(cov_rse)$name
cov_rse = cov_rse[bed$name,]

## divide by number of reads
assays(cov_rse)$counts = assays(cov_rse)$counts/100 # divide by read length

## make positive
assays(cov_rse)$counts = abs(assays(cov_rse)$counts) 

## add which ones are bonf
#bedBonf= import("/dcl01/lieber/ajaffe/lab/degradation_experiments/DLPFC_RiboZero/bed/DLPFC_RiboZero_degradation_regions_bonf.bed")
#rowData(cov_rse)$bonfSig = FALSE
#rowData(cov_rse)$bonfSig[rownames(cov_rse) %in% bedBonf$name] = TRUE

## The above piece is not needed since they are all bonf sig for DLPFC
rowData(cov_rse)$bonfSig <- TRUE

## save to final people
cov_rse_dlpfc = cov_rse
save(cov_rse_dlpfc, file = "/dcl01/lieber/ajaffe/lab/brainseq_phase2/count_data/degradation_rse_phase2_dlpfc.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
