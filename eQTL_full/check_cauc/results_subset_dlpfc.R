####
### libraries
library(SummarizedExperiment)
library(jaffelab)
library('data.table')
library('devtools')

dir.create('rdas', showWarnings = FALSE)

## in each region, and in interaction:
## what percent of eQTLs are in the same direction and then also marginally significant


#### load in, subset, resave

###################
###### DLPFC ######
###################

# CAUC
message(paste(Sys.time(), 'loading CAUC eQTL results'))
load("eqtl_tables/mergedEqtl_output_dlpfc_gene.rda", verbose=TRUE)

## To avoid changing the rest of the code
allEqtl <- geneEqtl
rm(geneEqtl)

message(paste(Sys.time(), 'checking for NAs on the CAUC eQTL table'))
na_vec <- !is.na(allEqtl$snps) & !is.na(allEqtl$gene)
table(na_vec)
if(any(!na_vec)) {
    message(paste(Sys.time(), 'removing NAs from the CAUC eQTL table'))
    allEqtl <- allEqtl[na_vec, ]
}

message(paste(Sys.time(), 'convert to a data.table'))
allEqtl <- data.table(as.data.frame(allEqtl))

# break up into pieces
message(paste(Sys.time(), 'breaking up by feature'))
dlpfc_cauc_genes = allEqtl[allEqtl$Type=="Gene",]
rm(allEqtl)

# BrainSeq
message(paste(Sys.time(), 'loading BrainSeq Phase II eQTL results'))
load("../eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)

# keep only significant
message(paste(Sys.time(), 'subsetting to significant results'))
d_sig = data.table(as.data.frame(allEqtl[allEqtl$FDR < 0.01,]))
rm(allEqtl)

message(paste(Sys.time(), 'breaking up by feature'))
proc_brainseq <- function(df) {
    message(paste(Sys.time(), 'setting keys'))
    setkey(df, snps, gene)
    return(df)
}
d_sig_genes = proc_brainseq(d_sig[d_sig$Type=="Gene",])
rm(d_sig)

## subset CAUC to our results
subset_cauc <- function(cauc, brainseq) {    
    message(paste(Sys.time(), 'create keys: cauc'))
    setkey(cauc, snps, gene)

    message(paste(Sys.time(), 'subset cauc by brainseq'))
    cauc[.(brainseq$snps, brainseq$gene)]
}

## Subset FDR<1% in CAUC
message(paste(Sys.time(), 'subsetting and saving FDR<1% CAUC-only results'))
d_sig_genes_cauc <- dlpfc_cauc_genes[dlpfc_cauc_genes$FDR < 0.01, ]
setkey(d_sig_genes_cauc, snps, gene)
save(d_sig_genes_cauc, file = 'rdas/d_sig_genes_cauc.Rdata')
rim(d_sig_genes_cauc)

message(paste(Sys.time(), 'matching gene results'))
dlpfc_cauc_genes <- subset_cauc(dlpfc_cauc_genes, d_sig_genes)
message(paste(Sys.time(), 'saving gene results'))
save(dlpfc_cauc_genes,d_sig_genes, file = "rdas/dlpfc_compare_genes.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()



	  