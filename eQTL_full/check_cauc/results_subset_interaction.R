####
### libraries
library(SummarizedExperiment)
library(jaffelab)
library('data.table')
library('devtools')

setDTthreads(1)

dir.create('rdas', showWarnings = FALSE)

## in each region, and in interaction:
## what percent of eQTLs are in the same direction and then also marginally significant


#### load in, subset, resave

#########################
###### interaction ######
#########################


# CAUC
message(paste(Sys.time(), 'loading CAUC eQTL results'))
load("eqtl_tables/mergedEqtl_output_interaction_gene.rda", verbose=TRUE)

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
inter_cauc_genes = allEqtl[allEqtl$Type=="Gene",]
rm(allEqtl)

# BrainSeq
message(paste(Sys.time(), 'loading BrainSeq Phase II eQTL results'))
load("../eqtl_tables/mergedEqtl_output_interaction_4features.rda", verbose=TRUE)

# keep only significant
message(paste(Sys.time(), 'subsetting to significant results'))
i_sig = data.table(as.data.frame(allEqtl[allEqtl$FDR < 0.01,]))
rm(allEqtl)

message(paste(Sys.time(), 'breaking up by feature'))
proc_brainseq <- function(df) {
    message(paste(Sys.time(), 'setting keys'))
    setkey(df, snps, gene)
    return(df)
}
i_sig_genes = proc_brainseq(i_sig[i_sig$Type=="Gene",])
rm(i_sig)

## subset CAUC to our results
subset_cauc <- function(cauc, brainseq) {    
    message(paste(Sys.time(), 'create keys: cauc'))
    setkey(cauc, snps, gene)

    message(paste(Sys.time(), 'subset cauc by brainseq'))
    cauc[.(brainseq$snps, brainseq$gene)]
}

message(paste(Sys.time(), 'matching gene results'))
inter_cauc_genes <- subset_cauc(inter_cauc_genes, i_sig_genes)
message(paste(Sys.time(), 'saving gene results'))
save(inter_cauc_genes,i_sig_genes, file = "rdas/inter_compare_genes.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
