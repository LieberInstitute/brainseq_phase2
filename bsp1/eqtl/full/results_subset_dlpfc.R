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

# BSP1
message(paste(Sys.time(), 'loading BSP1 eQTL results'))
load("eqtl_tables/mergedEqtl_output_dlpfc_4features_in_progress.rda", verbose=TRUE)

message(paste(Sys.time(), 'checking for NAs on the BSP1 eQTL table'))
na_vec <- !is.na(allEqtl$snps) & !is.na(allEqtl$gene)
table(na_vec)
if(any(!na_vec)) {
    message(paste(Sys.time(), 'removing NAs from the BSP1 eQTL table'))
    allEqtl <- allEqtl[na_vec, ]
}

message(paste(Sys.time(), 'convert to a data.table'))
allEqtl <- data.table(as.data.frame(allEqtl))

# break up into pieces
message(paste(Sys.time(), 'breaking up by feature'))
dlpfc_bsp1_genes = geneEqtl
dlpfc_bsp1_exons = exonEqtl
dlpfc_bsp1_jxns = jxnEqtl
dlpfc_bsp1_txs = txEqtl

# BrainSeq
message(paste(Sys.time(), 'loading BrainSeq Phase II eQTL results'))
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)

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
d_sig_exons = proc_brainseq(d_sig[d_sig$Type=="Exon",])
d_sig_jxns = proc_brainseq(d_sig[d_sig$Type=="Jxn",])
d_sig_txs = proc_brainseq(d_sig[d_sig$Type=="Tx",])
rm(d_sig)

## subset BSP1 to our results
subset_bsp1 <- function(bsp1, brainseq) {    
    message(paste(Sys.time(), 'create keys: bsp1'))
    setkey(bsp1, snps, gene)

    message(paste(Sys.time(), 'subset bsp1 by brainseq'))
    bsp1[.(brainseq$snps, brainseq$gene)]
}

message(paste(Sys.time(), 'matching gene results'))
dlpfc_bsp1_genes <- subset_bsp1(dlpfc_bsp1_genes, d_sig_genes)
message(paste(Sys.time(), 'saving gene results'))
save(dlpfc_bsp1_genes,d_sig_genes, file = "rdas/dlpfc_compare_genes.rda")

message(paste(Sys.time(), 'matching exon results'))
dlpfc_bsp1_exons <- subset_bsp1(dlpfc_bsp1_exons, d_sig_exons)
message(paste(Sys.time(), 'saving exon results'))
save(dlpfc_bsp1_exons,d_sig_exons, file = "rdas/dlpfc_compare_exons.rda")

message(paste(Sys.time(), 'matching jxn results'))
dlpfc_bsp1_jxns <- subset_bsp1(dlpfc_bsp1_jxns, d_sig_jxns)
message(paste(Sys.time(), 'saving jxn results'))
save(dlpfc_bsp1_jxns,d_sig_jxns, file = "rdas/dlpfc_compare_jxns.rda")

message(paste(Sys.time(), 'matching tx results'))
dlpfc_bsp1_txs <- subset_bsp1(dlpfc_bsp1_txs, d_sig_txs)
message(paste(Sys.time(), 'saving tx results'))
save(dlpfc_bsp1_txs,d_sig_txs, file = "rdas/dlpfc_compare_txs.rda")


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()



	  