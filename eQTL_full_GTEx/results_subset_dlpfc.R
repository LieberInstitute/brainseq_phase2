####
### libraries
library(SummarizedExperiment)
library(jaffelab)
library('data.table')
library('devtools')

## in each region, and in interaction:
## what percent of eQTLs are in the same direction and then also marginally significant


#### load in, subset, resave

###################
###### DLPFC ######
###################

# GTEx
message(paste(Sys.time(), 'loading GTEx eQTL results'))
load("eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)

message(paste(Sys.time(), 'checking for NAs on the GTEx eQTL table'))
na_vec <- !is.na(allEqtl$snps) & !is.na(allEqtl$gene)
table(na_vec)
if(any(!na_vec)) {
    message(paste(Sys.time(), 'removing NAs from the GTEx eQTL table'))
    allEqtl <- allEqtl[na_vec, ]
}

message(paste(Sys.time(), 'convert to a data.table'))
allEqtl <- data.table(as.data.frame(allEqtl))

# break up into pieces
message(paste(Sys.time(), 'breaking up by feature'))
dlpfc_gtex_genes = allEqtl[allEqtl$Type=="Gene",]
dlpfc_gtex_exons = allEqtl[allEqtl$Type=="Exon",]
dlpfc_gtex_jxns = allEqtl[allEqtl$Type=="Jxn",]
dlpfc_gtex_txs = allEqtl[allEqtl$Type=="Tx",]
rm(allEqtl)

# BrainSeq
message(paste(Sys.time(), 'loading BrainSeq Phase II eQTL results'))
load("../eQTL_full/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)

# keep only significant
message(paste(Sys.time(), 'subsetting to significant results'))
d_sig = allEqtl[allEqtl$FDR < 0.01,]
rm(allEqtl)

message(paste(Sys.time(), 'breaking up by feature'))
d_sig_genes = d_sig[d_sig$Type=="Gene",]
d_sig_exons = d_sig[d_sig$Type=="Exon",]
d_sig_jxns = d_sig[d_sig$Type=="Jxn",]
d_sig_txs = d_sig[d_sig$Type=="Tx",]
rm(d_sig)

## subset GTEx to our results
subset_gtex <- function(gtex, brainseq) {
    message(paste(Sys.time(), 'converting brainseq to data.table'))
    brainseq <- data.table(as.data.frame(brainseq))
    
    message(paste(Sys.time(), 'create keys: gtex'))
    setkey(gtex, snps, gene)
    message(paste(Sys.time(), 'create keys: brainseq'))
    setkey(brainseq, snps, gene)

    message(paste(Sys.time(), 'subset gtex by brainseq'))
    gtex[.(brainseq$snps, brainseq$gene)]
}

message(paste(Sys.time(), 'matching gene results'))
dlpfc_gtex_genes <- subset_gtex(dlpfc_gtex_genes, d_sig_genes)
message(paste(Sys.time(), 'saving gene results'))
save(dlpfc_gtex_genes,d_sig_genes, file = "dlpfc_compare_genes.rda")

message(paste(Sys.time(), 'matching exon results'))
dlpfc_gtex_exons = dlpfc_gtex_exons[rownames(d_sig_exons),]
message(paste(Sys.time(), 'saving exon results'))
save(dlpfc_gtex_exons,d_sig_exons, file = "dlpfc_compare_exons.rda")

message(paste(Sys.time(), 'matching jxn results'))
dlpfc_gtex_jxns = dlpfc_gtex_jxns[rownames(d_sig_jxns),]
message(paste(Sys.time(), 'saving jxn results'))
save(dlpfc_gtex_jxns,d_sig_jxns, file = "dlpfc_compare_jxns.rda")

message(paste(Sys.time(), 'matching tx results'))
dlpfc_gtex_txs = dlpfc_gtex_txs[rownames(d_sig_txs),]
message(paste(Sys.time(), 'saving tx results'))
save(dlpfc_gtex_txs,d_sig_txs, file = "dlpfc_compare_txs.rda")


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()



	  