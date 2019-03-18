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
###### HIPPO ######
###################

# GTEx
message(paste(Sys.time(), 'loading GTEx eQTL results'))
load("eqtl_tables/mergedEqtl_output_hippo_4features.rda", verbose=TRUE)

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
hippo_gtex_genes = allEqtl[allEqtl$Type=="Gene",]
hippo_gtex_exons = allEqtl[allEqtl$Type=="Exon",]
hippo_gtex_jxns = allEqtl[allEqtl$Type=="Jxn",]
hippo_gtex_txs = allEqtl[allEqtl$Type=="Tx",]
rm(allEqtl)

# BrainSeq
message(paste(Sys.time(), 'loading BrainSeq Phase II eQTL results'))
load("../eQTL_full/eqtl_tables/mergedEqtl_output_hippo_4features.rda", verbose=TRUE)

# keep only significant
message(paste(Sys.time(), 'subsetting to significant results'))
h_sig = data.table(as.data.frame(allEqtl[allEqtl$FDR < 0.01,]))
rm(allEqtl)

message(paste(Sys.time(), 'breaking up by feature'))
proc_brainseq <- function(df) {
    message(paste(Sys.time(), 'setting keys'))
    setkey(df, snps, gene)
    return(df)
}
h_sig_genes = proc_brainseq(h_sig[h_sig$Type=="Gene",])
h_sig_exons = proc_brainseq(h_sig[h_sig$Type=="Exon",])
h_sig_jxns = proc_brainseq(h_sig[h_sig$Type=="Jxn",])
h_sig_txs = proc_brainseq(h_sig[h_sig$Type=="Tx",])
rm(h_sig)

## subset GTEx to our results
subset_gtex <- function(gtex, brainseq) {    
    message(paste(Sys.time(), 'create keys: gtex'))
    setkey(gtex, snps, gene)

    message(paste(Sys.time(), 'subset gtex by brainseq'))
    gtex[.(brainseq$snps, brainseq$gene)]
}

message(paste(Sys.time(), 'matching gene results'))
hippo_gtex_genes <- subset_gtex(hippo_gtex_genes, h_sig_genes)
message(paste(Sys.time(), 'saving gene results'))
save(hippo_gtex_genes,h_sig_genes, file = "rdas/hippo_compare_genes.rda")

message(paste(Sys.time(), 'matching exon results'))
hippo_gtex_exons <- subset_gtex(hippo_gtex_exons, h_sig_exons)
message(paste(Sys.time(), 'saving exon results'))
save(hippo_gtex_exons,h_sig_exons, file = "rdas/hippo_compare_exons.rda")

message(paste(Sys.time(), 'matching jxn results'))
hippo_gtex_jxns <- subset_gtex(hippo_gtex_jxns, h_sig_jxns)
message(paste(Sys.time(), 'saving jxn results'))
save(hippo_gtex_jxns,h_sig_jxns, file = "rdas/hippo_compare_jxns.rda")

message(paste(Sys.time(), 'matching tx results'))
hippo_gtex_txs <- subset_gtex(hippo_gtex_txs, h_sig_txs)
message(paste(Sys.time(), 'saving tx results'))
save(hippo_gtex_txs,h_sig_txs, file = "rdas/hippo_compare_txs.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
	  