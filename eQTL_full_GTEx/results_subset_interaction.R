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

#########################
###### interaction ######
#########################


# GTEx
message(paste(Sys.time(), 'loading GTEx eQTL results'))
load("eqtl_tables/mergedEqtl_output_interaction_4features.rda", verbose=TRUE)

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
inter_gtex_genes = allEqtl[allEqtl$Type=="Gene",]
inter_gtex_exons = allEqtl[allEqtl$Type=="Exon",]
inter_gtex_jxns = allEqtl[allEqtl$Type=="Jxn",]
inter_gtex_txs = allEqtl[allEqtl$Type=="Tx",]
rm(allEqtl)

# BrainSeq
message(paste(Sys.time(), 'loading BrainSeq Phase II eQTL results'))
load("../eQTL_full/eqtl_tables/mergedEqtl_output_interaction_4features.rda", verbose=TRUE)

# keep only significant
message(paste(Sys.time(), 'subsetting to significant results'))
i_sig = allEqtl[allEqtl$FDR < 0.01,]
rm(allEqtl)

message(paste(Sys.time(), 'breaking up by feature'))
i_sig_genes = i_sig[i_sig$Type=="Gene",]
i_sig_exons = i_sig[i_sig$Type=="Exon",]
i_sig_jxns = i_sig[i_sig$Type=="Jxn",]
i_sig_txs = i_sig[i_sig$Type=="Tx",]
rm(i_sig)

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
inter_gtex_genes <- subset_gtex(inter_gtex_genes, i_sig_genes)
message(paste(Sys.time(), 'saving gene results'))
save(inter_gtex_genes,i_sig_genes, file = "rdas/inter_compare_genes.rda")

message(paste(Sys.time(), 'matching exon results'))
inter_gtex_exons = inter_gtex_exons[rownames(i_sig_exons),]
message(paste(Sys.time(), 'saving exon results'))
save(inter_gtex_exons,i_sig_exons, file = "rdas/inter_compare_exons.rda")

message(paste(Sys.time(), 'matching jxn results'))
inter_gtex_jxns = inter_gtex_jxns[rownames(i_sig_jxns),]
message(paste(Sys.time(), 'saving jxn results'))
save(inter_gtex_jxns,i_sig_jxns, file = "rdas/inter_compare_jxns.rda")

message(paste(Sys.time(), 'matching tx results'))
inter_gtex_txs = inter_gtex_txs[rownames(i_sig_txs),]
message(paste(Sys.time(), 'saving tx results'))
save(inter_gtex_txs,i_sig_txs, file = "rdas/inter_compare_txs.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
