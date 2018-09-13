####
### libraries
library(SummarizedExperiment)
library(jaffelab)
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
# break up into pieces
message(paste(Sys.time(), 'breaking up into pieces'))
dlpfc_gtex_genes = allEqtl[which(allEqtl$Type=="Gene"),]
dlpfc_gtex_exons = allEqtl[which(allEqtl$Type=="Exon"),]
dlpfc_gtex_jxns = allEqtl[which(allEqtl$Type=="Jxn"),]
dlpfc_gtex_txs = allEqtl[which(allEqtl$Type=="Tx"),]
rm(allEqtl)

message(paste(Sys.time(), 'setting row names'))
rownames(dlpfc_gtex_genes) = paste0(dlpfc_gtex_genes$snps,"_",dlpfc_gtex_genes$gene)
rownames(dlpfc_gtex_exons) = paste0(dlpfc_gtex_exons$snps,"_",dlpfc_gtex_exons$gene)
rownames(dlpfc_gtex_jxns) = paste0(dlpfc_gtex_jxns$snps,"_",dlpfc_gtex_jxns$gene)
rownames(dlpfc_gtex_txs) = paste0(dlpfc_gtex_txs$snps,"_",dlpfc_gtex_txs$gene)


# BrainSeq
message(paste(Sys.time(), 'loading BrainSeq Phase II eQTL results'))
load("../eQTL_full/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)
# keep only significant
d_sig = allEqtl[allEqtl$FDR < 0.01,]
rm(allEqtl)

message(paste(Sys.time(), 'setting row names'))
rownames(d_sig) = paste0(d_sig$snps,"_",d_sig$gene)

message(paste(Sys.time(), 'subsetting to significant results'))
d_sig_genes = d_sig[which(d_sig$Type=="Gene"),]
d_sig_exons = d_sig[which(d_sig$Type=="Exon"),]
d_sig_jxns = d_sig[which(d_sig$Type=="Jxn"),]
d_sig_txs = d_sig[which(d_sig$Type=="Tx"),]
rm(d_sig)

## subset GTEx to our results
message(paste(Sys.time(), 'matching gene results'))
dlpfc_gtex_genes = dlpfc_gtex_genes[rownames(i_sig_genes),]
message(paste(Sys.time(), 'saving gene results'))
save(dlpfc_gtex_genes,d_sig_genes, "dlpfc_compare_genes.rda")

message(paste(Sys.time(), 'matching exon results'))
dlpfc_gtex_exons = dlpfc_gtex_exons[rownames(i_sig_exons),]
message(paste(Sys.time(), 'saving exon results'))
save(dlpfc_gtex_exons,d_sig_exons, "dlpfc_compare_exons.rda")

message(paste(Sys.time(), 'matching jxn results'))
dlpfc_gtex_jxns = dlpfc_gtex_jxns[rownames(i_sig_jxns),]
message(paste(Sys.time(), 'saving jxn results'))
save(dlpfc_gtex_jxns,d_sig_jxns, "dlpfc_compare_jxns.rda")

message(paste(Sys.time(), 'matching tx results'))
dlpfc_gtex_txs = dlpfc_gtex_txs[rownames(i_sig_txs),]
message(paste(Sys.time(), 'saving tx results'))
save(dlpfc_gtex_txs,d_sig_txs, "dlpfc_compare_txs.rda")


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()



	  