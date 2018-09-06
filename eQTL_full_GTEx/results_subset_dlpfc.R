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
load("eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)
# break up into pieces
dlpfc_gtex_genes = allEqtl[which(allEqtl$Type=="Gene"),]
dlpfc_gtex_exons = allEqtl[which(allEqtl$Type=="Exon"),]
dlpfc_gtex_jxns = allEqtl[which(allEqtl$Type=="Jxn"),]
dlpfc_gtex_txs = allEqtl[which(allEqtl$Type=="Tx"),]
rm(allEqtl)

rownames(dlpfc_gtex_genes) = paste0(dlpfc_gtex_genes$snps,"_",dlpfc_gtex_genes$gene)
rownames(dlpfc_gtex_exons) = paste0(dlpfc_gtex_exons$snps,"_",dlpfc_gtex_exons$gene)
rownames(dlpfc_gtex_jxns) = paste0(dlpfc_gtex_jxns$snps,"_",dlpfc_gtex_jxns$gene)
rownames(dlpfc_gtex_txs) = paste0(dlpfc_gtex_txs$snps,"_",dlpfc_gtex_txs$gene)


# BrainSeq
load("../eQTL_full/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)
# keep only significant
d_sig = allEqtl[allEqtl$FDR < 0.01,]
rm(allEqtl)
rownames(d_sig) = paste0(d_sig$snps,"_",d_sig$gene)

d_sig_genes = d_sig[which(d_sig$Type=="Gene"),]
d_sig_exons = d_sig[which(d_sig$Type=="Exon"),]
d_sig_jxns = d_sig[which(d_sig$Type=="Jxn"),]
d_sig_txs = d_sig[which(d_sig$Type=="Tx"),]
rm(d_sig)

## subset GTEx to our results
print("DLPFC gene")
dlpfc_gtex_genes = dlpfc_gtex_genes[rownames(d_sig_genes),]
print("DLPFC exon")
dlpfc_gtex_exons = dlpfc_gtex_exons[rownames(d_sig_exons),]
print("DLPFC jxns")
dlpfc_gtex_jxns = dlpfc_gtex_jxns[rownames(d_sig_jxns),]
print("DLPFC tx")
dlpfc_gtex_txs = dlpfc_gtex_txs[rownames(d_sig_txs),]


### Save separately for easier loading
save(dlpfc_gtex_genes,d_sig_genes, "dlpfc_compare_genes.rda")
save(dlpfc_gtex_exons,d_sig_exons, "dlpfc_compare_exons.rda")
save(dlpfc_gtex_jxns,d_sig_jxns, "dlpfc_compare_jxns.rda")
save(dlpfc_gtex_txs,d_sig_txs, "dlpfc_compare_txs.rda")


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()



	  