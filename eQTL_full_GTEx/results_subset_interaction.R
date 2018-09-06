####
### libraries
library(SummarizedExperiment)
library(jaffelab)
library('devtools')


## in each region, and in interaction:
## what percent of eQTLs are in the same direction and then also marginally significant


#### load in, subset, resave

#########################
###### interaction ######
#########################


# GTEx
load("eqtl_tables/mergedEqtl_output_interaction_4features.rda", verbose=TRUE)
# break up into pieces
inter_gtex_genes = allEqtl[which(allEqtl$Type=="Gene"),]
inter_gtex_exons = allEqtl[which(allEqtl$Type=="Exon"),]
inter_gtex_jxns = allEqtl[which(allEqtl$Type=="Jxn"),]
inter_gtex_txs = allEqtl[which(allEqtl$Type=="Tx"),]
rm(allEqtl)

rownames(inter_gtex_genes) = paste0(inter_gtex_genes$snps,"_",inter_gtex_genes$gene)
rownames(inter_gtex_exons) = paste0(inter_gtex_exons$snps,"_",inter_gtex_exons$gene)
rownames(inter_gtex_jxns) = paste0(inter_gtex_jxns$snps,"_",inter_gtex_jxns$gene)
rownames(inter_gtex_txs) = paste0(inter_gtex_txs$snps,"_",inter_gtex_txs$gene)


# BrainSeq
load("../eQTL_full/eqtl_tables/mergedEqtl_output_interaction_4features.rda", verbose=TRUE)
# keep only significant
i_sig = allEqtl[allEqtl$FDR < 0.01,]
rm(allEqtl)
rownames(i_sig) = paste0(i_sig$snps,"_",i_sig$gene)

i_sig_genes = i_sig[which(i_sig$Type=="Gene"),]
i_sig_exons = i_sig[which(i_sig$Type=="Exon"),]
i_sig_jxns = i_sig[which(i_sig$Type=="Jxn"),]
i_sig_txs = i_sig[which(i_sig$Type=="Tx"),]
rm(i_sig)

## subset GTEx to our results
print("Interaction gene")
inter_gtex_genes = inter_gtex_genes[rownames(i_sig_genes),]
print("Interaction exon")
inter_gtex_exons = inter_gtex_exons[rownames(i_sig_exons),]
print("Interaction jxn")
inter_gtex_jxns = inter_gtex_jxns[rownames(i_sig_jxns),]
print("Interaction tx")
inter_gtex_txs = inter_gtex_txs[rownames(i_sig_txs),]

### Save separately for easier loading
save(inter_gtex_genes,i_sig_genes, "inter_compare_genes.rda")
save(inter_gtex_exons,i_sig_exons, "inter_compare_exons.rda")
save(inter_gtex_jxns,i_sig_jxns, "inter_compare_jxns.rda")
save(inter_gtex_txs,i_sig_txs, "inter_compare_txs.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
