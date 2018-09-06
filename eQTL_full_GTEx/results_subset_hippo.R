####
### libraries
library(SummarizedExperiment)
library(jaffelab)


## in each region, and in interaction:
## what percent of eQTLs are in the same direction and then also marginally significant


# ###################
# ###### DLPFC ######
# ###################

# # GTEx
# load("eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)
# # break up into pieces
# dlpfc_gtex_genes = allEqtl[which(allEqtl$Type=="Gene"),]
# dlpfc_gtex_exons = allEqtl[which(allEqtl$Type=="Exon"),]
# dlpfc_gtex_jxns = allEqtl[which(allEqtl$Type=="Jxn"),]
# dlpfc_gtex_txs = allEqtl[which(allEqtl$Type=="Tx"),]

# rownames(dlpfc_gtex_genes) = paste0(dlpfc_gtex_genes$snps,"_",dlpfc_gtex_genes$gene)
# rownames(dlpfc_gtex_exons) = paste0(dlpfc_gtex_exons$snps,"_",dlpfc_gtex_exons$gene)
# rownames(dlpfc_gtex_jxns) = paste0(dlpfc_gtex_jxns$snps,"_",dlpfc_gtex_jxns$gene)
# rownames(dlpfc_gtex_txs) = paste0(dlpfc_gtex_txs$snps,"_",dlpfc_gtex_txs$gene)


# # BrainSeq
# load("../eQTL_full/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)
# # keep only significant
# d_sig = allEqtl[allEqtl$FDR < 0.01,]
# rownames(d_sig) = paste0(d_sig$snps,"_",d_sig$gene)

# d_sig_genes = d_sig[which(d_sig$Type=="Gene"),]
# d_sig_exons = d_sig[which(d_sig$Type=="Exon"),]
# d_sig_jxns = d_sig[which(d_sig$Type=="Jxn"),]
# d_sig_txs = d_sig[which(d_sig$Type=="Tx"),]

# ## subset GTEx to our results
# print("DLPFC gene")
# dlpfc_gtex_genes = dlpfc_gtex_genes[rownames(d_sig_genes),]
# print("DLPFC exon")
# dlpfc_gtex_exons = dlpfc_gtex_exons[rownames(d_sig_exons),]
# print("DLPFC jxns")
# dlpfc_gtex_jxns = dlpfc_gtex_jxns[rownames(d_sig_jxns),]
# print("DLPFC tx")
# dlpfc_gtex_txs = dlpfc_gtex_txs[rownames(d_sig_txs),]


# ### Save separately for easier loading
# save(dlpfc_gtex_genes,d_sig_genes, "dlpfc_compare_genes.rda")
# save(dlpfc_gtex_exons,d_sig_exons, "dlpfc_compare_exons.rda")
# save(dlpfc_gtex_jxns,d_sig_jxns, "dlpfc_compare_jxns.rda")
# save(dlpfc_gtex_txs,d_sig_txs, "dlpfc_compare_txs.rda")

# rm(list=ls())



###################
###### HIPPO ######
###################

# GTEx
load("eqtl_tables/mergedEqtl_output_hippo_4features.rda", verbose=TRUE)
# break up into pieces
hippo_gtex_genes = allEqtl[which(allEqtl$Type=="Gene"),]
hippo_gtex_exons = allEqtl[which(allEqtl$Type=="Exon"),]
hippo_gtex_jxns = allEqtl[which(allEqtl$Type=="Jxn"),]
hippo_gtex_txs = allEqtl[which(allEqtl$Type=="Tx"),]

rownames(hippo_gtex_genes) = paste0(hippo_gtex_genes$snps,"_",hippo_gtex_genes$gene)
rownames(hippo_gtex_exons) = paste0(hippo_gtex_exons$snps,"_",hippo_gtex_exons$gene)
rownames(hippo_gtex_jxns) = paste0(hippo_gtex_jxns$snps,"_",hippo_gtex_jxns$gene)
rownames(hippo_gtex_txs) = paste0(hippo_gtex_txs$snps,"_",hippo_gtex_txs$gene)


# BrainSeq
load("../eQTL_full/eqtl_tables/mergedEqtl_output_hippo_4features.rda", verbose=TRUE)
# keep only significant
h_sig = allEqtl[allEqtl$FDR < 0.01,]
rownames(h_sig) = paste0(h_sig$snps,"_",h_sig$gene)

h_sig_genes = h_sig[which(h_sig$Type=="Gene"),]
h_sig_exons = h_sig[which(h_sig$Type=="Exon"),]
h_sig_jxns = h_sig[which(h_sig$Type=="Jxn"),]
h_sig_txs = h_sig[which(h_sig$Type=="Tx"),]

## subset GTEx to our results
print("HIPPO gene")
hippo_gtex_genes = hippo_gtex_genes[rownames(h_sig_genes),]
print("HIPPO exon")
hippo_gtex_exons = hippo_gtex_exons[rownames(h_sig_exons),]
print("HIPPO jxns")
hippo_gtex_jxns = hippo_gtex_jxns[rownames(h_sig_jxns),]
print("HIPPO tx")
hippo_gtex_txs = hippo_gtex_txs[rownames(h_sig_txs),]


### Save separately for easier loading
save(hippo_gtex_genes,h_sig_genes, "hippo_compare_genes.rda")
save(hippo_gtex_exons,h_sig_exons, "hippo_compare_exons.rda")
save(hippo_gtex_jxns,h_sig_jxns, "hippo_compare_jxns.rda")
save(hippo_gtex_txs,h_sig_txs, "hippo_compare_txs.rda")

rm(list=ls())






	  