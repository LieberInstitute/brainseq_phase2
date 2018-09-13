####
### libraries
library(SummarizedExperiment)
library(jaffelab)
library('devtools')

## in each region, and in interaction:
## what percent of eQTLs are in the same direction and then also marginally significant

#### load in, subset, resave

###################
###### HIPPO ######
###################

# GTEx
message(paste(Sys.time(), 'loading GTEx eQTL results'))
load("eqtl_tables/mergedEqtl_output_hippo_4features.rda", verbose=TRUE)
# break up into pieces
message(paste(Sys.time(), 'breaking up into pieces'))
hippo_gtex_genes = allEqtl[which(allEqtl$Type=="Gene"),]
hippo_gtex_exons = allEqtl[which(allEqtl$Type=="Exon"),]
hippo_gtex_jxns = allEqtl[which(allEqtl$Type=="Jxn"),]
hippo_gtex_txs = allEqtl[which(allEqtl$Type=="Tx"),]
rm(allEqtl)

message(paste(Sys.time(), 'setting row names'))
rownames(hippo_gtex_genes) = paste0(hippo_gtex_genes$snps,"_",hippo_gtex_genes$gene)
rownames(hippo_gtex_exons) = paste0(hippo_gtex_exons$snps,"_",hippo_gtex_exons$gene)
rownames(hippo_gtex_jxns) = paste0(hippo_gtex_jxns$snps,"_",hippo_gtex_jxns$gene)
rownames(hippo_gtex_txs) = paste0(hippo_gtex_txs$snps,"_",hippo_gtex_txs$gene)


# BrainSeq
message(paste(Sys.time(), 'loading BrainSeq Phase II eQTL results'))
load("../eQTL_full/eqtl_tables/mergedEqtl_output_hippo_4features.rda", verbose=TRUE)
# keep only significant
h_sig = allEqtl[allEqtl$FDR < 0.01,]
rm(allEqtl)

message(paste(Sys.time(), 'setting row names'))
rownames(h_sig) = paste0(h_sig$snps,"_",h_sig$gene)

message(paste(Sys.time(), 'subsetting to significant results'))
h_sig_genes = h_sig[which(h_sig$Type=="Gene"),]
h_sig_exons = h_sig[which(h_sig$Type=="Exon"),]
h_sig_jxns = h_sig[which(h_sig$Type=="Jxn"),]
h_sig_txs = h_sig[which(h_sig$Type=="Tx"),]
rm(h_sig)

## subset GTEx to our results
message(paste(Sys.time(), 'matching gene results'))
hippo_gtex_genes = hippo_gtex_genes[rownames(i_sig_genes),]
message(paste(Sys.time(), 'saving gene results'))
save(hippo_gtex_genes,h_sig_genes, "hippo_compare_genes.rda")

message(paste(Sys.time(), 'matching exon results'))
hippo_gtex_exons = hippo_gtex_exons[rownames(i_sig_exons),]
message(paste(Sys.time(), 'saving exon results'))
save(hippo_gtex_exons,h_sig_exons, "hippo_compare_exons.rda")

message(paste(Sys.time(), 'matching jxn results'))
hippo_gtex_jxns = hippo_gtex_jxns[rownames(i_sig_jxns),]
message(paste(Sys.time(), 'saving jxn results'))
save(hippo_gtex_jxns,h_sig_jxns, "hippo_compare_jxns.rda")

message(paste(Sys.time(), 'matching tx results'))
hippo_gtex_txs = hippo_gtex_txs[rownames(i_sig_txs),]
message(paste(Sys.time(), 'saving tx results'))
save(hippo_gtex_txs,h_sig_txs, "hippo_compare_txs.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
	  