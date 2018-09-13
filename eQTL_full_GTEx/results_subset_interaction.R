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
message(paste(Sys.time(), 'loading GTEx eQTL results'))
load("eqtl_tables/mergedEqtl_output_interaction_4features.rda", verbose=TRUE)
# break up into pieces
message(paste(Sys.time(), 'breaking up into pieces'))
inter_gtex_genes = allEqtl[which(allEqtl$Type=="Gene"),]
inter_gtex_exons = allEqtl[which(allEqtl$Type=="Exon"),]
inter_gtex_jxns = allEqtl[which(allEqtl$Type=="Jxn"),]
inter_gtex_txs = allEqtl[which(allEqtl$Type=="Tx"),]
rm(allEqtl)

message(paste(Sys.time(), 'setting row names'))
rownames(inter_gtex_genes) = paste0(inter_gtex_genes$snps,"_",inter_gtex_genes$gene)
rownames(inter_gtex_exons) = paste0(inter_gtex_exons$snps,"_",inter_gtex_exons$gene)
rownames(inter_gtex_jxns) = paste0(inter_gtex_jxns$snps,"_",inter_gtex_jxns$gene)
## Remove NAs: otherwise it crashes due to duplicated names (there are just a handful)
# > x <- paste0(inter_gtex_txs$snps,"_",inter_gtex_txs$gene)
# > length(x)
# [1] 50185567
# > length(unique(x))
# [1] 50185558
# > y <- duplicated(x)
# > sum(y)
# [1] 9
# > x[y]
# [1] "NA_NA" "NA_NA" "NA_NA" "NA_NA" "NA_NA" "NA_NA" "NA_NA" "NA_NA" "NA_NA"
message(paste(Sys.time(), 'removing NAs from the tx GTEx eQTL table'))
table(!is.na(inter_gtex_txs$snps) & !is.na(inter_gtex_txs$gene))
inter_gtex_txs <- inter_gtex_txs[!is.na(inter_gtex_txs$snps) & !is.na(inter_gtex_txs$gene), ]
rownames(inter_gtex_txs) = paste0(inter_gtex_txs$snps,"_",inter_gtex_txs$gene)


# BrainSeq
message(paste(Sys.time(), 'loading BrainSeq Phase II eQTL results'))
load("../eQTL_full/eqtl_tables/mergedEqtl_output_interaction_4features.rda", verbose=TRUE)
# keep only significant
i_sig = allEqtl[allEqtl$FDR < 0.01,]
rm(allEqtl)

message(paste(Sys.time(), 'setting row names'))
rownames(i_sig) = paste0(i_sig$snps,"_",i_sig$gene)

message(paste(Sys.time(), 'subsetting to significant results'))
i_sig_genes = i_sig[which(i_sig$Type=="Gene"),]
i_sig_exons = i_sig[which(i_sig$Type=="Exon"),]
i_sig_jxns = i_sig[which(i_sig$Type=="Jxn"),]
i_sig_txs = i_sig[which(i_sig$Type=="Tx"),]
rm(i_sig)

system.time( test_1 <- inter_gtex_genes[rownames(i_sig_genes),] )
system.time( test_2 <- inter_gtex_genes[inter_gtex_genes$snps == i_sig_genes$snps & inter_gtex_genes$gene == i_sig_genes$gene, ] ) 

## subset GTEx to our results
message(paste(Sys.time(), 'matching gene results'))
inter_gtex_genes = inter_gtex_genes[rownames(i_sig_genes),]
message(paste(Sys.time(), 'saving gene results'))
save(inter_gtex_genes,i_sig_genes, "inter_compare_genes.rda")

message(paste(Sys.time(), 'matching exon results'))
inter_gtex_exons = inter_gtex_exons[rownames(i_sig_exons),]
message(paste(Sys.time(), 'saving exon results'))
save(inter_gtex_exons,i_sig_exons, "inter_compare_exons.rda")

message(paste(Sys.time(), 'matching jxn results'))
inter_gtex_jxns = inter_gtex_jxns[rownames(i_sig_jxns),]
message(paste(Sys.time(), 'saving jxn results'))
save(inter_gtex_jxns,i_sig_jxns, "inter_compare_jxns.rda")

message(paste(Sys.time(), 'matching tx results'))
inter_gtex_txs = inter_gtex_txs[rownames(i_sig_txs),]
message(paste(Sys.time(), 'saving tx results'))
save(inter_gtex_txs,i_sig_txs, "inter_compare_txs.rda")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
