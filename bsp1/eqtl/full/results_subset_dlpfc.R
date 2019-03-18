####
### libraries
library(SummarizedExperiment)
library(jaffelab)
library('data.table')
library('devtools')

setDTthreads(threads = 1)

dir.create('rdas', showWarnings = FALSE)

## in each region, and in interaction:
## what percent of eQTLs are in the same direction and then also marginally significant


#### load in, subset, resave

###################
###### DLPFC ######
###################

# BSP1
message(paste(Sys.time(), 'loading BSP1 eQTL results'))
bsp1 <- fread(
    '/dcl01/lieber/ajaffe/lab/brainseq_phase2/browser/BrainSeqPhaseII_eQTL_dlpfc_replication_bsp1.txt',
    col.names = c('snps', 'gene', 'statistic', 'pvalue', 'FDR', 'beta', 'Type')
)

# break up into pieces
message(paste(Sys.time(), 'breaking up by feature'))

message(paste(Sys.time(), 'gene'))
dlpfc_bsp1_genes <- bsp1[Type == 'gene']
message(paste(Sys.time(), 'exon'))
dlpfc_bsp1_exons <- bsp1[Type == 'exon']
message(paste(Sys.time(), 'jxn'))
dlpfc_bsp1_jxns <- bsp1[Type == 'jxn']
message(paste(Sys.time(), 'tx'))
dlpfc_bsp1_txs <- bsp1[Type == 'tx']

rm(bsp1)

# BrainSeq
message(paste(Sys.time(), 'loading BrainSeq Phase II eQTL results'))
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda", verbose=TRUE)

# keep only significant
message(paste(Sys.time(), 'subsetting to significant results'))
d_sig = data.table(as.data.frame(allEqtl[allEqtl$FDR < 0.01,]))
rm(allEqtl)

## Fix exon ids so they'll match with those from
## /dcl01/lieber/ajaffe/lab/brainseq_phase2/browser/BrainSeqPhaseII_eQTL_dlpfc_replication_bsp1.txt
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/browser/rda/exon_name_map.Rdata', verbose = TRUE)
setkey(exon_name_map, libd_bsp2)
d_sig$gene[d_sig$Type == 'Exon'] <- exon_name_map[.(d_sig$gene[d_sig$Type == 'Exon']), gencode]
rm(exon_name_map)

## Change jxn ids so they'll match as well
d_sig$gene[d_sig$Type == 'Jxn'] <- gsub('\\(\\+\\)|\\(\\-\\)', '(*)', d_sig$gene[d_sig$Type == 'Jxn'])

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

## Also write to smaller files for the browser
export_subset <- function(DT) {
    ## So the column matches the snpAnno column name
    colnames(DT)[1] <- 'snp'

    ## So the column matches the feature annotation column
    colnames(DT)[2] <- 'feature_id'

    ## Make Type lowercase to match file names from other tables
    DT$Type <- tolower(DT$Type)
    
    feature <- unique(DT$Type)
    feature <- feature[!is.na(feature)]

    f_new <- paste0('/dcl01/lieber/ajaffe/lab/brainseq_phase2/browser/BrainSeqPhaseII_eQTL_dlpfc_replication_bsp1_', feature, '.txt')
    message(paste(Sys.time(), 'writing', f_new))
    fwrite(DT, file = f_new, sep = '\t', row.names = FALSE)
}
export_subset(dlpfc_bsp1_genes)
export_subset(dlpfc_bsp1_exons)
export_subset(dlpfc_bsp1_jxn)
export_subset(dlpfc_bsp1_txs)

message(paste(Sys.time(), 'lines for each file'))
system('wc -l /dcl01/lieber/ajaffe/lab/brainseq_phase2/browser/BrainSeqPhaseII_eQTL_dlpfc_replication_bsp1_*')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
