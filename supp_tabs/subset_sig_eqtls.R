library('data.table')
library('sessioninfo')
library('SummarizedExperiment')

## Load the significant results with GTEx replication
files_sub <- dir('../eQTL_full_GTEx/rdas', pattern = 'merged_GTEx_BrainSeq_QTLs', full.names = TRUE)
stopifnot(length(files_sub) == 3)
for(f in files_sub) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
}

## Now load the gene-level CAUC-only replication
files_sub <- dir('../eQTL_full/check_cauc/rdas', pattern = '_compare_', full.names = TRUE)
stopifnot(length(files_sub) == 3)

for(f in files_sub) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
}

## Merge with the CAUC-only replication
merge_qtl_cauc <- function(cauc, brainseq) {
    ## All should be in the same order
    stopifnot(identical(cauc$snps, brainseq$snps))
    stopifnot(identical(cauc$gene, brainseq$gene))
    
    # ## Use the latest Symbol info from CAUC files
    # brainseq$Symbol <- cauc$Symbol
    
    ## Keep only a few CAUC columns
    to_add <- cauc[, c('statistic', 'pvalue', 'FDR', 'beta')]
    colnames(to_add) <- paste0('cauc_', colnames(to_add))
    
    ## Combine
    cbind(brainseq, to_add)
}

interaction$gene <- merge_qtl_cauc(inter_cauc_genes, interaction$gene)
hippo$gene <- merge_qtl_cauc(hippo_cauc_genes, hippo$gene)
dlpfc$gene <- merge_qtl_cauc(dlpfc_cauc_genes, dlpfc$gene)

## Remove used objects
rm(dlpfc_cauc_genes, d_sig_genes, hippo_cauc_genes, h_sig_genes, inter_cauc_genes, i_sig_genes)


## Now load the BrainSeq Phase I DLPFC polyA+ replication
files_sub <- dir('../bsp1/eqtl/full/rdas', pattern = '_compare_', full.names = TRUE)
stopifnot(length(files_sub) == 4)
for(f in files_sub) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
}

## Free some space
rm(d_sig_exons, d_sig_genes, d_sig_jxns, d_sig_txs)

## Merge with BSP1
merge_qtl_bsp1 <- function(bsp1, brainseq) {    
    ## All should be in the same order
    stopifnot(identical(bsp1$snps, brainseq$snps))
    stopifnot(identical(bsp1$gene, brainseq$gene))
    
    # ## Use the latest Symbol info from BSP1 files
    # brainseq$Symbol <- bsp1$Symbol
    
    ## Keep only a few BSP1 columns
    to_add <- bsp1[, c('statistic', 'pvalue', 'FDR', 'beta')]
    colnames(to_add) <- paste0('bsp1_', colnames(to_add))
    
    ## Combine
    cbind(brainseq, to_add)
}

## Merge
message(paste(Sys.time(), 'merging dlpfc QTLs'))
x <- dlpfc

## Due to BSP1 being unstranded
original <- dlpfc$jxn$gene 
dlpfc$jxn$gene <- gsub('\\(\\+\\)|\\(\\-\\)', '(*)', dlpfc$jxn$gene)

## Fix exon ids so they'll match with those from
## /dcl01/lieber/ajaffe/lab/brainseq_phase2/browser/BrainSeqPhaseII_eQTL_dlpfc_replication_bsp1.txt
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/browser/rda/exon_name_map.Rdata', verbose = TRUE)
setkey(exon_name_map, libd_bsp2)

## Do the matching manually...
dlpfc_exon <- dlpfc$exon
dlpfc_exon$original_gene <- dlpfc_exon$gene
dlpfc_exon$gene <- exon_name_map[.(dlpfc_exon$gene), gencode]
setkey(dlpfc_exon, snps, gene)

## Merge
dlpfc_exon <- merge_qtl_bsp1(dlpfc_bsp1_exons, dlpfc_exon)

## Restore original ids
dlpfc_exon$gene <- dlpfc_exon$original_gene
# https://stackoverflow.com/questions/9202413/how-do-you-delete-a-column-by-name-in-data-table
dlpfc_exon <- dlpfc_exon[, original_gene:=NULL]
setkey(dlpfc_exon, snps, gene)

dlpfc <- list(
    'gene' = merge_qtl_bsp1(dlpfc_bsp1_genes, dlpfc$gene),
    'exon' = dlpfc_exon,
    'jxn' = merge_qtl_bsp1(dlpfc_bsp1_jxns, dlpfc$jxn),
    'tx' = merge_qtl_bsp1(dlpfc_bsp1_txs, dlpfc$tx)
)

## Restore the original jxn ids
dlpfc$jxn$gene <- original

## Clean up
rm(dlpfc_bsp1_genes, dlpfc_bsp1_exons, dlpfc_bsp1_jxns, dlpfc_bsp1_txs, dlpfc_exon)


## Check against
## https://github.com/LieberInstitute/brainseq_phase2/blob/master/bsp1/eqtl/full/bsp1_eqtl_replication.R#L156-L161
# addmargins(with(dlpfc$gene, table(sign(statistic) == sign(bsp1_statistic),  bsp1_pvalue < 0.05, dnn = c('Equal sign', paste0('BSP1 p<', 0.05)))))


## Change the first column names, and make CSV friendly
fix_tables <- function(DT) {
    ## Based on https://github.com/LieberInstitute/brainseq_phase2/blob/master/browser/extract_data.R#L384-L413
    
    colnames(DT)[1] <- 'snp'
    colnames(DT)[2] <- 'feature_id'
    
    ## Make Type lowercase to match file names from other tables
    DT$Type <- tolower(DT$Type)
    
    if('gencodeTx' %in% colnames(DT)) {
        ## This is too slow...
        # message(paste(Sys.time(), 'changing commas for semicolons for the gencodeTx column'))
        # DT$gencodeTx <- gsub(',', ';', DT$gencodeTx)
        ## Simply drop it
        # DT <- DT[, gencodeTx:=NULL]
        
        ## Ok, this is not too bad!
        message(paste(Sys.time(), 'collapsing gencodeTx column'))
        DT$gencodeTx <- sapply(DT$gencodeTx, paste, collapse = ';')
    }
    
    
    
    return(DT)
}

dlpfc <- lapply(dlpfc, fix_tables)
hippo <- lapply(hippo, fix_tables)
interaction <- lapply(interaction, fix_tables)

## Add Gencode Exon IDs
interaction$exon$exonGencodeID <- exon_name_map[.(interaction$exon$feature_id), gencode]
dlpfc$exon$exonGencodeID <- exon_name_map[.(dlpfc$exon$feature_id), gencode]
hippo$exon$exonGencodeID <- exon_name_map[.(hippo$exon$feature_id), gencode]


## Add gene type info too
load('../expr_cutoff/unfiltered/rse_gene_unfiltered.Rdata', verbose = TRUE)
gene_type_map <- data.table(gene_type = rowRanges(rse_gene)$gene_type, ensemblID = rowRanges(rse_gene)$ensemblID)
setkey(gene_type_map, ensemblID)

add_gene_type <- function(DT) {
    DT$gene_type <- gene_type_map[.(DT$EnsemblGeneID), gene_type]
    return(DT)
}
## Could have done this in fix_tables() but forgot so when I wrote this code...

dlpfc <- lapply(dlpfc, add_gene_type)
hippo <- lapply(hippo, add_gene_type)
interaction <- lapply(interaction, add_gene_type)


export_tables <- function(DT, feature, set) {
    f_new <- paste0('BrainSeqPhaseII_eQTL_FDR1perc_', set, '_', feature, '.txt')
    message(paste(Sys.time(), 'writing', f_new))
    fwrite(DT, file = f_new, sep = '\t', row.names = FALSE)
    return(file.exists(f_new))
}

mapply(export_tables, interaction, names(interaction), set = 'Interaction')
mapply(export_tables, hippo, names(hippo), set = 'HIPPO')
mapply(export_tables, dlpfc, names(dlpfc), set = 'DLPFC')

system('tar -cvzf SupplementaryTableXX_eQTL.tar.gz BrainSeqPhaseII_eQTL_FDR1perc_*.txt')
system('ls -lh SupplementaryTableXX_eQTL.tar.gz')
# -rw-r--r-- 1 lcollado lieber_jaffe 3.0G Mar 29 14:03 SupplementaryTableXX_eQTL.tar.gz

## Save R objects too
dir.create('rda', showWarnings = FALSE)
message(paste(Sys.time(), 'saving merged interaction QTLs'))
save(interaction, file = 'rda/merged_CAUC_GTEX_BrainSeq_QTLs_interaction.Rdata')

message(paste(Sys.time(), 'saving merged hippo QTLs'))
save(hippo, file = 'rda/merged_CAUC_GTEx_BrainSeq_QTLs_hippo.Rdata')

message(paste(Sys.time(), 'saving merged dlpfc QTLs'))
save(dlpfc, file = 'rda/merged_CAUC_GTEx_BSP1_BrainSeq_QTLs_dlpfc.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.5.1 Patched (2018-10-29 r75535)
#  os       Red Hat Enterprise Linux Server release 6.9 (Santiago)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-03-29
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.5.1)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.6    2019-02-10 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  data.table           * 1.12.0    2019-01-13 [1] CRAN (R 3.5.1)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                  0.8.0.1   2019-02-15 [1] CRAN (R 3.5.1)
#  GenomeInfoDb         * 1.18.2    2019-02-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2                3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.5.1)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.5.0     2019-03-15 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  later                  0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.5.1)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                  0.3.2     2019-03-15 [2] CRAN (R 3.5.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  servr                  0.13      2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.5       2019-02-20 [1] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
