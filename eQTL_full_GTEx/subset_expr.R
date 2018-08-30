library('SummarizedExperiment')
library('GenomicRanges')
library('devtools')

dir.create('rdas', showWarnings = FALSE)

## Load GTEx data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_both/rse_gtex_gene.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_both/rse_gtex_exon.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_both/rse_gtex_jxn.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_both/rse_gtex_tx.Rdata", verbose = TRUE)

keepInfo <- vector('list', 5)
names(keepInfo) <- c('DLPFC', 'HIPPO', 'interaction', 'DLPFC_raggr', 'HIPPO_raggr')

## Load sig results
## DLPFC full first
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda', verbose = TRUE)
keepInfo[['DLPFC']] <- allEqtl[allEqtl$FDR < 0.05, c('snps', 'gene', 'Type')]

## HIPPO full second
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_hippo_4features.rda', verbose = TRUE)
keepInfo[['HIPPO']] <- allEqtl[allEqtl$FDR < 0.05, c('snps', 'gene', 'Type')]

## Interaction next
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/matrixEqtl_output_interaction_4features.rda', verbose = TRUE)
keepInfo[['interaction']] <- DataFrame(rbind(
        cbind(meGene$cis$eqtls[meGene$cis$eqtls$FDR < 0.05, c('snps', 'gene')], 'Type' = 'Gene'),
        cbind(meExon$cis$eqtls[meExon$cis$eqtls$FDR < 0.05, c('snps', 'gene')], 'Type' = 'Exon'),
        cbind(meJxn$cis$eqtls[meJxn$cis$eqtls$FDR < 0.05, c('snps', 'gene')], 'Type' = 'Jxn'),
        cbind(meTx$cis$eqtls[meTx$cis$eqtls$FDR < 0.05, c('snps', 'gene')], 'Type' = 'Tx')
))
keepInfo[['interaction']]$snps <- as.character(keepInfo[['interaction']]$snps)
keepInfo[['interaction']]$gene <- as.character(keepInfo[['interaction']]$gene)
keepInfo[['interaction']]$Type <- as.character(keepInfo[['interaction']]$Type)
rm(meGene, meExon, meJxn, meTx)

## PGC2 rAggr DLPFC subset
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/eqtl_tables/mergedEqtl_output_dlpfc_raggr_4features.rda', verbose = TRUE)
keepInfo[['DLPFC_raggr']] <- allEqtl[allEqtl$FDR < 0.05, c('snps', 'gene', 'Type')]

## finally PGC2 rAggr HIPPO subset
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/eqtl_tables/mergedEqtl_output_hippo_raggr_4features.rda', verbose = TRUE)
keepInfo[['HIPPO_raggr']] <- allEqtl[allEqtl$FDR < 0.05, c('snps', 'gene', 'Type')]
rm(allEqtl)

## Save info for later if needed
save(keepInfo, file = 'rdas/keepInfo.Rdata')

## Find unique
find_uniq <- function(df, type, var) {
    unique(df[df$Type == type, var])
}
get_uniq <- function(var, df_list) {
    print(var)
    types <- c('Gene', 'Exon', 'Jxn', 'Tx')
    res <- lapply(types, function(type) {
        unique(unlist(lapply(df_list, find_uniq, type = type, var = var)))
    })
    names(res) <- types
    print(sapply(res, length))
    return(res)
}

uniqueInfo <- lapply(c('snps', 'gene'), get_uniq, df_list = keepInfo)
# [1] "snps"
#    Gene    Exon     Jxn      Tx
# 1397500 2132252 2298065 1600406
# [1] "gene"
#   Gene   Exon    Jxn     Tx
#  21131 276409 172858  64577
names(uniqueInfo) <- c('snps', 'gene')
save(uniqueInfo, file = 'rdas/uniqueInfo.Rdata')

sapply(uniqueInfo, function(x) { sum(sapply(x, length)) })
#    snps    gene
# 7428223  534975

length(unique(unlist(uniqueInfo$snps)))
# [1] 3197335

## Exon and jxn are special cases with the names and all that
## Here I finally fixed the names so they'll match the previous objects
rse_gtex_gene <- rse_gtex_gene[match(uniqueInfo$gene$Gene, rownames(rse_gtex_gene)), ]
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata', verbose = TRUE)
ov <- findOverlaps(rowRanges(rse_exon), rowRanges(rse_gtex_exon), type = 'equal', ignore.strand = FALSE)
m <- subjectHits(ov)
rse_gtex_exon <- rse_gtex_exon[m[match(uniqueInfo$gene$Exon, rownames(rse_exon))], ]
rownames(rse_gtex_exon) <- rownames(rse_exon)[match(uniqueInfo$gene$Exon, rownames(rse_exon))]
rse_gtex_jxn <- rse_gtex_jxn[match(uniqueInfo$gene$Jxn, with(rowRanges(rse_gtex_jxn), paste0(seqnames, ':', start, '-', end, '(', strand, ')'))), ]
rownames(rse_gtex_jxn) <-  with(rowRanges(rse_gtex_jxn), paste0(seqnames, ':', start, '-', end, '(', strand, ')'))
rse_gtex_tx <- rse_gtex_tx[match(uniqueInfo$gene$Tx, rownames(rse_gtex_tx)), ]

save(rse_gtex_gene, file = 'rdas/rse_gtex_gene_subset.Rdata')
save(rse_gtex_exon, file = 'rdas/rse_gtex_exon_subset.Rdata')
save(rse_gtex_jxn, file = 'rdas/rse_gtex_jxn_subset.Rdata')
save(rse_gtex_tx, file = 'rdas/rse_gtex_tx_subset.Rdata')


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.5.0 Patched (2018-04-30 r74679)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       US/Eastern
#  date     2018-08-30
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  assertthat             0.2.0     2017-04-11 CRAN (R 3.5.0)
#  base                 * 3.5.0     2018-05-02 local
#  bindr                  0.1.1     2018-03-13 CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 CRAN (R 3.5.0)
#  Biobase              * 2.40.0    2018-05-02 Bioconductor
#  BiocGenerics         * 0.26.0    2018-05-03 Bioconductor
#  BiocParallel         * 1.14.2    2018-07-08 Bioconductor
#  bitops                 1.0-6     2013-08-17 CRAN (R 3.5.0)
#  colorout             * 1.2-0     2018-05-02 Github (jalvesaq/colorout@c42088d)
#  colorspace             1.3-2     2016-12-14 CRAN (R 3.5.0)
#  compiler               3.5.0     2018-05-02 local
#  crayon                 1.3.4     2017-09-16 CRAN (R 3.5.0)
#  datasets             * 3.5.0     2018-05-02 local
#  DelayedArray         * 0.6.2     2018-07-23 Bioconductor
#  devtools             * 1.13.6    2018-06-27 CRAN (R 3.5.0)
#  digest                 0.6.15    2018-01-28 CRAN (R 3.5.0)
#  dplyr                  0.7.6     2018-06-29 CRAN (R 3.5.0)
#  GenomeInfoDb         * 1.16.0    2018-05-03 Bioconductor
#  GenomeInfoDbData       1.1.0     2018-04-17 Bioconductor
#  GenomicRanges        * 1.32.6    2018-07-20 Bioconductor
#  ggplot2                3.0.0     2018-07-03 CRAN (R 3.5.0)
#  glue                   1.3.0     2018-07-17 CRAN (R 3.5.0)
#  graphics             * 3.5.0     2018-05-02 local
#  grDevices            * 3.5.0     2018-05-02 local
#  grid                   3.5.0     2018-05-02 local
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 CRAN (R 3.5.0)
#  htmlwidgets            1.2       2018-04-19 CRAN (R 3.5.0)
#  httpuv                 1.4.5     2018-07-19 CRAN (R 3.5.0)
#  IRanges              * 2.14.10   2018-05-17 Bioconductor
#  later                  0.7.3     2018-06-08 CRAN (R 3.5.0)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.5.0)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 CRAN (R 3.5.0)
#  Matrix                 1.2-14    2018-04-13 CRAN (R 3.5.0)
#  matrixStats          * 0.54.0    2018-07-23 CRAN (R 3.5.0)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.5.0)
#  methods              * 3.5.0     2018-05-02 local
#  munsell                0.5.0     2018-06-12 CRAN (R 3.5.0)
#  parallel             * 3.5.0     2018-05-02 local
#  pillar                 1.3.0     2018-07-14 CRAN (R 3.5.0)
#  pkgconfig              2.0.1     2017-03-21 CRAN (R 3.5.0)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 CRAN (R 3.5.0)
#  purrr                  0.2.5     2018-05-29 CRAN (R 3.5.0)
#  R6                     2.2.2     2017-06-17 CRAN (R 3.5.0)
#  Rcpp                   0.12.18   2018-07-23 CRAN (R 3.5.0)
#  RCurl                  1.95-4.11 2018-07-15 CRAN (R 3.5.0)
#  rlang                  0.2.1     2018-05-30 cran (@0.2.1)
#  rmote                * 0.3.4     2018-05-02 deltarho (R 3.5.0)
#  S4Vectors            * 0.18.3    2018-06-13 Bioconductor
#  scales                 0.5.0     2017-08-24 CRAN (R 3.5.0)
#  servr                  0.10      2018-05-30 CRAN (R 3.5.0)
#  stats                * 3.5.0     2018-05-02 local
#  stats4               * 3.5.0     2018-05-02 local
#  SummarizedExperiment * 1.10.1    2018-05-17 Bioconductor
#  tibble                 1.4.2     2018-01-22 CRAN (R 3.5.0)
#  tidyselect             0.2.4     2018-02-26 CRAN (R 3.5.0)
#  tools                  3.5.0     2018-05-02 local
#  utils                * 3.5.0     2018-05-02 local
#  withr                  2.1.2     2018-03-15 CRAN (R 3.5.0)
#  xfun                   0.3       2018-07-06 CRAN (R 3.5.0)
#  XVector                0.20.0    2018-05-03 Bioconductor
#  zlibbioc               1.26.0    2018-05-02 Bioconductor
