library('data.table')
library('devtools')
library('SummarizedExperiment')

setDTthreads(1)

## Load subsets of data
files_sub <- dir('rdas', pattern = '_compare_', full.names = TRUE)
stopifnot(length(files_sub) == 3)
for(f in files_sub) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
}

## For testing
# cauc <- inter_cauc_genes
# brainseq <- i_sig_genes

merge_qtl <- function(cauc, brainseq) {
    ## All should be in the same order
    stopifnot(identical(cauc$snps, brainseq$snps))
    stopifnot(identical(cauc$gene, brainseq$gene))
    
    ## Use the latest Symbol info from CAUC files
    brainseq$Symbol <- cauc$Symbol
    
    ## Keep only a few CAUC columns
    to_add <- cauc[, c('statistic', 'pvalue', 'FDR', 'beta')]
    colnames(to_add) <- paste0('cauc_', colnames(to_add))
    
    ## Combine
    cbind(brainseq, to_add)
}

## Merge
message(paste(Sys.time(), 'merging interaction QTLs'))
interaction <- list(
    'gene' = merge_qtl(inter_cauc_genes, i_sig_genes)
)
message(paste(Sys.time(), 'saving merged QTLs'))
save(interaction, file = 'rdas/merged_CAUC_BrainSeq_QTLs_interaction.Rdata')

message(paste(Sys.time(), 'merging hippo QTLs'))
hippo <- list(
    'gene' = merge_qtl(hippo_cauc_genes, h_sig_genes)
)
message(paste(Sys.time(), 'saving merged QTLs'))
save(hippo, file = 'rdas/merged_CAUC_BrainSeq_QTLs_hippo.Rdata')


message(paste(Sys.time(), 'merging dlpfc QTLs'))
dlpfc <- list(
    'gene' = merge_qtl(dlpfc_cauc_genes, d_sig_genes)
)
message(paste(Sys.time(), 'saving merged QTLs'))
save(dlpfc, file = 'rdas/merged_CAUC_BrainSeq_QTLs_dlpfc.Rdata')


## Explore
comp_qtl <- function(type, dfs, perc = FALSE, cutde = 0.01) {
    df <- dfs[[type]]
    if(any(is.na(df$cauc_statistic)) & !perc) {
        message(paste(Sys.time(), 'removing some NAs from CAUC (TRUEs below) for type', type))
        print(table(is.na(df$cauc_statistic)))
    }
    
    res <- addmargins(table(sign(df$statistic) == sign(df$cauc_statistic),  df$cauc_pvalue < cutde, dnn = c('Equal sign', paste0('CAUC p<', cutde))))
    if(!perc) return(res)
    
    ## Calculate percent over all of brainseq
    ## the total marginal will not be 100% unless there were no NAs
    res / nrow(df) * 100
}
comp_qtl_short <- function(dfs, perc = FALSE, cutde = 0.01) {
    res <- lapply(names(dfs), comp_qtl, dfs = dfs, perc = perc, cutde = cutde)
    names(res) <- names(dfs)
    return(res)
}

comp_qtl_short(interaction)
# $gene
#           CAUC p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE  1780     0  1780
#      TRUE   7418 30894 38312
#      Sum    9198 30894 40092

comp_qtl_short(interaction, perc = TRUE)
# $gene
#           CAUC p<0.01
# Equal sign      FALSE       TRUE        Sum
#      FALSE   4.439788   0.000000   4.439788
#      TRUE   18.502444  77.057767  95.560212
#      Sum    22.942233  77.057767 100.000000

comp_qtl_short(interaction, cutde = 0.05)
# $gene
#           CAUC p<0.05
# Equal sign FALSE  TRUE   Sum
#      FALSE  1780     0  1780
#      TRUE   3351 34961 38312
#      Sum    5131 34961 40092

comp_qtl_short(interaction, perc = TRUE, cutde = 0.05)
# $gene
#           CAUC p<0.05
# Equal sign      FALSE       TRUE        Sum
#      FALSE   4.439788   0.000000   4.439788
#      TRUE    8.358276  87.201936  95.560212
#      Sum    12.798064  87.201936 100.000000

comp_qtl_short(hippo)
# $gene
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE   59441     123   59564
#      TRUE   221859  794310 1016169
#      Sum    281300  794433 1075733

comp_qtl_short(hippo, perc = TRUE)
# $gene
#           CAUC p<0.01
# Equal sign        FALSE         TRUE          Sum
#      FALSE   5.52562764   0.01143406   5.53706171
#      TRUE   20.62398383  73.83895446  94.46293829
#      Sum    26.14961147  73.85038853 100.00000000

comp_qtl_short(hippo, cutde = 0.05)
# $gene
#           CAUC p<0.05
# Equal sign   FALSE    TRUE     Sum
#      FALSE   59092     472   59564
#      TRUE   122322  893847 1016169
#      Sum    181414  894319 1075733

comp_qtl_short(hippo, perc = TRUE, cutde = 0.05)
# $gene
#           CAUC p<0.05
# Equal sign        FALSE         TRUE          Sum
#      FALSE   5.49318465   0.04387706   5.53706171
#      TRUE   11.37103724  83.09190106  94.46293829
#      Sum    16.86422188  83.13577812 100.00000000

comp_qtl_short(dlpfc)
# $gene
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE   92908     245   93153
#      TRUE   338138 1146672 1484810
#      Sum    431046 1146917 1577963

comp_qtl_short(dlpfc, perc = TRUE)
# $gene
#           CAUC p<0.01
# Equal sign        FALSE         TRUE          Sum
#      FALSE   5.88784401   0.01552635   5.90337036
#      TRUE   21.42876607  72.66786357  94.09662964
#      Sum    27.31661009  72.68338991 100.00000000

comp_qtl_short(dlpfc, cutde = 0.05)
# $gene
#           CAUC p<0.05
# Equal sign   FALSE    TRUE     Sum
#      FALSE   92405     748   93153
#      TRUE   185490 1299320 1484810
#      Sum    277895 1300068 1577963

comp_qtl_short(dlpfc, perc = TRUE, cutde = 0.05)
# $gene
#           CAUC p<0.05
# Equal sign        FALSE         TRUE          Sum
#      FALSE   5.85596747   0.04740289   5.90337036
#      TRUE   11.75502848  82.34160117  94.09662964
#      Sum    17.61099595  82.38900405 100.00000000


## Get original and new sample sizes
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata")


get_n <- function(index, group) {
    
    data.frame(
        n = sum(index),
        AA = sum(rse_gene$Race[index] == 'AA'),
        CAUC = sum(rse_gene$Race[index] == 'CAUC'),
        HISP = sum(rse_gene$Race[index] == 'HISP'),
        AS = sum(rse_gene$Race[index] == 'AS'),
        CAUC_percent = sum(rse_gene$Race[index] == 'CAUC') / sum(index) * 100,
        model = group,
        stringsAsFactors = FALSE
    )
    
}

sample_sizes <- rbind(
    get_n(colData(rse_gene)$Age > 13, 'interaction'),
    get_n(colData(rse_gene)$Age > 13 & colData(rse_gene)$Region == "DLPFC", 'DLPFC'),
    get_n(colData(rse_gene)$Age > 13 & colData(rse_gene)$Region == "HIPPO", 'HIPPO')
    ## Don't need the subsets
    # ,
    # get_n(colData(rse_gene)$Age > 13 & colData(rse_gene)$Race == 'CAUC', 'interaction - CAUC'),
    # get_n(colData(rse_gene)$Age > 13 & colData(rse_gene)$Region == "DLPFC" & colData(rse_gene)$Race == 'CAUC', 'DLPFC - CAUC'),
    # get_n(colData(rse_gene)$Age > 13 & colData(rse_gene)$Region == "HIPPO" & colData(rse_gene)$Race == 'CAUC', 'HIPPO - CAUC')
)
sample_sizes
#     n  AA CAUC HISP AS CAUC_percent       model
# 1 792 417  356   10  9     44.94949 interaction
# 2 397 204  174   10  9     43.82872       DLPFC
# 3 395 213  182    0  0     46.07595       HIPPO

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
#  date     2019-03-18
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  backports              1.1.3     2018-12-14 [2] CRAN (R 3.5.1)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.6    2019-02-10 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  callr                  3.1.1     2018-12-21 [2] CRAN (R 3.5.1)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  data.table           * 1.12.0    2019-01-13 [1] CRAN (R 3.5.1)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  desc                   1.2.0     2018-05-01 [2] CRAN (R 3.5.1)
#  devtools             * 2.0.1     2018-10-26 [1] CRAN (R 3.5.1)
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                  0.8.0.1   2019-02-15 [1] CRAN (R 3.5.1)
#  fs                     1.2.6     2018-08-23 [2] CRAN (R 3.5.1)
#  GenomeInfoDb         * 1.18.2    2019-02-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2                3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  later                  0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  memoise                1.1.0     2017-04-21 [2] CRAN (R 3.5.0)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgbuild               1.0.2     2018-10-16 [2] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  pkgload                1.0.2     2018-10-29 [2] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  prettyunits            1.0.2     2015-07-13 [1] CRAN (R 3.5.0)
#  processx               3.3.0     2019-03-10 [1] CRAN (R 3.5.1)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  ps                     1.3.0     2018-12-21 [2] CRAN (R 3.5.1)
#  purrr                  0.3.1     2019-03-03 [2] CRAN (R 3.5.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  remotes                2.0.2     2018-10-30 [1] CRAN (R 3.5.1)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  rprojroot              1.3-2     2018-01-03 [2] CRAN (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  servr                  0.13      2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo            1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  testthat               2.0.1     2018-10-13 [1] CRAN (R 3.5.1)
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  usethis              * 1.4.0     2018-08-14 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.5       2019-02-20 [1] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
