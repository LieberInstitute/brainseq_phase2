library('data.table')
library('devtools')

## Load subsets of data
files_sub <- dir('rdas', pattern = '_compare_', full.names = TRUE)
stopifnot(length(files_sub) == 12)
for(f in files_sub) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
}

## For testing
# bsp1 <- inter_bsp1_genes
# brainseq <- i_sig_genes

merge_qtl <- function(bsp1, brainseq) {
    ## All should be in the same order
    stopifnot(identical(bsp1$snps, brainseq$snps))
    stopifnot(identical(bsp1$gene, brainseq$gene))
    
    ## Use the latest Symbol info from BSP1 files
    brainseq$Symbol <- bsp1$Symbol
    
    ## Keep only a few BSP1 columns
    to_add <- bsp1[, c('statistic', 'pvalue', 'FDR', 'beta')]
    colnames(to_add) <- paste0('bsp1_', colnames(to_add))
    
    ## Combine
    cbind(brainseq, to_add)
}

## Merge
message(paste(Sys.time(), 'merging dlpfc QTLs'))
dlpfc <- list(
    'gene' = merge_qtl(dlpfc_bsp1_genes, d_sig_genes),
    'exon' = merge_qtl(dlpfc_bsp1_exons, d_sig_exons),
    'jxn' = merge_qtl(dlpfc_bsp1_jxns, d_sig_jxns),
    'tx' = merge_qtl(dlpfc_bsp1_txs, d_sig_txs)
)
#interaction_all <- do.call(rbind, dlpfc)
message(paste(Sys.time(), 'saving merged QTLs'))
save(dlpfc, file = 'rdas/merged_BSP1_BrainSeq_QTLs_dlpfc.Rdata')


## Explore
comp_qtl <- function(type, dfs, perc = FALSE, cutde = 0.01) {
    df <- dfs[[type]]
    if(any(is.na(df$bsp1_statistic)) & !perc) {
        message(paste(Sys.time(), 'removing some NAs from BSP1 (TRUEs below) for type', type))
        print(table(is.na(df$bsp1_statistic)))
    }
    
    res <- addmargins(table('Equal sign' = sign(df$statistic) == sign(df$bsp1_statistic), 'BSP1 p<0.01' = df$bsp1_pvalue < cutde))
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

comp_qtl_short(dlpfc)
# 2018-09-18 21:56:41 removing some NAs from BSP1 (TRUEs below) for type gene
#
#   FALSE    TRUE
# 1371570  206393
# 2018-09-18 21:56:44 removing some NAs from BSP1 (TRUEs below) for type exon
#
#   FALSE    TRUE
# 7805138 1123171
# 2018-09-18 21:57:00 removing some NAs from BSP1 (TRUEs below) for type jxn
#
#   FALSE    TRUE
# 4630635  629491
# 2018-09-18 21:57:09 removing some NAs from BSP1 (TRUEs below) for type tx
#
#   FALSE    TRUE
# 1994482  276782
# $gene
#           BSP1 p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE  244982   22214  267196
#      TRUE   750876  353498 1104374
#      Sum    995858  375712 1371570
#
# $exon
#           BSP1 p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE 1701101  140661 1841762
#      TRUE  4187885 1775491 5963376
#      Sum   5888986 1916152 7805138
#
# $jxn
#           BSP1 p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE 1277449   47132 1324581
#      TRUE  2248301 1057753 3306054
#      Sum   3525750 1104885 4630635
#
# $tx
#           BSP1 p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE  348925   28295  377220
#      TRUE  1044562  572700 1617262
#      Sum   1393487  600995 1994482
#

comp_qtl_short(dlpfc, perc = TRUE)
# $gene
#           BSP1 p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 15.525206  1.407764 16.932970
#      TRUE  47.585146 22.402173 69.987319
#      Sum   63.110352 23.809937 86.920289
#
# $exon
#           BSP1 p<0.01
# Equal sign    FALSE     TRUE      Sum
#      FALSE 19.05289  1.57545 20.62834
#      TRUE  46.90569 19.88608 66.79177
#      Sum   65.95858 21.46153 87.42012
#
# $jxn
#           BSP1 p<0.01
# Equal sign      FALSE       TRUE        Sum
#      FALSE 24.2855209  0.8960242 25.1815451
#      TRUE  42.7423412 20.1088909 62.8512321
#      Sum   67.0278621 21.0049151 88.0327772
#
# $tx
#           BSP1 p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 15.362591  1.245782 16.608373
#      TRUE  45.990338 25.215034 71.205373
#      Sum   61.352929 26.460817 87.813746

comp_qtl_short(dlpfc, cutde = 0.05)
# $gene
#           BSP1 p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE  231205   35991  267196
#      TRUE   592282  512092 1104374
#      Sum    823487  548083 1371570
#
# $exon
#           BSP1 p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE 1584302  257460 1841762
#      TRUE  3374285 2589091 5963376
#      Sum   4958587 2846551 7805138
#
# $jxn
#           BSP1 p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE 1246207   78374 1324581
#      TRUE  1799967 1506087 3306054
#      Sum   3046174 1584461 4630635
#
# $tx
#           BSP1 p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE  329602   47618  377220
#      TRUE   827367  789895 1617262
#      Sum   1156969  837513 1994482

comp_qtl_short(dlpfc, perc = TRUE, cutde = 0.05)
# $gene
#           BSP1 p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 14.652118  2.280852 16.932970
#      TRUE  37.534594 32.452725 69.987319
#      Sum   52.186712 34.733577 86.920289
#
# $exon
#           BSP1 p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 17.744704  2.883637 20.628341
#      TRUE  37.793103 28.998672 66.791774
#      Sum   55.537807 31.882308 87.420115
#
# $jxn
#           BSP1 p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 23.691581  1.489964 25.181545
#      TRUE  34.219085 28.632147 62.851232
#      Sum   57.910666 30.122111 88.032777
#
# $tx
#           BSP1 p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 14.511831  2.096542 16.608373
#      TRUE  36.427602 34.777771 71.205373
#      Sum   50.939433 36.874313 87.813746


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
#  date     2018-09-18
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package     * version date       source
#  assertthat    0.2.0   2017-04-11 CRAN (R 3.5.0)
#  base        * 3.5.0   2018-05-02 local
#  bindr         0.1.1   2018-03-13 CRAN (R 3.5.0)
#  bindrcpp      0.2.2   2018-03-29 CRAN (R 3.5.0)
#  colorout    * 1.2-0   2018-05-02 Github (jalvesaq/colorout@c42088d)
#  colorspace    1.3-2   2016-12-14 CRAN (R 3.5.0)
#  compiler      3.5.0   2018-05-02 local
#  crayon        1.3.4   2017-09-16 CRAN (R 3.5.0)
#  data.table  * 1.11.4  2018-05-27 cran (@1.11.4)
#  datasets    * 3.5.0   2018-05-02 local
#  devtools    * 1.13.6  2018-06-27 CRAN (R 3.5.0)
#  digest        0.6.15  2018-01-28 CRAN (R 3.5.0)
#  dplyr         0.7.6   2018-06-29 CRAN (R 3.5.0)
#  ggplot2       3.0.0   2018-07-03 CRAN (R 3.5.0)
#  glue          1.3.0   2018-07-17 CRAN (R 3.5.0)
#  graphics    * 3.5.0   2018-05-02 local
#  grDevices   * 3.5.0   2018-05-02 local
#  grid          3.5.0   2018-05-02 local
#  gtable        0.2.0   2016-02-26 CRAN (R 3.5.0)
#  htmltools     0.3.6   2017-04-28 CRAN (R 3.5.0)
#  htmlwidgets   1.2     2018-04-19 CRAN (R 3.5.0)
#  httpuv        1.4.5   2018-07-19 CRAN (R 3.5.0)
#  later         0.7.4   2018-08-31 CRAN (R 3.5.0)
#  lattice       0.20-35 2017-03-25 CRAN (R 3.5.0)
#  lazyeval      0.2.1   2017-10-29 CRAN (R 3.5.0)
#  magrittr      1.5     2014-11-22 CRAN (R 3.5.0)
#  memoise       1.1.0   2017-04-21 CRAN (R 3.5.0)
#  methods     * 3.5.0   2018-05-02 local
#  munsell       0.5.0   2018-06-12 CRAN (R 3.5.0)
#  pillar        1.3.0   2018-07-14 CRAN (R 3.5.0)
#  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.5.0)
#  plyr          1.8.4   2016-06-08 CRAN (R 3.5.0)
#  png           0.1-7   2013-12-03 CRAN (R 3.5.0)
#  promises      1.0.1   2018-04-13 CRAN (R 3.5.0)
#  purrr         0.2.5   2018-05-29 CRAN (R 3.5.0)
#  R6            2.2.2   2017-06-17 CRAN (R 3.5.0)
#  Rcpp          0.12.18 2018-07-23 CRAN (R 3.5.0)
#  rlang         0.2.1   2018-05-30 cran (@0.2.1)
#  rmote       * 0.3.4   2018-05-02 deltarho (R 3.5.0)
#  scales        1.0.0   2018-08-09 CRAN (R 3.5.0)
#  servr         0.10    2018-05-30 CRAN (R 3.5.0)
#  stats       * 3.5.0   2018-05-02 local
#  tibble        1.4.2   2018-01-22 CRAN (R 3.5.0)
#  tidyselect    0.2.4   2018-02-26 CRAN (R 3.5.0)
#  tools         3.5.0   2018-05-02 local
#  utils       * 3.5.0   2018-05-02 local
#  withr         2.1.2   2018-03-15 CRAN (R 3.5.0)
#  xfun          0.3     2018-07-06 CRAN (R 3.5.0)

