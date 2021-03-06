library('data.table')
library('sessioninfo')

## Load subsets of data
files_sub <- dir('rdas', pattern = '_compare_', full.names = TRUE)
stopifnot(length(files_sub) == 4)
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
message(paste(Sys.time(), 'saving merged QTLs'))
save(dlpfc, file = 'rdas/merged_BSP1_BrainSeq_QTLs_dlpfc.Rdata')


## Explore
comp_qtl <- function(type, dfs, perc = FALSE, cutde = 0.01) {
    df <- dfs[[type]]
    if(any(is.na(df$bsp1_statistic)) & !perc) {
        message(paste(Sys.time(), 'removing some NAs from BSP1 (TRUEs below) for type', type))
        print(table(is.na(df$bsp1_statistic)))
    }
    
    res <- addmargins(table(sign(df$statistic) == sign(df$bsp1_statistic),  df$bsp1_pvalue < cutde, dnn = c('Equal sign', paste0('BSP1 p<', cutde))))
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
# 2019-03-19 16:11:38 removing some NAs from BSP1 (TRUEs below) for type gene
#
#   FALSE    TRUE
# 1514240   63723
# 2019-03-19 16:11:39 removing some NAs from BSP1 (TRUEs below) for type exon
#
#   FALSE    TRUE
# 8620980  307329
# 2019-03-19 16:11:47 removing some NAs from BSP1 (TRUEs below) for type jxn
#
#   FALSE    TRUE
# 4639633  620493
# 2019-03-19 16:11:52 removing some NAs from BSP1 (TRUEs below) for type tx
#
#   FALSE    TRUE
# 2191553   79711
# $gene
#           BSP1 p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE   57902   21286   79188
#      TRUE   218108 1216944 1435052
#      Sum    276010 1238230 1514240
#
# $exon
#           BSP1 p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE  353861   96471  450332
#      TRUE  1237307 6933341 8170648
#      Sum   1591168 7029812 8620980
#
# $jxn
#           BSP1 p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE  177156   20938  198094
#      TRUE   811865 3629674 4441539
#      Sum    989021 3650612 4639633
#
# $tx
#           BSP1 p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE   90782   18695  109477
#      TRUE   340420 1741656 2082076
#      Sum    431202 1760351 2191553

comp_qtl_short(dlpfc, perc = TRUE)
# $gene
#           BSP1 p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE  3.669414  1.348954  5.018369
#      TRUE  13.822124 77.121200 90.943324
#      Sum   17.491538 78.470154 95.961692
#
# $tx
#           BSP1 p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE  3.996981  0.823110  4.820091
#      TRUE  14.988130 76.682235 91.670365
#      Sum   18.985111 77.505345 96.490456
#
# > comp_qtl_short(dlpfc, perc = TRUE)
# $gene
#           BSP1 p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE  3.669414  1.348954  5.018369
#      TRUE  13.822124 77.121200 90.943324
#      Sum   17.491538 78.470154 95.961692
#
# $exon
#           BSP1 p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE  3.963360  1.080507  5.043867
#      TRUE  13.858246 77.655702 91.513947
#      Sum   17.821605 78.736209 96.557814
#
# $jxn
#           BSP1 p<0.01
# Equal sign      FALSE       TRUE        Sum
#      FALSE  3.3679041  0.3980513  3.7659554
#      TRUE  15.4343261 69.0035562 84.4378823
#      Sum   18.8022302 69.4016075 88.2038377
#
# $tx
#           BSP1 p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE  3.996981  0.823110  4.820091
#      TRUE  14.988130 76.682235 91.670365
#      Sum   18.985111 77.505345 96.490456


comp_qtl_short(dlpfc, cutde = 0.05)
# $gene
#           BSP1 p<0.05
# Equal sign   FALSE    TRUE     Sum
#      FALSE   51983   27205   79188
#      TRUE   144519 1290533 1435052
#      Sum    196502 1317738 1514240
#
# $exon
#           BSP1 p<0.05
# Equal sign   FALSE    TRUE     Sum
#      FALSE  314920  135412  450332
#      TRUE   837198 7333450 8170648
#      Sum   1152118 7468862 8620980
#
# $jxn
#           BSP1 p<0.05
# Equal sign   FALSE    TRUE     Sum
#      FALSE  165798   32296  198094
#      TRUE   571094 3870445 4441539
#      Sum    736892 3902741 4639633
#
# $tx
#           BSP1 p<0.05
# Equal sign   FALSE    TRUE     Sum
#      FALSE   84053   25424  109477
#      TRUE   237753 1844323 2082076
#      Sum    321806 1869747 2191553

comp_qtl_short(dlpfc, perc = TRUE, cutde = 0.05)
# $gene
#           BSP1 p<0.05
# Equal sign     FALSE      TRUE       Sum
#      FALSE  3.294310  1.724058  5.018369
#      TRUE   9.158580 81.784744 90.943324
#      Sum   12.452890 83.508802 95.961692
#
# $exon
#           BSP1 p<0.05
# Equal sign     FALSE      TRUE       Sum
#      FALSE  3.527208  1.516659  5.043867
#      TRUE   9.376893 82.137054 91.513947
#      Sum   12.904101 83.653713 96.557814
#
# $jxn
#           BSP1 p<0.05
# Equal sign      FALSE       TRUE        Sum
#      FALSE  3.1519777  0.6139777  3.7659554
#      TRUE  10.8570403 73.5808420 84.4378823
#      Sum   14.0090180 74.1948197 88.2038377
#
# $tx
#           BSP1 p<0.05
# Equal sign     FALSE      TRUE       Sum
#      FALSE  3.700715  1.119377  4.820091
#      TRUE  10.467872 81.202493 91.670365
#      Sum   14.168586 82.321870 96.490456


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
#  date     2019-03-19
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date       lib source
#  assertthat    0.2.0   2017-04-11 [2] CRAN (R 3.5.0)
#  cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.1)
#  colorout    * 1.2-0   2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace    1.4-0   2019-01-13 [2] CRAN (R 3.5.1)
#  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)
#  data.table  * 1.12.0  2019-01-13 [1] CRAN (R 3.5.1)
#  digest        0.6.18  2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr         0.8.0.1 2019-02-15 [1] CRAN (R 3.5.1)
#  ggplot2       3.1.0   2018-10-25 [1] CRAN (R 3.5.1)
#  glue          1.3.1   2019-03-12 [1] CRAN (R 3.5.1)
#  gtable        0.2.0   2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools     0.3.6   2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets   1.3     2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv        1.4.5.1 2018-12-18 [2] CRAN (R 3.5.1)
#  jsonlite      1.6     2018-12-07 [2] CRAN (R 3.5.1)
#  later         0.8.0   2019-02-11 [2] CRAN (R 3.5.1)
#  lattice       0.20-38 2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval      0.2.1   2017-10-29 [2] CRAN (R 3.5.0)
#  magrittr      1.5     2014-11-22 [1] CRAN (R 3.5.0)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 3.5.1)
#  pillar        1.3.1   2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig     2.0.2   2018-08-16 [1] CRAN (R 3.5.1)
#  plyr          1.8.4   2016-06-08 [2] CRAN (R 3.5.0)
#  png           0.1-7   2013-12-03 [2] CRAN (R 3.5.0)
#  promises      1.0.1   2018-04-13 [2] CRAN (R 3.5.0)
#  purrr         0.3.1   2019-03-03 [2] CRAN (R 3.5.1)
#  R6            2.4.0   2019-02-14 [2] CRAN (R 3.5.1)
#  Rcpp          1.0.0   2018-11-07 [1] CRAN (R 3.5.1)
#  rlang         0.3.1   2019-01-08 [1] CRAN (R 3.5.1)
#  rmote       * 0.3.4   2018-05-02 [1] deltarho (R 3.5.0)
#  scales        1.0.0   2018-08-09 [2] CRAN (R 3.5.1)
#  servr         0.13    2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo * 1.1.1   2018-11-05 [1] CRAN (R 3.5.1)
#  tibble        2.0.1   2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect    0.2.5   2018-10-11 [2] CRAN (R 3.5.1)
#  withr         2.1.2   2018-03-15 [2] CRAN (R 3.5.0)
#  xfun          0.5     2019-02-20 [1] CRAN (R 3.5.1)
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library

