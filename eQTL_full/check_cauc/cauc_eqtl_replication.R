library('data.table')
library('devtools')

setDTthreads(1)

## Load subsets of data
files_sub <- dir('rdas', pattern = '_compare_', full.names = TRUE)
stopifnot(length(files_sub) == 12)
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
    'gene' = merge_qtl(inter_cauc_genes, i_sig_genes),
    'exon' = merge_qtl(inter_cauc_exons, i_sig_exons),
    'jxn' = merge_qtl(inter_cauc_jxns, i_sig_jxns),
    'tx' = merge_qtl(inter_cauc_txs, i_sig_txs)
)
#interaction_all <- do.call(rbind, interaction)
message(paste(Sys.time(), 'saving merged QTLs'))
save(interaction, file = 'rdas/merged_CAUC_BrainSeq_QTLs_interaction.Rdata')

message(paste(Sys.time(), 'merging hippo QTLs'))
hippo <- list(
    'gene' = merge_qtl(hippo_cauc_genes, h_sig_genes),
    'exon' = merge_qtl(hippo_cauc_exons, h_sig_exons),
    'jxn' = merge_qtl(hippo_cauc_jxns, h_sig_jxns),
    'tx' = merge_qtl(hippo_cauc_txs, h_sig_txs)
)
#hippo_all <- do.call(rbind, hippo)
message(paste(Sys.time(), 'saving merged QTLs'))
save(hippo, file = 'rdas/merged_CAUC_BrainSeq_QTLs_hippo.Rdata')


message(paste(Sys.time(), 'merging dlpfc QTLs'))
dlpfc <- list(
    'gene' = merge_qtl(dlpfc_cauc_genes, d_sig_genes),
    'exon' = merge_qtl(dlpfc_cauc_exons, d_sig_exons),
    'jxn' = merge_qtl(dlpfc_cauc_jxns, d_sig_jxns),
    'tx' = merge_qtl(dlpfc_cauc_txs, d_sig_txs)
)
#interaction_all <- do.call(rbind, dlpfc)
message(paste(Sys.time(), 'saving merged QTLs'))
save(dlpfc, file = 'rdas/merged_CAUC_BrainSeq_QTLs_dlpfc.Rdata')


## Explore
comp_qtl <- function(type, dfs, perc = FALSE, cutde = 0.01) {
    df <- dfs[[type]]
    if(any(is.na(df$cauc_statistic)) & !perc) {
        message(paste(Sys.time(), 'removing some NAs from CAUC (TRUEs below) for type', type))
        print(table(is.na(df$cauc_statistic)))
    }
    
    res <- addmargins(table('Equal sign' = sign(df$statistic) == sign(df$cauc_statistic), 'CAUC p<0.01' = df$cauc_pvalue < cutde))
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
# 2018-09-18 17:13:37 removing some NAs from CAUC (TRUEs below) for type gene
#
# FALSE  TRUE
# 35958  4134
# 2018-09-18 17:13:37 removing some NAs from CAUC (TRUEs below) for type exon
#
# FALSE  TRUE
# 79435 10488
# 2018-09-18 17:13:37 removing some NAs from CAUC (TRUEs below) for type jxn
#
# FALSE  TRUE
# 66829  8774
# 2018-09-18 17:13:37 removing some NAs from CAUC (TRUEs below) for type tx
#
# FALSE  TRUE
# 19658  2237
# $gene
#           CAUC p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE  8613   334  8947
#      TRUE  20852  6159 27011
#      Sum   29465  6493 35958
#
# $exon
#           CAUC p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE 17026   759 17785
#      TRUE  39580 22070 61650
#      Sum   56606 22829 79435
#
# $jxn
#           CAUC p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE 24361   639 25000
#      TRUE  32179  9650 41829
#      Sum   56540 10289 66829
#
# $tx
#           CAUC p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE  2774   234  3008
#      TRUE  11077  5573 16650
#      Sum   13851  5807 19658

comp_qtl_short(interaction, perc = TRUE)
# $gene
#           CAUC p<0.01
# Equal sign      FALSE       TRUE        Sum
#      FALSE 21.4830889  0.8330839 22.3161728
#      TRUE  52.0103761 15.3621670 67.3725432
#      Sum   73.4934650 16.1952509 89.6887160
#
# $exon
#           CAUC p<0.01
# Equal sign      FALSE       TRUE        Sum
#      FALSE 18.9339768  0.8440555 19.7780323
#      TRUE  44.0154354 24.5432203 68.5586557
#      Sum   62.9494123 25.3872758 88.3366881
#
# $jxn
#           CAUC p<0.01
# Equal sign      FALSE       TRUE        Sum
#      FALSE 32.2222663  0.8452046 33.0674709
#      TRUE  42.5631258 12.7640438 55.3271696
#      Sum   74.7853921 13.6092483 88.3946404
#
# $tx
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 12.669559  1.068737 13.738296
#      TRUE  50.591459 25.453300 76.044759
#      Sum   63.261018 26.522037 89.783055

comp_qtl_short(interaction, cutde = 0.05)
# $gene
#           CAUC p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE  8425   522  8947
#      TRUE  17284  9727 27011
#      Sum   25709 10249 35958
#
# $exon
#           CAUC p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE 16158  1627 17785
#      TRUE  33357 28293 61650
#      Sum   49515 29920 79435
#
# $jxn
#           CAUC p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE 23903  1097 25000
#      TRUE  26352 15477 41829
#      Sum   50255 16574 66829
#
# $tx
#           CAUC p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE  2539   469  3008
#      TRUE   7855  8795 16650
#      Sum   10394  9264 19658

comp_qtl_short(interaction, perc = TRUE, cutde = 0.05)
# $gene
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 21.014167  1.302005 22.316173
#      TRUE  43.110845 24.261698 67.372543
#      Sum   64.125012 25.563703 89.688716
#
# $exon
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 17.968707  1.809326 19.778032
#      TRUE  37.095070 31.463586 68.558656
#      Sum   55.063777 33.272911 88.336688
#
# $jxn
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 31.616470  1.451001 33.067471
#      TRUE  34.855760 20.471410 55.327170
#      Sum   66.472230 21.922410 88.394640
#
# $tx
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 11.596255  2.142042 13.738296
#      TRUE  35.875771 40.168988 76.044759
#      Sum   47.472026 42.311030 89.783055


comp_qtl_short(hippo)
# 2018-09-18 21:55:25 removing some NAs from CAUC (TRUEs below) for type gene
#
#  FALSE   TRUE
# 937821 137912
# 2018-09-18 21:55:27 removing some NAs from CAUC (TRUEs below) for type exon
#
#   FALSE    TRUE
# 5425045  781330
# 2018-09-18 21:55:38 removing some NAs from CAUC (TRUEs below) for type jxn
#
#   FALSE    TRUE
# 3478432  476817
# 2018-09-18 21:55:44 removing some NAs from CAUC (TRUEs below) for type tx
#
#   FALSE    TRUE
# 1497731  204844
# $gene
#           CAUC p<0.01
# Equal sign  FALSE   TRUE    Sum
#      FALSE 171170  15204 186374
#      TRUE  518619 232828 751447
#      Sum   689789 248032 937821
#
# $exon
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE 1209677   82436 1292113
#      TRUE  2947392 1185540 4132932
#      Sum   4157069 1267976 5425045
#
# $jxn
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE 1058807   29985 1088792
#      TRUE  1672659  716981 2389640
#      Sum   2731466  746966 3478432
#
# $tx
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE  284771   19350  304121
#      TRUE   808193  385417 1193610
#      Sum   1092964  404767 1497731

comp_qtl_short(hippo, perc = TRUE)
# $gene
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 15.911941  1.413362 17.325303
#      TRUE  48.210755 21.643661 69.854416
#      Sum   64.122696 23.057023 87.179718
#
# $exon
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 19.490878  1.328247 20.819125
#      TRUE  47.489750 19.101972 66.591722
#      Sum   66.980629 20.430219 87.410848
#
# $jxn
#           CAUC p<0.01
# Equal sign      FALSE       TRUE        Sum
#      FALSE 26.7696673  0.7581065 27.5277739
#      TRUE  42.2896005 18.1273290 60.4169295
#      Sum   69.0592678 18.8854355 87.9447034
#
# $tx
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 16.725900  1.136514 17.862414
#      TRUE  47.468863 22.637299 70.106163
#      Sum   64.194764 23.773813 87.968577

comp_qtl_short(hippo, cutde = 0.05)
# $gene
#           CAUC p<0.01
# Equal sign  FALSE   TRUE    Sum
#      FALSE 160756  25618 186374
#      TRUE  408772 342675 751447
#      Sum   569528 368293 937821
#
# $exon
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE 1140304  151809 1292113
#      TRUE  2400563 1732369 4132932
#      Sum   3540867 1884178 5425045
#
# $jxn
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE 1029699   59093 1088792
#      TRUE  1333462 1056178 2389640
#      Sum   2363161 1115271 3478432
#
# $tx
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE  271088   33033  304121
#      TRUE   645543  548067 1193610
#      Sum    916631  581100 1497731

comp_qtl_short(hippo, perc = TRUE, cutde = 0.05)
# $gene
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 14.943857  2.381446 17.325303
#      TRUE  37.999392 31.855024 69.854416
#      Sum   52.943249 34.236469 87.179718
#
# $exon
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 18.373108  2.446017 20.819125
#      TRUE  38.678987 27.912735 66.591722
#      Sum   57.052096 30.358752 87.410848
#
# $jxn
#           CAUC p<0.01
# Equal sign    FALSE     TRUE      Sum
#      FALSE 26.03373  1.49404 27.52777
#      TRUE  33.71373 26.70320 60.41693
#      Sum   59.74746 28.19724 87.94470
#
# $tx
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 15.922235  1.940179 17.862414
#      TRUE  37.915687 32.190476 70.106163
#      Sum   53.837922 34.130655 87.968577

comp_qtl_short(dlpfc)
# 2018-09-18 21:56:41 removing some NAs from CAUC (TRUEs below) for type gene
#
#   FALSE    TRUE
# 1371570  206393
# 2018-09-18 21:56:44 removing some NAs from CAUC (TRUEs below) for type exon
#
#   FALSE    TRUE
# 7805138 1123171
# 2018-09-18 21:57:00 removing some NAs from CAUC (TRUEs below) for type jxn
#
#   FALSE    TRUE
# 4630635  629491
# 2018-09-18 21:57:09 removing some NAs from CAUC (TRUEs below) for type tx
#
#   FALSE    TRUE
# 1994482  276782
# $gene
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE  244982   22214  267196
#      TRUE   750876  353498 1104374
#      Sum    995858  375712 1371570
#
# $exon
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE 1701101  140661 1841762
#      TRUE  4187885 1775491 5963376
#      Sum   5888986 1916152 7805138
#
# $jxn
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE 1277449   47132 1324581
#      TRUE  2248301 1057753 3306054
#      Sum   3525750 1104885 4630635
#
# $tx
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE  348925   28295  377220
#      TRUE  1044562  572700 1617262
#      Sum   1393487  600995 1994482
#

comp_qtl_short(dlpfc, perc = TRUE)
# $gene
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 15.525206  1.407764 16.932970
#      TRUE  47.585146 22.402173 69.987319
#      Sum   63.110352 23.809937 86.920289
#
# $exon
#           CAUC p<0.01
# Equal sign    FALSE     TRUE      Sum
#      FALSE 19.05289  1.57545 20.62834
#      TRUE  46.90569 19.88608 66.79177
#      Sum   65.95858 21.46153 87.42012
#
# $jxn
#           CAUC p<0.01
# Equal sign      FALSE       TRUE        Sum
#      FALSE 24.2855209  0.8960242 25.1815451
#      TRUE  42.7423412 20.1088909 62.8512321
#      Sum   67.0278621 21.0049151 88.0327772
#
# $tx
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 15.362591  1.245782 16.608373
#      TRUE  45.990338 25.215034 71.205373
#      Sum   61.352929 26.460817 87.813746

comp_qtl_short(dlpfc, cutde = 0.05)
# $gene
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE  231205   35991  267196
#      TRUE   592282  512092 1104374
#      Sum    823487  548083 1371570
#
# $exon
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE 1584302  257460 1841762
#      TRUE  3374285 2589091 5963376
#      Sum   4958587 2846551 7805138
#
# $jxn
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE 1246207   78374 1324581
#      TRUE  1799967 1506087 3306054
#      Sum   3046174 1584461 4630635
#
# $tx
#           CAUC p<0.01
# Equal sign   FALSE    TRUE     Sum
#      FALSE  329602   47618  377220
#      TRUE   827367  789895 1617262
#      Sum   1156969  837513 1994482

comp_qtl_short(dlpfc, perc = TRUE, cutde = 0.05)
# $gene
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 14.652118  2.280852 16.932970
#      TRUE  37.534594 32.452725 69.987319
#      Sum   52.186712 34.733577 86.920289
#
# $exon
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 17.744704  2.883637 20.628341
#      TRUE  37.793103 28.998672 66.791774
#      Sum   55.537807 31.882308 87.420115
#
# $jxn
#           CAUC p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 23.691581  1.489964 25.181545
#      TRUE  34.219085 28.632147 62.851232
#      Sum   57.910666 30.122111 88.032777
#
# $tx
#           CAUC p<0.01
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

