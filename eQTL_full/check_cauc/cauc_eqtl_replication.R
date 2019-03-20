library('data.table')
library('sessioninfo')
library('SummarizedExperiment')

setDTthreads(1)

## Load subsets of data
files_sub <- dir('rdas', pattern = '_compare_|_sig_genes_', full.names = TRUE)
stopifnot(length(files_sub) == 6)

## If resuming the job:
if(FALSE) {
    files_sub <- dir('rdas', pattern = 'merged_CAUC_BrainSeq_QTLs_', full.names = TRUE)
    stopifnot(length(files_sub) == 4)
}

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
    
    res <- addmargins(table(sign(df$statistic) == sign(df$cauc_statistic),  df$cauc_pvalue < cutde, dnn = c('Equal sign', paste0('CAUC p<', cutde)), useNA = 'ifany'))
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

comp_qtl_short(interaction, perc = TRUE, cutde = 0.001)
# $gene
#           CAUC p<0.001
# Equal sign      FALSE       TRUE        Sum
#      FALSE   4.439788   0.000000   4.439788
#      TRUE   37.690811  57.869400  95.560212
#      Sum    42.130600  57.869400 100.000000

comp_qtl_short(interaction, perc = TRUE, cutde = 0.0001)
# $gene
#           CAUC p<1e-04
# Equal sign      FALSE       TRUE        Sum
#      FALSE   4.439788   0.000000   4.439788
#      TRUE   57.340617  38.219595  95.560212
#      Sum    61.780405  38.219595 100.000000

comp_qtl_short(interaction, perc = TRUE, cutde = max(interaction$gene$pvalue))
# $gene
#           CAUC p<1.30632983923964e-05
# Equal sign      FALSE       TRUE        Sum
#      FALSE   4.439788   0.000000   4.439788
#      TRUE   67.901327  27.658885  95.560212
#      Sum    72.341115  27.658885 100.000000

or_chisq <- function(x) {
    x <- x[1:2, 1:2]
    list(
        'OR' = jaffelab::getOR(x),
        'chisq.test' = chisq.test(x),
        'p.value' = chisq.test(x)$p.value
    )
}

or_chisq(comp_qtl_short(interaction, cutde = 0.05)[[1]])
# $OR
# [1] Inf
#
# $chisq.test
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  x
# X-squared = 12684, df = 1, p-value < 2.2e-16
#
#
# $p.value
# [1] 0


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

or_chisq(comp_qtl_short(hippo, cutde = 0.05)[[1]])
# $OR
# [1] 914.8403
#
# $chisq.test
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  x
# X-squared = 304940, df = 1, p-value < 2.2e-16
#
#
# $p.value
# [1] 0

comp_qtl_short(hippo, perc = TRUE, cutde = 0.001)
# $gene
#           CAUC p<0.001
# Equal sign        FALSE         TRUE          Sum
#      FALSE 5.531670e+00 5.391672e-03 5.537062e+00
#      TRUE  3.535245e+01 5.911049e+01 9.446294e+01
#      Sum   4.088412e+01 5.911588e+01 1.000000e+02

comp_qtl_short(hippo, perc = TRUE, cutde = 0.0001)
# $gene
#           CAUC p<1e-04
# Equal sign        FALSE         TRUE          Sum
#      FALSE 5.534180e+00 2.881756e-03 5.537062e+00
#      TRUE  4.913347e+01 4.532946e+01 9.446294e+01
#      Sum   5.466765e+01 4.533235e+01 1.000000e+02

comp_qtl_short(hippo, perc = TRUE, cutde = max(hippo$gene$pvalue))
# $gene
#           CAUC p<0.000182528411297916
# Equal sign        FALSE         TRUE          Sum
#      FALSE 5.534087e+00 2.974716e-03 5.537062e+00
#      TRUE  4.592803e+01 4.853491e+01 9.446294e+01
#      Sum   5.146212e+01 4.853788e+01 1.000000e+02




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

or_chisq(comp_qtl_short(dlpfc, cutde = 0.05)[[1]])
# $OR
# [1] 865.3454
#
# $chisq.test
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  x
# X-squared = 454150, df = 1, p-value < 2.2e-16
#
#
# $p.value
# [1] 0

comp_qtl_short(dlpfc, perc = TRUE, cutde = 0.001)
# $gene
#           CAUC p<0.001
# Equal sign        FALSE         TRUE          Sum
#      FALSE   5.89779355   0.00557681   5.90337036
#      TRUE   36.60459719  57.49203245  94.09662964
#      Sum    42.50239074  57.49760926 100.00000000

comp_qtl_short(dlpfc, perc = TRUE, cutde = 0.0001)
# $gene
#           CAUC p<1e-04
# Equal sign        FALSE         TRUE          Sum
#      FALSE 5.900265e+00 3.105269e-03 5.903370e+00
#      TRUE  5.029123e+01 4.380540e+01 9.409663e+01
#      Sum   5.619149e+01 4.380851e+01 1.000000e+02

comp_qtl_short(dlpfc, perc = TRUE, cutde = max(dlpfc$gene$pvalue))
# $gene
#           CAUC p<0.000267745820149456
# Equal sign        FALSE         TRUE          Sum
#      FALSE 5.899948e+00 3.422133e-03 5.903370e+00
#      TRUE  4.478451e+01 4.931212e+01 9.409663e+01
#      Sum   5.068446e+01 4.931554e+01 1.000000e+02


## For testing
# cauc <- i_sig_genes_cauc
# brainseq <- i_sig_genes

## This is the reverse of subset_cauc() defined
## in the results_subset_*.R scripts
##
## subset BrainSeq based on CAUC
subset_brainseq <- function(cauc, brainseq) {    
    message(paste(Sys.time(), 'create keys: brainseq'))
    setkey(brainseq, snps, gene)

    message(paste(Sys.time(), 'subset brainseq by cauc'))
    brainseq[.(cauc$snps, cauc$gene)]
}


merge_qtl_rev <- function(cauc, brainseq) {
    ## Match CAUC to BrainSeq
    brainseq <- subset_brainseq(cauc, brainseq)
    
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

## Merge (reverse)
message(paste(Sys.time(), 'merging interaction QTLs (from CAUC)'))
cauc_sig <- list(
    'interaction' = merge_qtl_rev(i_sig_genes_cauc, i_sig_genes),
    'dlpfc' = merge_qtl_rev(d_sig_genes_cauc, d_sig_genes),
    'hippo' = merge_qtl_rev(h_sig_genes_cauc, h_sig_genes)
)
sapply(cauc_sig, nrow)
# interaction       dlpfc       hippo
#       11062      844726      573296

## How many of the FDR<1% CAUC eQTLs are FDR<1% BrainSeq Phase II eQTLs
add_prefix <- function(x, prefix) {
    y <- t(x)
    colnames(y) <- paste0(prefix, '_', rownames(x))
    return(y)
}
cauc_sig_in_brainseq <- sapply(cauc_sig, function(x) table(!is.na(x$statistic)))
cauc_sig_in_brainseq_perc <- sweep(cauc_sig_in_brainseq, 2, 1 / sapply(cauc_sig, nrow) * 100, '*')

(cauc_sig_in_brainseq_tab <- cbind(add_prefix(cauc_sig_in_brainseq, 'N'), add_prefix(cauc_sig_in_brainseq_perc, 'Percent')))
#             N_FALSE N_TRUE Percent_FALSE Percent_TRUE
# interaction    2393   8669      21.63262     78.36738
# dlpfc        124627 720099      14.75354     85.24646
# hippo         87323 485973      15.23175     84.76825

message(paste(Sys.time(), 'saving merged QTLs (from CAUC)'))
save(cauc_sig, file = 'rdas/merged_CAUC_BrainSeq_QTLs_CAUC_sig.Rdata')



comp_qtl_rev <- function(type, dfs, perc = FALSE) {
    df <- dfs[[type]]
    
    res <- addmargins(table(sign(df$statistic) == sign(df$cauc_statistic),  !is.na(df$statistic), dnn = c('Equal sign', 'BrainSeq FDR<1%'), useNA = 'ifany'))
    if(!perc) return(res)
    
    ## Calculate percent over all of brainseq
    ## the total marginal will not be 100% unless there were no NAs
    round(res / nrow(df) * 100, 5)
}
comp_qtl_short_rev <- function(dfs, perc = FALSE) {
    res <- lapply(names(dfs), comp_qtl_rev, dfs = dfs, perc = perc)
    names(res) <- names(dfs)
    return(res)
}



## These run, but are not informative
# comp_qtl_short(cauc_sig, cutde = 0.05)
# comp_qtl_short(cauc_sig, cutde = 0.05, perc = TRUE)

comp_qtl_short_rev(cauc_sig)
# $interaction
#           BrainSeq FDR<1%
# Equal sign FALSE  TRUE   Sum
#       TRUE     0  8669  8669
#       <NA>  2393     0  2393
#       Sum   2393  8669 11062
#
# $dlpfc
#           BrainSeq FDR<1%
# Equal sign  FALSE   TRUE    Sum
#      FALSE      0     50     50
#      TRUE       0 720049 720049
#      <NA>  124627      0 124627
#      Sum   124627 720099 844726
#
# $hippo
#           BrainSeq FDR<1%
# Equal sign  FALSE   TRUE    Sum
#      FALSE      0     31     31
#      TRUE       0 485942 485942
#      <NA>   87323      0  87323
#      Sum    87323 485973 573296
comp_qtl_short_rev(cauc_sig, perc = TRUE)
# $interaction
#           BrainSeq FDR<1%
# Equal sign     FALSE      TRUE       Sum
#       TRUE   0.00000  78.36738  78.36738
#       <NA>  21.63262   0.00000  21.63262
#       Sum   21.63262  78.36738 100.00000
#
# $dlpfc
#           BrainSeq FDR<1%
# Equal sign     FALSE      TRUE       Sum
#      FALSE   0.00000   0.00592   0.00592
#      TRUE    0.00000  85.24054  85.24054
#      <NA>   14.75354   0.00000  14.75354
#      Sum    14.75354  85.24646 100.00000
#
# $hippo
#           BrainSeq FDR<1%
# Equal sign     FALSE      TRUE       Sum
#      FALSE   0.00000   0.00541   0.00541
#      TRUE    0.00000  84.76285  84.76285
#      <NA>   15.23175   0.00000  15.23175
#      Sum    15.23175  84.76825 100.00000

## Percet of those that are FDR<1% in both that have either unequal t-stat signs (first row) or equal signs
e_sign_perc <- sapply(comp_qtl_short_rev(cauc_sig)[2:3], function(x) { y  <- x[1, 2] / x[4, 2] * 100 ; return(c(y, 100 - y))})
rownames(e_sign_perc) <- c('FALSE', 'TRUE')
e_sign_perc
#             dlpfc        hippo
# FALSE  0.00694349  0.006378955
# TRUE  99.99305651 99.993621045


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
#  date     2019-03-19
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.6    2019-02-10 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
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
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab               0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  later                  0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma                  3.38.3    2018-12-02 [1] Bioconductor
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                  0.3.1     2019-03-03 [2] CRAN (R 3.5.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  rafalib                1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
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
