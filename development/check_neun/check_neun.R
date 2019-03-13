library('sessioninfo')
library('purrr')
library('jaffelab')
library('ggplot2')

message(paste(Sys.time(), 'loading ../rda/pcheck_both.Rdata'))
load('../rda/pcheck_both.Rdata', verbose = TRUE)

get_de <- function(x) {
    sign(x$F) == sign(x$span_F) & x$span_P.Value < 0.05 & x$P.Bonf < 0.01
}
pcheck_both$de <- get_de(pcheck_both)

## Rename for simplicity
dev <- pcheck_both

features <- c('gene', 'exon', 'jxn', 'tx')
top <- lapply(features, function(type) {
    f <- paste0('../rda/limma_dev_interaction_adjNeunProp_', type, '.Rdata')
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    top$type <- type
    top$P.Bonf <- p.adjust(top$P.Value, 'bonf')
    return(top)
})
names(top) <- features


neun <- do.call(rbind, map(top, function(x) {
    colnames(x) <- paste0('neun_', colnames(x))
    x
}))

stopifnot(identical(nrow(dev), nrow(neun)))

m <- match(rownames(dev), rownames(neun))
stopifnot(!any(is.na(m)))
neun <- neun[m, ]

dev <- cbind(dev, neun)
dev_type <- split(dev, factor(dev$type, levels = c('gene', 'exon', 'jxn', 'tx')))

tab_pbonf <- map(
    dev_type,
    ~ with(.x, table('Original Bonf<1%' = P.Bonf < 0.01, 'NeuN Bonf<1%' = neun_P.Bonf < 0.01))
)
map(tab_pbonf, addmargins)
# $gene
#                 NeuN Bonf<1%
# Original Bonf<1% FALSE  TRUE   Sum
#            FALSE  6905     2  6907
#            TRUE  17061   684 17745
#            Sum   23966   686 24652
#
# $exon
#                 NeuN Bonf<1%
# Original Bonf<1%  FALSE   TRUE    Sum
#            FALSE 109716     28 109744
#            TRUE  283760   3075 286835
#            Sum   393476   3103 396579
#
# $jxn
#                 NeuN Bonf<1%
# Original Bonf<1%  FALSE   TRUE    Sum
#            FALSE 104732     62 104794
#            TRUE  191304   1083 192387
#            Sum   296036   1145 297181
#
# $tx
#                 NeuN Bonf<1%
# Original Bonf<1% FALSE  TRUE   Sum
#            FALSE 88902    15 88917
#            TRUE   3345   470  3815
#            Sum   92247   485 92732
map_dbl(tab_pbonf, getOR)
#       gene       exon        jxn         tx
# 138.415685  42.462531   9.562955 832.764126
map_dbl(tab_pbonf, ~ chisq.test(.x)$p.value)
#         gene          exon           jxn            tx
# 3.859747e-60 3.197824e-245  2.838502e-99  0.000000e+00


tab_pbonf_span <- map(
    dev_type,
    ~ with(.x, table('Original Bonf<1% & Rep BrainSpan' = P.Bonf < 0.01 & span_P.Value < 0.05, 'NeuN Bonf<1% & Rep BrainSpan' = neun_P.Bonf < 0.01  & span_P.Value < 0.05))
)
map(tab_pbonf_span , addmargins)
# $gene
#                                 NeuN Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan FALSE  TRUE   Sum
#                            FALSE 13813     0 13813
#                            TRUE  10293   546 10839
#                            Sum   24106   546 24652
#
# $exon
#                                 NeuN Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan  FALSE   TRUE    Sum
#                            FALSE 227311     15 227326
#                            TRUE  166744   2509 169253
#                            Sum   394055   2524 396579
#
# $jxn
#                                 NeuN Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan  FALSE   TRUE    Sum
#                            FALSE 153227     59 153286
#                            TRUE  142957    938 143895
#                            Sum   296184    997 297181
#
# $tx
#                                 NeuN Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan FALSE  TRUE   Sum
#                            FALSE 91012     5 91017
#                            TRUE   1415   300  1715
#                            Sum   92427   305 92732
map_dbl(tab_pbonf_span , getOR)
# gene       exon        jxn         tx
#  Inf  228.02352   17.04044 3859.16608
map_dbl(tab_pbonf_span , ~ chisq.test(.x)$p.value)
#          gene          exon           jxn            tx
# 2.916268e-156  0.000000e+00 3.084074e-183  0.000000e+00

make_table <- function(ov) {
    
    res <- map_dfr(ov,
        ~ as.data.frame(matrix(as.vector(.x[1:2, 1:2]), nrow = 1, dimnames = list(1, c('Null_both', 'Original_only', 'NeuN_only', 'Both'))))
    )
    res$feature <- names(ov)
    res$OR <- map_dbl(ov, ~ getOR(.x[1:2, 1:2]))
    res$pval <- map_dbl(ov, ~ chisq.test(.x[1:2, 1:2])$p.value)
    res$pval_bonf <- p.adjust(res$pval, 'bonf')
    return(res)
}

options(width = 120)
make_table(tab_pbonf)
#   Null_both Original_only NeuN_only Both feature         OR          pval     pval_bonf
# 1      6905         17061         2  684    gene 138.415685  3.859747e-60  1.543899e-59
# 2    109716        283760        28 3075    exon  42.462531 3.197824e-245 1.279130e-244
# 3    104732        191304        62 1083     jxn   9.562955  2.838502e-99  1.135401e-98
# 4     88902          3345        15  470      tx 832.764126  0.000000e+00  0.000000e+00
make_table(tab_pbonf_span)
  Null_both Original_only NeuN_only Both feature         OR          pval     pval_bonf
# 1     13813         10293         0  546    gene        Inf 2.916268e-156 1.166507e-155
# 2    227311        166744        15 2509    exon  228.02352  0.000000e+00  0.000000e+00
# 3    153227        142957        59  938     jxn   17.04044 3.084074e-183 1.233629e-182
# 4     91012          1415         5  300      tx 3859.16608  0.000000e+00  0.000000e+00

map_dbl(dev_type, ~ cor(.x$F, .x$neun_F))
#       gene        jxn         tx
# 0.08688325 0.16308335 0.70752450

## Compute the correlation on the scale that I'm actually plotting below
map_dbl(dev_type, ~ cor(log10(.x$F), log10(.x$neun_F)))
#       gene       exon        jxn         tx
# 0.08688325 0.09072702 0.16308335 0.70752450

corrs <- cbind(map_dfr(dev_type, ~ map_dbl(split(.x, .x$de), ~ cor(log10(.x$F), log10(.x$neun_F)))), DE = c('FALSE', 'TRUE'))
corrs
#         gene      exon       jxn        tx    DE
# 1 0.16026235 0.1181208 0.3079739 0.5696692 FALSE
# 2 0.06215385 0.0578365 0.1700112 0.5925908  TRUE


# ggplot(dev_type$gene, aes(x = F, y = neun_F, color = de)) + geom_point() + scale_x_log10() + scale_y_log10() + facet_grid( ~ de)

pdf('f_original_vs_f_adjNeuN_by_feature.pdf', width = 12, useDingbats = FALSE)
map2(dev_type, names(dev_type), function(df, type) {
    print(
        ggplot(df, aes(x = F, y = neun_F)) +
            geom_hex(aes(fill=..density..), bins = 100) +
            scale_x_log10() +
            scale_y_log10() +
            facet_grid( ~ de) +
            theme_bw(base_size = 30) +
            xlab('F-statistic: original') + 
            ylab('F-statistic: adj. NeuN prop') +
            labs(caption = 'Separated by DE status', title =  paste(type, 'corr =', paste(signif(corrs[, type], 3), collapse = ', ')))
    )
    return(invisible(NULL))
})
dev.off()

pdf('f_original_vs_f_adjNeuN.pdf', width = 12, useDingbats = FALSE, height = 18)
ggplot(dev, aes(x = F, y = neun_F)) +
    geom_hex(aes(fill=..density..), bins = 100) +
    scale_x_log10() +
    scale_y_log10() +
    facet_grid(type ~ de) +
    theme_bw(base_size = 30) +
    xlab('F-statistic: original') + 
    ylab('F-statistic: adj. NeuN prop') +
    labs(caption = 'Separated by DE status')
dev.off()

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
#  date     2019-03-13
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package          * version   date       lib source
#  assertthat         0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  bindr              0.1.1     2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp           0.2.2     2018-03-29 [1] CRAN (R 3.5.0)
#  BiocGenerics       0.28.0    2018-10-30 [1] Bioconductor
#  bitops             1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  cli                1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout         * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace         1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon             1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  digest             0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr              0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
#  GenomeInfoDb       1.18.1    2018-11-12 [1] Bioconductor
#  GenomeInfoDbData   1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges      1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2          * 3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue               1.3.0     2018-07-17 [1] CRAN (R 3.5.1)
#  gtable             0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  hexbin           * 1.27.2    2018-01-15 [2] CRAN (R 3.5.0)
#  htmltools          0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets        1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv             1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges            2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab         * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  labeling           0.3       2014-08-23 [2] CRAN (R 3.5.0)
#  later              0.7.5     2018-09-18 [2] CRAN (R 3.5.1)
#  lattice            0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval           0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma              3.38.3    2018-12-02 [1] Bioconductor
#  magrittr           1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  munsell            0.5.0     2018-06-12 [2] CRAN (R 3.5.0)
#  pillar             1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig          2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr               1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises           1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr            * 0.2.5     2018-05-29 [2] CRAN (R 3.5.0)
#  R6                 2.3.0     2018-10-04 [2] CRAN (R 3.5.1)
#  rafalib          * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer       1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp               1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl              1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  reshape2           1.4.3     2017-12-11 [2] CRAN (R 3.5.0)
#  rlang              0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote            * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors          0.20.1    2018-11-09 [1] Bioconductor
#  scales             1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented          0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr              0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo      * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  stringi            1.2.4     2018-07-20 [2] CRAN (R 3.5.1)
#  stringr            1.3.1     2018-05-10 [1] CRAN (R 3.5.0)
#  tibble             2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect         0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr              2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun               0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  XVector            0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc           1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
