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
#            FALSE  4566  2341  6907
#            TRUE   1241 16504 17745
#            Sum    5807 18845 24652
#
# $exon
#                 NeuN Bonf<1%
# Original Bonf<1%  FALSE   TRUE    Sum
#            FALSE  77518  32226 109744
#            TRUE   19697 267138 286835
#            Sum    97215 299364 396579
#
# $jxn
#                 NeuN Bonf<1%
# Original Bonf<1%  FALSE   TRUE    Sum
#            FALSE  89159  15635 104794
#            TRUE   15911 176476 192387
#            Sum   105070 192111 297181
#
# $tx
#                 NeuN Bonf<1%
# Original Bonf<1% FALSE  TRUE   Sum
#            FALSE 88313   604 88917
#            TRUE    662  3153  3815
#            Sum   88975  3757 92732
map_dbl(tab_pbonf, getOR)
#     gene      exon       jxn        tx
# 25.93892  32.62359  63.24926 696.39185
map_dbl(tab_pbonf, ~ chisq.test(.x)$p.value)
# gene exon  jxn   tx
#    0    0    0    0


tab_pbonf_span <- map(
    dev_type,
    ~ with(.x, table('Original Bonf<1% & Rep BrainSpan' = P.Bonf < 0.01 & span_P.Value < 0.05, 'NeuN Bonf<1% & Rep BrainSpan' = neun_P.Bonf < 0.01  & span_P.Value < 0.05))
)
map(tab_pbonf_span , addmargins)
# $gene
#                                 NeuN Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan FALSE  TRUE   Sum
#                            FALSE 12743  1070 13813
#                            TRUE    647 10192 10839
#                            Sum   13390 11262 24652
#
# $exon
#                                 NeuN Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan  FALSE   TRUE    Sum
#                            FALSE 211284  16042 227326
#                            TRUE   10722 158531 169253
#                            Sum   222006 174573 396579
#
# $jxn
#                                 NeuN Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan  FALSE   TRUE    Sum
#                            FALSE 142089  11197 153286
#                            TRUE   12034 131861 143895
#                            Sum   154123 143058 297181
#
# $tx
#                                 NeuN Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan FALSE  TRUE   Sum
#                            FALSE 90846   171 91017
#                            TRUE    224  1491  1715
#                            Sum   91070  1662 92732
map_dbl(tab_pbonf_span , getOR)
#     gene      exon       jxn        tx
# 187.6044  194.7361  139.0481 3536.2204
map_dbl(tab_pbonf_span , ~ chisq.test(.x)$p.value)
# gene exon  jxn   tx
#    0    0    0    0

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
#   Null_both Original_only NeuN_only   Both feature        OR pval pval_bonf
# 1      4566          1241      2341  16504    gene  25.93892    0         0
# 2     77518         19697     32226 267138    exon  32.62359    0         0
# 3     89159         15911     15635 176476     jxn  63.24926    0         0
# 4     88313           662       604   3153      tx 696.39185    0         0
make_table(tab_pbonf_span)
#   Null_both Original_only NeuN_only   Both feature        OR pval pval_bonf
# 1     12743           647      1070  10192    gene  187.6044    0         0
# 2    211284         10722     16042 158531    exon  194.7361    0         0
# 3    142089         12034     11197 131861     jxn  139.0481    0         0
# 4     90846           224       171   1491      tx 3536.2204    0         0

map_dbl(dev_type, ~ cor(.x$F, .x$neun_F))
#      gene      exon       jxn        tx
# 0.9088368 0.9150995 0.9132362 0.9671942

## Compute the correlation on the scale that I'm actually plotting below
map_dbl(dev_type, ~ cor(log10(.x$F), log10(.x$neun_F)))
#      gene      exon       jxn        tx
# 0.8823008 0.8932323 0.9192852 0.9490049

corrs <- cbind(map_dfr(dev_type, ~ map_dbl(split(.x, .x$de), ~ cor(log10(.x$F), log10(.x$neun_F)))), DE = c('FALSE', 'TRUE'))
corrs
#        gene      exon       jxn        tx    DE
# 1 0.8517315 0.8760462 0.9010913 0.9420888 FALSE
# 2 0.8547177 0.8640188 0.8863927 0.9478907  TRUE


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
#  date     2019-03-19
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package          * version   date       lib source
#  assertthat         0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  BiocGenerics       0.28.0    2018-10-30 [1] Bioconductor
#  bitops             1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  cli                1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout         * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace         1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon             1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  digest             0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr              0.8.0.1   2019-02-15 [1] CRAN (R 3.5.1)
#  GenomeInfoDb       1.18.2    2019-02-12 [1] Bioconductor
#  GenomeInfoDbData   1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges      1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2          * 3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue               1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
#  gtable             0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  hexbin           * 1.27.2    2018-01-15 [2] CRAN (R 3.5.0)
#  htmltools          0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets        1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv             1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges            2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab         * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  jsonlite           1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  labeling           0.3       2014-08-23 [2] CRAN (R 3.5.0)
#  later              0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice            0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval           0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma              3.38.3    2018-12-02 [1] Bioconductor
#  magrittr           1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  munsell            0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  pillar             1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig          2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr               1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises           1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr            * 0.3.1     2019-03-03 [2] CRAN (R 3.5.1)
#  R6                 2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  rafalib          * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer       1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp               1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl              1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  reshape2           1.4.3     2017-12-11 [2] CRAN (R 3.5.0)
#  rlang              0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote            * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors          0.20.1    2018-11-09 [1] Bioconductor
#  scales             1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented          0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr              0.13      2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo      * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  stringi            1.4.3     2019-03-12 [2] CRAN (R 3.5.1)
#  stringr            1.4.0     2019-02-10 [1] CRAN (R 3.5.1)
#  tibble             2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect         0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr              2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun               0.5       2019-02-20 [1] CRAN (R 3.5.1)
#  XVector            0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc           1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
