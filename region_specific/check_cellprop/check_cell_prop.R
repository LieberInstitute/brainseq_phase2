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
    f <- paste0('../rda/limma_dev_interaction_adjCellProp_', type, '.Rdata')
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
dev_type <- split(dev, factor(dev$type, levels = features))

tab_pbonf <- map(
    dev_type,
    ~ with(.x, table('Original Bonf<1%' = P.Bonf < 0.01, 'CellProp Bonf<1%' = neun_P.Bonf < 0.01))
)
map(tab_pbonf, addmargins)
# $gene
#                 CellProp Bonf<1%
# Original Bonf<1% FALSE  TRUE   Sum
#            FALSE  3878  3029  6907
#            TRUE   1382 16363 17745
#            Sum    5260 19392 24652
#
# $exon
#                 CellProp Bonf<1%
# Original Bonf<1%  FALSE   TRUE    Sum
#            FALSE  67292  42452 109744
#            TRUE   24156 262679 286835
#            Sum    91448 305131 396579
#
# $jxn
#                 CellProp Bonf<1%
# Original Bonf<1%  FALSE   TRUE    Sum
#            FALSE  85447  19347 104794
#            TRUE   26737 165650 192387
#            Sum   112184 184997 297181
#
# $tx
#                 CellProp Bonf<1%
# Original Bonf<1% FALSE  TRUE   Sum
#            FALSE 88040   877 88917
#            TRUE   1049  2766  3815
#            Sum   89089  3643 92732
map_dbl(tab_pbonf, getOR)
#     gene      exon       jxn        tx
# 15.15875  17.23716  27.36289 264.70194
map_dbl(tab_pbonf, ~ chisq.test(.x)$p.value)
# gene exon  jxn   tx
#    0    0    0    0


tab_pbonf_span <- map(
    dev_type,
    ~ with(.x, table('Original Bonf<1% & Rep BrainSpan' = P.Bonf < 0.01 & span_P.Value < 0.05, 'CellProp Bonf<1% & Rep BrainSpan' = neun_P.Bonf < 0.01  & span_P.Value < 0.05))
)
map(tab_pbonf_span , addmargins)
# $gene
#                                 CellProp Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan FALSE  TRUE   Sum
#                            FALSE 12501  1312 13813
#                            TRUE    748 10091 10839
#                            Sum   13249 11403 24652
#
# $exon
#                                 CellProp Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan  FALSE   TRUE    Sum
#                            FALSE 206644  20682 227326
#                            TRUE   13333 155920 169253
#                            Sum   219977 176602 396579
#
# $jxn
#                                 CellProp Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan  FALSE   TRUE    Sum
#                            FALSE 139690  13596 153286
#                            TRUE   21195 122700 143895
#                            Sum   160885 136296 297181
#
# $tx
#                                 CellProp Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan FALSE  TRUE   Sum
#                            FALSE 90783   234 91017
#                            TRUE    408  1307  1715
#                            Sum   91191  1541 92732
map_dbl(tab_pbonf_span , getOR)
#      gene       exon        jxn         tx
# 128.54155  116.84341   59.47923 1242.80816
map_dbl(tab_pbonf_span , ~ chisq.test(.x)$p.value)
# gene exon  jxn   tx
#    0    0    0    0

make_table <- function(ov) {
    
    res <- map_dfr(ov,
        ~ as.data.frame(matrix(as.vector(.x[1:2, 1:2]), nrow = 1, dimnames = list(1, c('Null_both', 'Original_only', 'CellProp_only', 'Both'))))
    )
    res$feature <- names(ov)
    res$OR <- map_dbl(ov, ~ getOR(.x[1:2, 1:2]))
    res$pval <- map_dbl(ov, ~ chisq.test(.x[1:2, 1:2])$p.value)
    res$pval_bonf <- p.adjust(res$pval, 'bonf')
    return(res)
}

options(width = 120)
make_table(tab_pbonf)
#   Null_both Original_only CellProp_only   Both feature        OR pval pval_bonf
# 1      3878          1382          3029  16363    gene  15.15875    0         0
# 2     67292         24156         42452 262679    exon  17.23716    0         0
# 3     85447         26737         19347 165650     jxn  27.36289    0         0
# 4     88040          1049           877   2766      tx 264.70194    0         0
make_table(tab_pbonf_span)
#   Null_both Original_only CellProp_only   Both feature         OR pval pval_bonf
# 1     12501           748          1312  10091    gene  128.54155    0         0
# 2    206644         13333         20682 155920    exon  116.84341    0         0
# 3    139690         21195         13596 122700     jxn   59.47923    0         0
# 4     90783           408           234   1307      tx 1242.80816    0         0

map_dbl(dev_type, ~ cor(.x$F, .x$neun_F))
#      gene      exon       jxn        tx
# 0.8702397 0.8715382 0.9115000 0.9275750

## Compute the correlation on the scale that I'm actually plotting below
map_dbl(dev_type, ~ cor(log10(.x$F), log10(.x$neun_F)))
#      gene      exon       jxn        tx
# 0.8227183 0.8367599 0.8740846 0.8821472

corrs <- cbind(map_dfr(dev_type, ~ map_dbl(split(.x, .x$de), ~ cor(log10(.x$F), log10(.x$neun_F)))), DE = c('FALSE', 'TRUE'))
corrs
#        gene      exon       jxn        tx    DE
# 1 0.7959842 0.8225326 0.8574005 0.8673563 FALSE
# 2 0.7670905 0.7819954 0.8324641 0.8729942  TRUE


pdf('f_original_vs_f_adjCellProp_by_feature.pdf', width = 12, useDingbats = FALSE)
map2(dev_type, names(dev_type), function(df, type) {
    print(
        ggplot(df, aes(x = F, y = neun_F)) +
            geom_hex(aes(fill=..density..), bins = 100) +
            scale_x_log10() +
            scale_y_log10() +
            facet_grid( ~ de) +
            theme_bw(base_size = 30) +
            xlab('F-statistic: original') + 
            ylab('F-statistic: adj. cell fraction') +
            labs(caption = 'Separated by DE status', title =  paste(type, 'corr =', paste(signif(corrs[, type], 3), collapse = ', ')))
    )
    return(invisible(NULL))
})
dev.off()

pdf('f_original_vs_f_adjCellProp.pdf', width = 12, useDingbats = FALSE, height = 18)
ggplot(dev, aes(x = F, y = neun_F)) +
    geom_hex(aes(fill=..density..), bins = 100) +
    scale_x_log10() +
    scale_y_log10() +
    facet_grid(type ~ de) +
    theme_bw(base_size = 30) +
    xlab('F-statistic: original') + 
    ylab('F-statistic: adj. cell fraction') +
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
#  date     2019-03-25
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
