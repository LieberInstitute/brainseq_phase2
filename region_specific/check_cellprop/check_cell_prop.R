library('sessioninfo')
library('purrr')
library('jaffelab')
library('ggplot2')

message(paste(Sys.time(), 'loading ../rda/pcheck_both.Rdata'))
load('../rda/pcheck_both.Rdata', verbose = TRUE)

## Subset to only adults
dim(pcheck_both)
# [1] 1622288      20
pcheck_both <- subset(pcheck_both, age == 'adult')
dim(pcheck_both)
# [1] 811144     20

get_de <- function(x) {
    sign(x$t) == sign(x$span_t) & x$span_P.Value < 0.05 & x$P.Bonf < 0.01
}
pcheck_both$de <- get_de(pcheck_both)

## Rename for simplicity
ptab <- pcheck_both

features <- c('gene', 'exon', 'jxn', 'tx')
top <- lapply(features, function(type) {
    f <- paste0('../rda/limma_region_specific_adult_', type, '_adjCellProp.Rdata')
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


stopifnot(identical(nrow(ptab), nrow(neun)))

## Drop the age part from the names
rownames(ptab) <- gsub('adult_', '', rownames(ptab))

m <- match(rownames(ptab), rownames(neun))
stopifnot(!any(is.na(m)))
neun <- neun[m, ]

ptab <- cbind(ptab, neun)
ptab_type <- split(ptab, factor(ptab$type, levels = features))

tab_pbonf <- map(
    ptab_type,
    ~ with(.x, table('Original Bonf<1%' = P.Bonf < 0.01, 'CellProp Bonf<1%' = neun_P.Bonf < 0.01))
)
map(tab_pbonf, addmargins)
# $gene
#                 CellProp Bonf<1%
# Original Bonf<1% FALSE  TRUE   Sum
#            FALSE 15560  1197 16757
#            TRUE   2000  5895  7895
#            Sum   17560  7092 24652
#
# $exon
#                 CellProp Bonf<1%
# Original Bonf<1%  FALSE   TRUE    Sum
#            FALSE 301499  12744 314243
#            TRUE   24048  58288  82336
#            Sum   325547  71032 396579
#
# $jxn
#                 CellProp Bonf<1%
# Original Bonf<1%  FALSE   TRUE    Sum
#            FALSE 256768   4549 261317
#            TRUE   10988  24876  35864
#            Sum   267756  29425 297181
#
# $tx
#                 CellProp Bonf<1%
# Original Bonf<1% FALSE  TRUE   Sum
#            FALSE 76355  2337 78692
#            TRUE   4182  9858 14040
#            Sum   80537 12195 92732
map_dbl(tab_pbonf, getOR)
#     gene      exon       jxn        tx
# 38.31504  57.34299 127.78705  77.01646
map_dbl(tab_pbonf, ~ chisq.test(.x)$p.value)
# gene exon  jxn   tx
#    0    0    0    0


tab_pbonf_span <- map(
    ptab_type,
    ~ with(.x, table('Original Bonf<1% & Rep BrainSpan' = P.Bonf < 0.01 & span_P.Value < 0.05 & span_P.Value < 0.05  & sign(t) == sign(span_t), 'CellProp Bonf<1% & Rep BrainSpan' = neun_P.Bonf < 0.01  & span_P.Value < 0.05  & sign(neun_t) == sign(span_t)))
)
map(tab_pbonf_span , addmargins)
# $gene
#                                 CellProp Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan FALSE  TRUE   Sum
#                            FALSE 22844   196 23040
#                            TRUE    107  1505  1612
#                            Sum   22951  1701 24652
#
# $exon
#                                 CellProp Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan  FALSE   TRUE    Sum
#                            FALSE 379246   1891 381137
#                            TRUE    2048  13394  15442
#                            Sum   381294  15285 396579
#
# $jxn
#                                 CellProp Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan  FALSE   TRUE    Sum
#                            FALSE 290973    647 291620
#                            TRUE     887   4674   5561
#                            Sum   291860   5321 297181
#
# $tx
#                                 CellProp Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan FALSE  TRUE   Sum
#                            FALSE 90850   143 90993
#                            TRUE    372  1367  1739
#                            Sum   91222  1510 92732
map_dbl(tab_pbonf_span , getOR)
#     gene     exon      jxn       tx
# 1639.339 1311.625 2369.810 2334.611
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
#   Null_both Original_only CellProp_only  Both feature        OR pval pval_bonf
# 1     15560          2000          1197  5895    gene  38.31504    0         0
# 2    301499         24048         12744 58288    exon  57.34299    0         0
# 3    256768         10988          4549 24876     jxn 127.78705    0         0
# 4     76355          4182          2337  9858      tx  77.01646    0         0
make_table(tab_pbonf_span)
#   Null_both Original_only CellProp_only  Both feature       OR pval pval_bonf
# 1     22844           107           196  1505    gene 1639.339    0         0
# 2    379246          2048          1891 13394    exon 1311.625    0         0
# 3    290973           887           647  4674     jxn 2369.810    0         0
# 4     90850           372           143  1367      tx 2334.611    0         0

map_dbl(ptab_type, ~ cor(.x$t, .x$neun_t))
#      gene      exon       jxn        tx
# 0.9355817 0.9352117 0.9412695 0.9490653

corrs <- cbind(map_dfr(ptab_type, ~ map_dbl(split(.x, .x$de), ~ cor(.x$t, .x$neun_t))), DE = c('FALSE', 'TRUE'))
corrs
#        gene      exon       jxn        tx    DE
# 1 0.9136217 0.9213504 0.9335811 0.9420496 FALSE
# 2 0.9873258 0.9894706 0.9899145 0.9930099  TRUE


pdf('t_original_vs_t_adjCellProp_by_feature.pdf', width = 12, useDingbats = FALSE)
map2(ptab_type, names(ptab_type), function(df, type) {
    print(
        ggplot(df, aes(x = t, y = neun_t)) +
            geom_hex(aes(fill=..density..), bins = 100) +
            # scale_x_log10() +
            # scale_y_log10() +
            facet_grid( ~ de) +
            theme_bw(base_size = 30) +
            xlab('t-statistic: original') + 
            ylab('t-statistic: adj. cell fraction') +
            labs(caption = 'Separated by DE status', title =  paste(type, 'corr =', paste(signif(corrs[, type], 3), collapse = ', ')))
    )
    return(invisible(NULL))
})
dev.off()

pdf('t_original_vs_t_adjCellProp.pdf', width = 12, useDingbats = FALSE, height = 18)
ggplot(ptab, aes(x = t, y = neun_t)) +
    geom_hex(aes(fill=..density..), bins = 100) +
    # scale_x_log10() +
    # scale_y_log10() +
    facet_grid(type ~ de) +
    theme_bw(base_size = 30) +
    xlab('t-statistic: original') + 
    ylab('t-statistic: adj. cell fraction') +
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
#  date     2019-03-26
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
