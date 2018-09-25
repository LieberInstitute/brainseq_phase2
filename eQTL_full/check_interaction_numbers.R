library('data.table')
library('sessioninfo')

load('rdas/merged_GTEx_BrainSeq_QTLs_interaction.Rdata', verbose = TRUE)

lapply(interaction, function(x) {
    table(x$FDR < 0.01)
})
# $gene
#
#  TRUE
# 40092
#
# $exon
#
#  TRUE
# 89923
#
# $jxn
#
#  TRUE
# 75603
#
# $tx
#
#  TRUE
# 21895

sum(sapply(interaction[c('gene', 'exon', 'jxn')], function(x) {
    sum(x$FDR < 0.1)
}))
# [1] 205618

## Check that the data matches what we used previously in development/explore_limma_dev.R
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/matrixEqtl_output_interaction_4features.rda', verbose = TRUE)
me <- list('gene'= meGene$cis$eqtls, 'exon' = meExon$cis$eqtls, 'jxn' = meJxn$cis$eqtls, 'tx' = meTx$cis$eqtls)
sapply(me, function(x) table(x$FDR < 0.01) )
stopifnot(identical(sapply(me, function(x) sum(x$FDR < 0.01) ), sapply(interaction, nrow)))

load('rdas/merged_GTEx_BrainSeq_QTLs_hippo.Rdata', verbose = TRUE)
sapply(hippo, nrow)
stopifnot(sum(sapply(hippo[c('gene', 'exon', 'jxn')], nrow)) == 11237357)

load('rdas/merged_GTEx_BrainSeq_QTLs_dlpfc.Rdata', verbose = TRUE)
sapply(dlpfc, nrow)
stopifnot(sum(sapply(dlpfc[c('gene', 'exon', 'jxn')], nrow)) == 15766398)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.5.0 Patched (2018-04-30 r74679)
#  os       Red Hat Enterprise Linux Server release 6.9 (Santiago)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       US/Eastern
#  date     2018-09-25
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date       source
#  assertthat    0.2.0   2017-04-11 CRAN (R 3.5.0)
#  bindr         0.1.1   2018-03-13 CRAN (R 3.5.0)
#  bindrcpp      0.2.2   2018-03-29 CRAN (R 3.5.0)
#  clisymbols    1.2.0   2017-05-21 cran (@1.2.0)
#  colorout    * 1.2-0   2018-05-02 Github (jalvesaq/colorout@c42088d)
#  colorspace    1.3-2   2016-12-14 CRAN (R 3.5.0)
#  crayon        1.3.4   2017-09-16 CRAN (R 3.5.0)
#  data.table  * 1.11.6  2018-09-19 CRAN (R 3.5.0)
#  digest        0.6.17  2018-09-12 CRAN (R 3.5.0)
#  dplyr         0.7.6   2018-06-29 CRAN (R 3.5.0)
#  ggplot2       3.0.0   2018-07-03 CRAN (R 3.5.0)
#  glue          1.3.0   2018-07-17 CRAN (R 3.5.0)
#  gtable        0.2.0   2016-02-26 CRAN (R 3.5.0)
#  htmltools     0.3.6   2017-04-28 CRAN (R 3.5.0)
#  htmlwidgets   1.2     2018-04-19 CRAN (R 3.5.0)
#  httpuv        1.4.5   2018-07-19 CRAN (R 3.5.0)
#  later         0.7.5   2018-09-18 CRAN (R 3.5.0)
#  lattice       0.20-35 2017-03-25 CRAN (R 3.5.0)
#  lazyeval      0.2.1   2017-10-29 CRAN (R 3.5.0)
#  magrittr      1.5     2014-11-22 CRAN (R 3.5.0)
#  munsell       0.5.0   2018-06-12 CRAN (R 3.5.0)
#  pillar        1.3.0   2018-07-14 CRAN (R 3.5.0)
#  pkgconfig     2.0.2   2018-08-16 CRAN (R 3.5.0)
#  plyr          1.8.4   2016-06-08 CRAN (R 3.5.0)
#  png           0.1-7   2013-12-03 CRAN (R 3.5.0)
#  promises      1.0.1   2018-04-13 CRAN (R 3.5.0)
#  purrr         0.2.5   2018-05-29 CRAN (R 3.5.0)
#  R6            2.2.2   2017-06-17 CRAN (R 3.5.0)
#  Rcpp          0.12.18 2018-07-23 CRAN (R 3.5.0)
#  rlang         0.2.2   2018-08-16 CRAN (R 3.5.0)
#  rmote       * 0.3.4   2018-05-02 deltarho (R 3.5.0)
#  scales        1.0.0   2018-08-09 CRAN (R 3.5.0)
#  servr         0.10    2018-05-30 CRAN (R 3.5.0)
#  sessioninfo * 1.0.0   2017-06-21 CRAN (R 3.5.0)
#  tibble        1.4.2   2018-01-22 CRAN (R 3.5.0)
#  tidyselect    0.2.4   2018-02-26 CRAN (R 3.5.0)
#  withr         2.1.2   2018-03-15 CRAN (R 3.5.0)
#  xfun          0.3     2018-07-06 CRAN (R 3.5.0)
