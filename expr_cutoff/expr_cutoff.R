library('SummarizedExperiment')
library('recount')
library('jaffelab')
library('devtools')

## Load the data
rse <- lapply(c('../count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseGene_merged_n449.rda', '../count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n442.rda'), function(f) {
    load(f)
    rowRanges(rse_gene)$meanExprs <- NA
    return(rse_gene)
})

## Calculate RPKMs
rpkm <- lapply(rse, recount::getRPKM, 'Length')
expr <- do.call(cbind, rpkm)

## Identify potential cutoffs
pdf('suggested_expr_cutoffs.pdf', width = 12)
expression_cutoff(expr, seed = 20171026)
# 2017-10-26 16:35:47 the suggested expression cutoff is 0.21
# percent_features_cut  samples_nonzero_cut
#                 0.25                 0.17
dev.off()


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.4.2 (2017-09-28)
#  system   x86_64, darwin15.6.0
#  ui       AQUA
#  language (EN)
#  collate  en_US.UTF-8
#  tz       America/New_York
#  date     2017-10-26
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version  date       source
#  acepack                1.4.1    2016-10-29 CRAN (R 3.4.0)
#  AnnotationDbi          1.39.4   2017-10-12 Bioconductor
#  assertthat             0.2.0    2017-04-11 cran (@0.2.0)
#  backports              1.1.1    2017-09-25 CRAN (R 3.4.2)
#  base                 * 3.4.2    2017-10-23 local
#  base64enc              0.1-3    2015-07-28 cran (@0.1-3)
#  bindr                  0.1      2016-11-13 CRAN (R 3.4.0)
#  bindrcpp               0.2      2017-06-17 CRAN (R 3.4.0)
#  Biobase              * 2.37.2   2017-05-05 Bioconductor
#  BiocGenerics         * 0.23.4   2017-10-21 Bioconductor
#  BiocParallel           1.11.13  2017-10-22 Bioconductor
#  biomaRt                2.33.4   2017-08-02 Bioconductor
#  Biostrings             2.45.4   2017-08-28 Bioconductor
#  bit                    1.1-12   2014-04-09 CRAN (R 3.4.0)
#  bit64                  0.9-7    2017-05-08 CRAN (R 3.4.0)
#  bitops                 1.0-6    2013-08-17 cran (@1.0-6)
#  blob                   1.1.0    2017-06-17 CRAN (R 3.4.0)
#  BSgenome               1.45.3   2017-09-15 Bioconductor
#  bumphunter             1.18.0   2017-10-05 Bioconductor
#  checkmate              1.8.4    2017-09-25 CRAN (R 3.4.2)
#  cluster                2.0.6    2017-03-10 CRAN (R 3.4.2)
#  codetools              0.2-15   2016-10-05 CRAN (R 3.4.2)
#  colorspace             1.3-2    2016-12-14 cran (@1.3-2)
#  compiler               3.4.2    2017-10-23 local
#  data.table             1.10.4-2 2017-10-12 CRAN (R 3.4.2)
#  datasets             * 3.4.2    2017-10-23 local
#  DBI                    0.7      2017-06-18 CRAN (R 3.4.0)
#  DelayedArray         * 0.3.21   2017-09-27 Bioconductor
#  derfinder              1.11.10  2017-10-12 Bioconductor
#  derfinderHelper        1.11.4   2017-10-12 Bioconductor
#  devtools             * 1.13.3   2017-08-02 CRAN (R 3.4.1)
#  digest                 0.6.12   2017-01-27 CRAN (R 3.4.0)
#  doRNG                  1.6.6    2017-04-10 CRAN (R 3.4.0)
#  downloader             0.4      2015-07-09 CRAN (R 3.4.0)
#  dplyr                  0.7.4    2017-09-28 CRAN (R 3.4.2)
#  foreach                1.4.3    2015-10-13 CRAN (R 3.4.0)
#  foreign                0.8-69   2017-06-22 CRAN (R 3.4.2)
#  Formula                1.2-2    2017-07-10 CRAN (R 3.4.1)
#  GenomeInfoDb         * 1.13.5   2017-10-05 Bioconductor
#  GenomeInfoDbData       0.99.1   2017-07-17 Bioconductor
#  GenomicAlignments      1.13.6   2017-09-15 Bioconductor
#  GenomicFeatures        1.29.13  2017-10-13 Bioconductor
#  GenomicFiles           1.13.10  2017-07-17 Bioconductor
#  GenomicRanges        * 1.29.15  2017-10-05 Bioconductor
#  GEOquery               2.45.2   2017-09-26 Bioconductor
#  ggplot2                2.2.1    2016-12-30 CRAN (R 3.4.0)
#  glue                   1.1.1    2017-06-21 CRAN (R 3.4.1)
#  graphics             * 3.4.2    2017-10-23 local
#  grDevices            * 3.4.2    2017-10-23 local
#  grid                   3.4.2    2017-10-23 local
#  gridExtra              2.3      2017-09-09 CRAN (R 3.4.1)
#  gtable                 0.2.0    2016-02-26 CRAN (R 3.4.0)
#  Hmisc                  4.0-3    2017-05-02 CRAN (R 3.4.0)
#  hms                    0.3      2016-11-22 cran (@0.3)
#  htmlTable              1.9      2017-01-26 CRAN (R 3.4.0)
#  htmltools              0.3.6    2017-04-28 CRAN (R 3.4.0)
#  htmlwidgets            0.9      2017-07-10 CRAN (R 3.4.1)
#  httr                   1.3.1    2017-08-20 CRAN (R 3.4.1)
#  IRanges              * 2.11.19  2017-10-08 Bioconductor
#  iterators              1.0.8    2015-10-13 CRAN (R 3.4.0)
#  jaffelab             * 0.99.15  2017-10-26 Github (LieberInstitute/jaffelab@94307b0)
#  jsonlite               1.5      2017-06-01 CRAN (R 3.4.0)
#  knitr                  1.17     2017-08-10 CRAN (R 3.4.1)
#  lattice                0.20-35  2017-03-25 CRAN (R 3.4.2)
#  latticeExtra           0.6-28   2016-02-09 CRAN (R 3.4.0)
#  lazyeval               0.2.0    2016-06-12 cran (@0.2.0)
#  limma                  3.33.14  2017-10-12 Bioconductor
#  locfit                 1.5-9.1  2013-04-20 CRAN (R 3.4.0)
#  magrittr               1.5      2014-11-22 cran (@1.5)
#  Matrix                 1.2-11   2017-08-21 CRAN (R 3.4.2)
#  matrixStats          * 0.52.2   2017-04-14 CRAN (R 3.4.0)
#  memoise                1.1.0    2017-04-21 CRAN (R 3.4.0)
#  methods              * 3.4.2    2017-10-23 local
#  munsell                0.4.3    2016-02-13 cran (@0.4.3)
#  nnet                   7.3-12   2016-02-02 CRAN (R 3.4.2)
#  parallel             * 3.4.2    2017-10-23 local
#  pkgconfig              2.0.1    2017-03-21 CRAN (R 3.4.0)
#  pkgmaker               0.22     2014-05-14 CRAN (R 3.4.0)
#  plyr                   1.8.4    2016-06-08 cran (@1.8.4)
#  prettyunits            1.0.2    2015-07-13 CRAN (R 3.4.0)
#  progress               1.1.2    2016-12-14 CRAN (R 3.4.0)
#  purrr                  0.2.4    2017-10-18 CRAN (R 3.4.2)
#  qvalue                 2.9.0    2017-05-04 Bioconductor
#  R6                     2.2.2    2017-06-17 CRAN (R 3.4.0)
#  rafalib              * 1.0.0    2015-08-09 cran (@1.0.0)
#  RColorBrewer           1.1-2    2014-12-07 cran (@1.1-2)
#  Rcpp                   0.12.13  2017-09-28 CRAN (R 3.4.1)
#  RCurl                  1.95-4.8 2016-03-01 cran (@1.95-4.)
#  readr                  1.1.1    2017-05-16 CRAN (R 3.4.0)
#  recount              * 1.3.12   2017-10-13 Bioconductor
#  registry               0.3      2015-07-08 CRAN (R 3.4.0)
#  rentrez                1.1.0    2017-06-01 CRAN (R 3.4.0)
#  reshape2               1.4.2    2016-10-22 cran (@1.4.2)
#  rlang                  0.1.2    2017-08-09 CRAN (R 3.4.1)
#  rngtools               1.2.4    2014-03-06 CRAN (R 3.4.0)
#  rpart                  4.1-11   2017-03-13 CRAN (R 3.4.2)
#  Rsamtools              1.29.1   2017-08-19 Bioconductor
#  RSQLite                2.0      2017-06-19 CRAN (R 3.4.1)
#  rtracklayer            1.37.3   2017-07-22 Bioconductor
#  S4Vectors            * 0.15.14  2017-10-16 Bioconductor
#  scales                 0.5.0    2017-08-24 CRAN (R 3.4.1)
#  segmented              0.5-2.2  2017-09-19 CRAN (R 3.4.2)
#  splines                3.4.2    2017-10-23 local
#  stats                * 3.4.2    2017-10-23 local
#  stats4               * 3.4.2    2017-10-23 local
#  stringi                1.1.5    2017-04-07 cran (@1.1.5)
#  stringr                1.2.0    2017-02-18 cran (@1.2.0)
#  SummarizedExperiment * 1.7.10   2017-09-29 Bioconductor
#  survival               2.41-3   2017-04-04 CRAN (R 3.4.2)
#  tibble                 1.3.4    2017-08-22 CRAN (R 3.4.1)
#  tidyr                  0.7.2    2017-10-16 CRAN (R 3.4.2)
#  tools                  3.4.2    2017-10-23 local
#  utils                * 3.4.2    2017-10-23 local
#  VariantAnnotation      1.23.9   2017-10-12 Bioconductor
#  withr                  2.0.0    2017-07-28 CRAN (R 3.4.1)
#  XML                    3.98-1.9 2017-06-19 CRAN (R 3.4.1)
#  xml2                   1.1.1    2017-01-24 cran (@1.1.1)
#  xtable                 1.8-2    2016-02-05 cran (@1.8-2)
#  XVector                0.17.2   2017-10-18 Bioconductor
#  zlibbioc               1.23.0   2017-05-04 Bioconductor
