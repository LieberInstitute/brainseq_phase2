library('recount')
library('devtools')

## Get the DLPFC GTEx samples
m <- all_metadata('gtex')
m2 <- subset(m, smafrze == 'USE ME' & smts == 'Brain')
table(m2$smtsd)
dlpfc <- subset(m2, smtsd == 'Brain - Frontal Cortex (BA9)')
# write.table(dlpfc$run, file = '~/Desktop/dlpfc_gtex_run_ids.txt', row.names = FALSE, quote = FALSE, col.names = FALSE)

## Find FASTQ files
p1 <- paste0('/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/FASTQ/', dlpfc$run, '_1.fastq.gz')
p2 <- paste0('/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/FASTQ/', dlpfc$run, '_2.fastq.gz')
stopifnot(all(file.exists(p1)))
stopifnot(all(file.exists(p2)))

## Write manifest files for the RNA-seq pipeline
manifest <- data.frame(
    p1,
    0,
    p2,
    0,
    dlpfc$run,
    stringsAsFactors = FALSE
)

write.table(manifest, file = 'samples.manifest', sep = '\t',
            row.names = FALSE, quote = FALSE, col.names = FALSE)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.4.3 Patched (2018-01-20 r74142)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       US/Eastern
#  date     2018-04-06
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  acepack                1.4.1     2016-10-29 CRAN (R 3.4.1)
#  AnnotationDbi          1.40.0    2017-11-29 Bioconductor
#  assertthat             0.2.0     2017-04-11 CRAN (R 3.4.1)
#  backports              1.1.2     2017-12-13 CRAN (R 3.4.2)
#  base                 * 3.4.3     2018-01-20 local
#  base64enc              0.1-3     2015-07-28 CRAN (R 3.4.1)
#  bindr                  0.1       2016-11-13 CRAN (R 3.4.1)
#  bindrcpp               0.2       2017-06-17 CRAN (R 3.4.1)
#  Biobase              * 2.38.0    2017-11-07 Bioconductor
#  BiocGenerics         * 0.24.0    2017-11-29 Bioconductor
#  BiocParallel           1.12.0    2017-11-29 Bioconductor
#  biomaRt                2.34.2    2018-02-17 Bioconductor
#  Biostrings             2.46.0    2017-11-29 Bioconductor
#  bit                    1.1-12    2014-04-09 CRAN (R 3.4.1)
#  bit64                  0.9-7     2017-05-08 CRAN (R 3.4.1)
#  bitops                 1.0-6     2013-08-17 CRAN (R 3.4.1)
#  blob                   1.1.0     2017-06-17 CRAN (R 3.4.1)
#  BSgenome               1.46.0    2017-11-29 Bioconductor
#  bumphunter             1.20.0    2017-11-29 Bioconductor
#  checkmate              1.8.5     2017-10-24 CRAN (R 3.4.2)
#  cluster                2.0.6     2017-03-10 CRAN (R 3.4.3)
#  codetools              0.2-15    2016-10-05 CRAN (R 3.4.3)
#  colorout             * 1.2-0     2018-02-19 Github (jalvesaq/colorout@2f01173)
#  colorspace             1.3-2     2016-12-14 CRAN (R 3.4.1)
#  compiler               3.4.3     2018-01-20 local
#  data.table             1.10.4-3  2017-10-27 CRAN (R 3.4.2)
#  datasets             * 3.4.3     2018-01-20 local
#  DBI                    0.8       2018-03-02 CRAN (R 3.4.3)
#  DelayedArray         * 0.4.1     2017-11-07 Bioconductor
#  derfinder              1.12.6    2018-02-17 Bioconductor
#  derfinderHelper        1.12.0    2017-11-29 Bioconductor
#  devtools             * 1.13.4    2017-11-09 CRAN (R 3.4.2)
#  digest                 0.6.15    2018-01-28 cran (@0.6.15)
#  doRNG                  1.6.6     2017-04-10 CRAN (R 3.4.1)
#  downloader             0.4       2015-07-09 CRAN (R 3.4.1)
#  dplyr                  0.7.4     2017-09-28 CRAN (R 3.4.1)
#  foreach                1.4.4     2017-12-12 CRAN (R 3.4.2)
#  foreign                0.8-69    2017-06-22 CRAN (R 3.4.3)
#  Formula                1.2-2     2017-07-10 CRAN (R 3.4.1)
#  GenomeInfoDb         * 1.14.0    2017-11-29 Bioconductor
#  GenomeInfoDbData       1.0.0     2018-01-09 Bioconductor
#  GenomicAlignments      1.14.1    2017-11-29 Bioconductor
#  GenomicFeatures        1.30.3    2018-02-17 Bioconductor
#  GenomicFiles           1.14.0    2017-11-29 Bioconductor
#  GenomicRanges        * 1.30.2    2018-02-17 Bioconductor
#  GEOquery               2.46.15   2018-03-06 Bioconductor
#  ggplot2                2.2.1     2016-12-30 CRAN (R 3.4.1)
#  glue                   1.2.0     2017-10-29 CRAN (R 3.4.2)
#  graphics             * 3.4.3     2018-01-20 local
#  grDevices            * 3.4.3     2018-01-20 local
#  grid                   3.4.3     2018-01-20 local
#  gridExtra              2.3       2017-09-09 CRAN (R 3.4.1)
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.4.1)
#  Hmisc                  4.1-1     2018-01-03 CRAN (R 3.4.2)
#  hms                    0.4.1     2018-01-24 CRAN (R 3.4.3)
#  htmlTable              1.11.2    2018-01-20 CRAN (R 3.4.3)
#  htmltools              0.3.6     2017-04-28 CRAN (R 3.4.1)
#  htmlwidgets            0.9       2017-07-10 CRAN (R 3.4.1)
#  httpuv                 1.3.6.2   2018-03-02 CRAN (R 3.4.3)
#  httr                   1.3.1     2017-08-20 CRAN (R 3.4.1)
#  IRanges              * 2.12.0    2017-11-29 Bioconductor
#  iterators              1.0.9     2017-12-12 CRAN (R 3.4.2)
#  jsonlite               1.5       2017-06-01 CRAN (R 3.4.1)
#  knitr                  1.20      2018-02-20 CRAN (R 3.4.3)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.4.3)
#  latticeExtra           0.6-28    2016-02-09 CRAN (R 3.4.1)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.4.2)
#  limma                  3.34.8    2018-02-17 Bioconductor
#  locfit                 1.5-9.1   2013-04-20 CRAN (R 3.4.1)
#  magrittr               1.5       2014-11-22 CRAN (R 3.4.1)
#  Matrix                 1.2-12    2017-11-30 CRAN (R 3.4.3)
#  matrixStats          * 0.53.1    2018-02-11 CRAN (R 3.4.3)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.4.1)
#  methods              * 3.4.3     2018-01-20 local
#  munsell                0.4.3     2016-02-13 CRAN (R 3.4.1)
#  nnet                   7.3-12    2016-02-02 CRAN (R 3.4.3)
#  parallel             * 3.4.3     2018-01-20 local
#  pillar                 1.1.0     2018-01-14 CRAN (R 3.4.2)
#  pkgconfig              2.0.1     2017-03-21 CRAN (R 3.4.1)
#  pkgmaker               0.22      2014-05-14 CRAN (R 3.4.1)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.4.1)
#  png                    0.1-7     2013-12-03 CRAN (R 3.4.1)
#  prettyunits            1.0.2     2015-07-13 CRAN (R 3.4.1)
#  progress               1.1.2     2016-12-14 CRAN (R 3.4.1)
#  purrr                  0.2.4     2017-10-18 CRAN (R 3.4.2)
#  qvalue                 2.10.0    2017-11-29 Bioconductor
#  R6                     2.2.2     2017-06-17 CRAN (R 3.4.1)
#  RColorBrewer           1.1-2     2014-12-07 CRAN (R 3.4.1)
#  Rcpp                   0.12.14   2017-11-23 CRAN (R 3.4.2)
#  RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.4.2)
#  readr                  1.1.1     2017-05-16 CRAN (R 3.4.1)
#  recount              * 1.4.4     2018-02-17 Bioconductor
#  registry               0.5       2017-12-03 CRAN (R 3.4.2)
#  rentrez                1.2.0     2018-02-12 CRAN (R 3.4.3)
#  reshape2               1.4.3     2017-12-11 CRAN (R 3.4.2)
#  rlang                  0.1.6     2017-12-21 CRAN (R 3.4.2)
#  rmote                * 0.3.4     2018-02-16 deltarho (R 3.4.3)
#  RMySQL                 0.10.14   2018-02-26 CRAN (R 3.4.3)
#  rngtools               1.2.4     2014-03-06 CRAN (R 3.4.1)
#  rpart                  4.1-12    2018-01-12 CRAN (R 3.4.3)
#  Rsamtools              1.30.0    2017-11-29 Bioconductor
#  RSQLite                2.0       2017-06-19 CRAN (R 3.4.1)
#  rstudioapi             0.7       2017-09-07 CRAN (R 3.4.1)
#  rtracklayer            1.38.3    2018-02-17 Bioconductor
#  S4Vectors            * 0.16.0    2017-11-29 Bioconductor
#  scales                 0.5.0     2017-08-24 CRAN (R 3.4.1)
#  servr                  0.8       2017-11-06 CRAN (R 3.4.3)
#  splines                3.4.3     2018-01-20 local
#  stats                * 3.4.3     2018-01-20 local
#  stats4               * 3.4.3     2018-01-20 local
#  stringi                1.1.6     2017-11-17 CRAN (R 3.4.2)
#  stringr                1.2.0     2017-02-18 CRAN (R 3.4.1)
#  SummarizedExperiment * 1.8.1     2018-01-09 Bioconductor
#  survival               2.41-3    2017-04-04 CRAN (R 3.4.3)
#  tibble                 1.4.1     2017-12-25 CRAN (R 3.4.2)
#  tidyr                  0.8.0     2018-01-29 CRAN (R 3.4.3)
#  tools                  3.4.3     2018-01-20 local
#  utils                * 3.4.3     2018-01-20 local
#  VariantAnnotation      1.24.5    2018-01-16 Bioconductor
#  withr                  2.1.1     2017-12-19 CRAN (R 3.4.2)
#  XML                    3.98-1.10 2018-02-19 CRAN (R 3.4.3)
#  xml2                   1.2.0     2018-01-24 CRAN (R 3.4.3)
#  xtable                 1.8-2     2016-02-05 CRAN (R 3.4.1)
#  XVector                0.18.0    2017-11-29 Bioconductor
#  zlibbioc               1.24.0    2017-11-07 Bioconductor
