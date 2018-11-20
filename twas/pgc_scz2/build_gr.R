library('GenomicRanges')
library('data.table')
library('sessioninfo')

## Check https://www.med.unc.edu/pgc/files/resultfiles/scz2.readme.pdf
## for file descriptions
scz2 <- fread('ckqny.scz2snpres')
scz2_gr_hg19 <- GRanges(
    seqnames = scz2$hg19chr,
    IRanges(start = scz2$bp, width = 1)
)
mcols(scz2_gr_hg19) <- scz2[, - c('hg19chrc', 'bp')]
save(scz2_gr_hg19, file = 'scz2_gr_hg19.Rdata')

## I'm using the "start of LD association interval (defined by r2>0.6)"
## saved under 'six1' and the end (saved under 'six2')
## to define the regions instead of 'bp' which
## is the "hg19 base pair position of SNP with minimum p-value"
scz2_rep <- fread('scz2.rep.128.txt')
scz2_rep_gr_hg19 <- GRanges(
    seqnames = scz2_rep$hg19chrc,
    IRanges(start = scz2_rep$six1, end = scz2_rep$six2)
)
mcols(scz2_rep_gr_hg19) <- scz2_rep[, - c('hg19chrc', 'six1', 'six2', 'chr')]
save(scz2_rep_gr_hg19, file = 'scz2_rep_gr_hg19.Rdata')

## Read in the "anneal" file
## then use the snp file to build a GRanges
scz2_anneal <- fread('scz2.anneal.108.txt')
length(unique(scz2_anneal$bestsnp))
m <- match(unique(scz2_anneal$bestsnp), scz2_gr_hg19$snpid)
stopifnot(length(m) == nrow(scz2_anneal))
scz2_anneal_gr_hg19 <- scz2_gr_hg19[m]
mcols(scz2_anneal_gr_hg19) <- cbind(mcols(scz2_anneal_gr_hg19), scz2_anneal[, -c('hg19chrc')])
save(scz2_anneal_gr_hg19, file = 'scz2_anneal_gr_hg19.Rdata')

## Save them all together just in case
save(scz2_gr_hg19, scz2_rep_gr_hg19, scz2_anneal_gr_hg19, file = 'scz2_all.Rdata')

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
#  date     2018-11-20
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package          * version   date       source
#  assertthat         0.2.0     2017-04-11 CRAN (R 3.5.0)
#  bindr              0.1.1     2018-03-13 CRAN (R 3.5.0)
#  bindrcpp           0.2.2     2018-03-29 CRAN (R 3.5.0)
#  BiocGenerics     * 0.26.0    2018-05-03 Bioconductor
#  bitops             1.0-6     2013-08-17 CRAN (R 3.5.0)
#  clisymbols         1.2.0     2017-05-21 cran (@1.2.0)
#  colorout         * 1.2-0     2018-05-02 Github (jalvesaq/colorout@c42088d)
#  colorspace         1.3-2     2016-12-14 CRAN (R 3.5.0)
#  crayon             1.3.4     2017-09-16 CRAN (R 3.5.0)
#  data.table       * 1.11.6    2018-09-19 CRAN (R 3.5.0)
#  digest             0.6.17    2018-09-12 CRAN (R 3.5.0)
#  dplyr              0.7.6     2018-06-29 CRAN (R 3.5.0)
#  GenomeInfoDb     * 1.16.0    2018-05-03 Bioconductor
#  GenomeInfoDbData   1.1.0     2018-04-17 Bioconductor
#  GenomicRanges    * 1.32.7    2018-09-20 Bioconductor
#  ggplot2            3.0.0     2018-07-03 CRAN (R 3.5.0)
#  glue               1.3.0     2018-07-17 CRAN (R 3.5.0)
#  gtable             0.2.0     2016-02-26 CRAN (R 3.5.0)
#  htmltools          0.3.6     2017-04-28 CRAN (R 3.5.0)
#  htmlwidgets        1.2       2018-04-19 CRAN (R 3.5.0)
#  httpuv             1.4.5     2018-07-19 CRAN (R 3.5.0)
#  IRanges          * 2.14.12   2018-09-20 Bioconductor
#  later              0.7.5     2018-09-18 CRAN (R 3.5.0)
#  lattice            0.20-35   2017-03-25 CRAN (R 3.5.0)
#  lazyeval           0.2.1     2017-10-29 CRAN (R 3.5.0)
#  magrittr           1.5       2014-11-22 CRAN (R 3.5.0)
#  munsell            0.5.0     2018-06-12 CRAN (R 3.5.0)
#  pillar             1.3.0     2018-07-14 CRAN (R 3.5.0)
#  pkgconfig          2.0.2     2018-08-16 CRAN (R 3.5.0)
#  plyr               1.8.4     2016-06-08 CRAN (R 3.5.0)
#  png                0.1-7     2013-12-03 CRAN (R 3.5.0)
#  promises           1.0.1     2018-04-13 CRAN (R 3.5.0)
#  purrr              0.2.5     2018-05-29 CRAN (R 3.5.0)
#  R6                 2.2.2     2017-06-17 CRAN (R 3.5.0)
#  Rcpp               0.12.18   2018-07-23 CRAN (R 3.5.0)
#  RCurl              1.95-4.11 2018-07-15 CRAN (R 3.5.0)
#  rlang              0.2.2     2018-08-16 CRAN (R 3.5.0)
#  rmote            * 0.3.4     2018-05-02 deltarho (R 3.5.0)
#  S4Vectors        * 0.18.3    2018-06-13 Bioconductor
#  scales             1.0.0     2018-08-09 CRAN (R 3.5.0)
#  servr              0.10      2018-05-30 CRAN (R 3.5.0)
#  sessioninfo      * 1.0.0     2017-06-21 CRAN (R 3.5.0)
#  tibble             1.4.2     2018-01-22 CRAN (R 3.5.0)
#  tidyselect         0.2.4     2018-02-26 CRAN (R 3.5.0)
#  withr              2.1.2     2018-03-15 CRAN (R 3.5.0)
#  xfun               0.3       2018-07-06 CRAN (R 3.5.0)
#  XVector            0.20.0    2018-05-03 Bioconductor
#  zlibbioc           1.26.0    2018-05-02 Bioconductor
