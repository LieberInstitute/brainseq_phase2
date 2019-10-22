library('readxl')
library('SummarizedExperiment')
library('sgejobs')
library('sessioninfo')

## Ligases database from https://www.physiology.org/doi/full/10.1152/physiolgenomics.00031.2016
## Downloaded from https://hpcwebapps.cit.nih.gov/ESBL/Database/E3-ligases/
## wget https://hpcwebapps.cit.nih.gov/ESBL/Database/E3-ligases/Definite%20Ligase%20List.xlsx
ligases <- read_excel('Definite Ligase List.xlsx')
load('../expr_cutoff/rse_gene.Rdata', verbose = TRUE)
rse_gene_filt <- rse_gene
load('../expr_cutoff/unfiltered/rse_gene_unfiltered.Rdata', verbose = TRUE)
addmargins(table(
    pass_filter = tolower(ligases[[1]]) %in% tolower(rowRanges(rse_gene_filt)$Symbol),
    present = tolower(ligases[[1]]) %in% tolower(rowRanges(rse_gene)$Symbol)
))
#            present
# pass_filter TRUE Sum
#       FALSE   61  61
#       TRUE   316 316
#       Sum    377 377


ligases_found <- rowRanges(rse_gene_filt)$Symbol[
    tolower(rowRanges(rse_gene_filt)$Symbol) %in% tolower(ligases[[1]])
]
length(ligases_found)
# [1] 316

write.table(ligases_found, file = 'ligases_found.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

job_single('ligases_bsp2',
    create_shell = TRUE, email = 'a', task_num = length(ligases_found),
    command = "SYMBOL=$(awk '{ print $0}' ligases_found.txt | awk \"NR==${SGE_TASK_ID}\")\nRscript age_custom_plot.R -t gene -s ${SYMBOL} -o ligases_pdf")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.6.1 Patched (2019-09-06 r77160)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-10-22
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.5     2019-10-02 [1] CRAN (R 3.6.1)
#  Biobase              * 2.44.0    2019-05-02 [2] Bioconductor
#  BiocGenerics         * 0.30.0    2019-05-02 [1] Bioconductor
#  BiocParallel         * 1.18.1    2019-08-06 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.6.1)
#  cellranger             1.1.0     2016-07-27 [1] CRAN (R 3.6.1)
#  cli                    1.1.0     2019-03-19 [1] CRAN (R 3.6.1)
#  codetools              0.2-16    2018-12-24 [3] CRAN (R 3.6.1)
#  colorout             * 1.2-2     2019-09-26 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.6.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.1)
#  DelayedArray         * 0.10.0    2019-05-02 [2] Bioconductor
#  digest                 0.6.21    2019-09-20 [1] CRAN (R 3.6.1)
#  dplyr                  0.8.3     2019-07-04 [1] CRAN (R 3.6.1)
#  fansi                  0.4.0     2018-10-05 [1] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.20.0    2019-05-02 [1] Bioconductor
#  GenomeInfoDbData       1.2.1     2019-09-09 [2] Bioconductor
#  GenomicRanges        * 1.36.1    2019-09-06 [1] Bioconductor
#  ggplot2                3.2.1     2019-08-10 [1] CRAN (R 3.6.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.6.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.6.1)
#  hms                    0.5.1     2019-08-23 [2] CRAN (R 3.6.1)
#  htmltools              0.4.0     2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5       2019-10-04 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2     2019-09-11 [1] CRAN (R 3.6.1)
#  IRanges              * 2.18.3    2019-09-24 [1] Bioconductor
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.6.1)
#  later                  1.0.0     2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.6.1)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.6.1)
#  lifecycle              0.1.0     2019-08-01 [1] CRAN (R 3.6.1)
#  lubridate              1.7.4     2018-04-11 [1] CRAN (R 3.6.1)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.6.1)
#  pillar                 1.4.2     2019-06-29 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.6.1)
#  promises               1.1.0     2019-10-04 [1] CRAN (R 3.6.1)
#  pryr                   0.1.4     2018-02-18 [2] CRAN (R 3.6.1)
#  purrr                  0.3.2     2019-03-15 [2] CRAN (R 3.6.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.2     2019-07-25 [1] CRAN (R 3.6.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.6.1)
#  readr                  1.3.1     2018-12-21 [1] CRAN (R 3.6.1)
#  readxl               * 1.3.1     2019-03-13 [2] CRAN (R 3.6.1)
#  rlang                  0.4.0     2019-06-25 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4     2019-09-26 [1] Github (cloudyr/rmote@fbce611)
#  S4Vectors            * 0.22.1    2019-09-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.6.1)
#  servr                  0.15      2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.6.1)
#  sgejobs              * 0.99.0    2019-09-26 [1] Github (LieberInstitute/sgejobs@94a6151)
#  stringi                1.4.3     2019-03-12 [2] CRAN (R 3.6.1)
#  stringr                1.4.0     2019-02-10 [1] CRAN (R 3.6.1)
#  SummarizedExperiment * 1.14.1    2019-07-31 [1] Bioconductor
#  tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.1)
#  tidyr                  1.0.0     2019-09-11 [1] CRAN (R 3.6.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.6.1)
#  utf8                   1.1.4     2018-05-24 [1] CRAN (R 3.6.1)
#  vctrs                  0.2.0     2019-07-05 [1] CRAN (R 3.6.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.6.1)
#  xfun                   0.10      2019-10-01 [1] CRAN (R 3.6.1)
#  XVector                0.24.0    2019-05-02 [1] Bioconductor
#  zeallot                0.1.0     2018-01-28 [1] CRAN (R 3.6.1)
#  zlibbioc               1.30.0    2019-05-02 [2] Bioconductor
#
# [1] /users/lcollado/R/3.6
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6/R/3.6/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6/R/3.6/lib64/R/library
#
