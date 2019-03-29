library('clusterProfiler')
library('sessioninfo')

go_files <- c(
    'sczd' = '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/go_de_genes.Rdata',
    'high_corr' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/correlation/rda/go_corr_genes.Rdata',
    'dev' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/development/rda/go_de_genes_novenn.Rdata',
    'reg_adult' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/region_specific/rda/go_de_genes_brain.Rdata',
    'reg_prenatal' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/region_specific/rda/go_de_genes_brain.Rdata'
)
stopifnot(all(file.exists(go_files)))

go_obj <- c(
    'sczd' = 'go_de_genes',
    'high_corr' = 'go_corr_genes',
    'dev' = 'go_de_genes_novenn',
    'reg_adult' = 'go_de_genes_brain_adult',
    'reg_prenatal' = 'go_de_genes_brain_fetal'
)

load_go <- function(file, obj) {
    load(file, verbose = TRUE)
    get(obj)
}

extract_go_table <- function(go, set) {
    ## Keep only full results
    go <- go[sapply(go, class) == 'compareClusterResult']
    
    res <- do.call(rbind, mapply(function(x, type) { cbind(x@compareClusterResult, type = type) }, go, names(go), SIMPLIFY = FALSE))
    res$set <- set
    rownames(res) <- NULL
    return(res)
}

go_list <- mapply(load_go, go_files, go_obj, SIMPLIFY = FALSE)


set_names <- c(
    'sczd' = 'SCZD case-control differences',
    'high_corr' = 'DLPFC and HIPPO high correlation',
    'dev' = 'DLPFC and HIPPO differential development',
    'reg_adult' = 'DLPFC and HIPPO adult differences',
    'reg_prenatal' = 'DLPFC and HIPPO prenatal differences'
)

go_table <- do.call(rbind, mapply(extract_go_table, go_list, set_names[names(go_files)], SIMPLIFY = FALSE))

dim(go_table)
# [1] 2774   12
head(go_table)

write.table(go_table, file = 'TableSxx_go_table.txt', sep = '\t', quote = FALSE, row.names = FALSE)
save(go_table, file = 'go_table.Rdata')


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
#  date     2019-03-29
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package         * version  date       lib source
#  AnnotationDbi     1.44.0   2018-10-30 [1] Bioconductor
#  assertthat        0.2.1    2019-03-21 [2] CRAN (R 3.5.1)
#  Biobase           2.42.0   2018-10-30 [2] Bioconductor
#  BiocGenerics      0.28.0   2018-10-30 [1] Bioconductor
#  BiocParallel      1.16.6   2019-02-10 [1] Bioconductor
#  bit               1.1-14   2018-05-29 [2] CRAN (R 3.5.1)
#  bit64             0.9-7    2017-05-08 [2] CRAN (R 3.5.0)
#  blob              1.1.1    2018-03-25 [2] CRAN (R 3.5.0)
#  cli               1.0.1    2018-09-25 [1] CRAN (R 3.5.1)
#  clusterProfiler * 3.10.1   2018-12-20 [1] Bioconductor
#  colorout        * 1.2-0    2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace        1.4-1    2019-03-18 [2] CRAN (R 3.5.1)
#  cowplot           0.9.4    2019-01-08 [1] CRAN (R 3.5.1)
#  crayon            1.3.4    2017-09-16 [1] CRAN (R 3.5.0)
#  data.table        1.12.0   2019-01-13 [1] CRAN (R 3.5.1)
#  DBI               1.0.0    2018-05-02 [2] CRAN (R 3.5.0)
#  digest            0.6.18   2018-10-10 [1] CRAN (R 3.5.1)
#  DO.db             2.9      2018-05-03 [1] Bioconductor
#  DOSE              3.8.2    2019-01-14 [1] Bioconductor
#  dplyr             0.8.0.1  2019-02-15 [1] CRAN (R 3.5.1)
#  enrichplot        1.2.0    2018-10-30 [2] Bioconductor
#  europepmc         0.3      2018-04-20 [2] CRAN (R 3.5.1)
#  farver            1.1.0    2018-11-20 [1] CRAN (R 3.5.1)
#  fastmatch         1.1-0    2017-01-28 [1] CRAN (R 3.5.0)
#  fgsea             1.8.0    2018-10-30 [1] Bioconductor
#  ggforce           0.2.1    2019-03-12 [2] CRAN (R 3.5.1)
#  ggplot2           3.1.0    2018-10-25 [1] CRAN (R 3.5.1)
#  ggplotify         0.0.3    2018-08-03 [2] CRAN (R 3.5.1)
#  ggraph            1.0.2    2018-07-07 [2] CRAN (R 3.5.1)
#  ggrepel           0.8.0    2018-05-09 [1] CRAN (R 3.5.0)
#  ggridges          0.5.1    2018-09-27 [1] CRAN (R 3.5.1)
#  glue              1.3.1    2019-03-12 [1] CRAN (R 3.5.1)
#  GO.db             3.7.0    2018-11-02 [1] Bioconductor
#  GOSemSim          2.8.0    2018-10-30 [1] Bioconductor
#  gridExtra         2.3      2017-09-09 [2] CRAN (R 3.5.0)
#  gridGraphics      0.3-0    2018-04-03 [2] CRAN (R 3.5.1)
#  gtable            0.3.0    2019-03-25 [2] CRAN (R 3.5.1)
#  hms               0.4.2    2018-03-10 [2] CRAN (R 3.5.0)
#  htmltools         0.3.6    2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets       1.3      2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv            1.5.0    2019-03-15 [2] CRAN (R 3.5.1)
#  httr              1.4.0    2018-12-11 [1] CRAN (R 3.5.1)
#  igraph            1.2.4    2019-02-13 [2] CRAN (R 3.5.1)
#  IRanges           2.16.0   2018-10-30 [1] Bioconductor
#  jsonlite          1.6      2018-12-07 [2] CRAN (R 3.5.1)
#  later             0.8.0    2019-02-11 [2] CRAN (R 3.5.1)
#  lattice           0.20-38  2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval          0.2.2    2019-03-15 [2] CRAN (R 3.5.1)
#  magrittr          1.5      2014-11-22 [1] CRAN (R 3.5.0)
#  MASS              7.3-51.1 2018-11-01 [3] CRAN (R 3.5.1)
#  Matrix            1.2-17   2019-03-22 [3] CRAN (R 3.5.1)
#  memoise           1.1.0    2017-04-21 [2] CRAN (R 3.5.0)
#  munsell           0.5.0    2018-06-12 [2] CRAN (R 3.5.1)
#  pillar            1.3.1    2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig         2.0.2    2018-08-16 [1] CRAN (R 3.5.1)
#  plyr              1.8.4    2016-06-08 [2] CRAN (R 3.5.0)
#  png               0.1-7    2013-12-03 [2] CRAN (R 3.5.0)
#  polyclip          1.10-0   2019-03-14 [2] CRAN (R 3.5.1)
#  prettyunits       1.0.2    2015-07-13 [1] CRAN (R 3.5.0)
#  progress          1.2.0    2018-06-14 [1] CRAN (R 3.5.1)
#  promises          1.0.1    2018-04-13 [2] CRAN (R 3.5.0)
#  purrr             0.3.2    2019-03-15 [2] CRAN (R 3.5.1)
#  qvalue            2.14.1   2019-01-10 [1] Bioconductor
#  R6                2.4.0    2019-02-14 [2] CRAN (R 3.5.1)
#  RColorBrewer      1.1-2    2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp              1.0.0    2018-11-07 [1] CRAN (R 3.5.1)
#  reshape2          1.4.3    2017-12-11 [2] CRAN (R 3.5.0)
#  rlang             0.3.1    2019-01-08 [1] CRAN (R 3.5.1)
#  rmote           * 0.3.4    2018-05-02 [1] deltarho (R 3.5.0)
#  RSQLite           2.1.1    2018-05-06 [2] CRAN (R 3.5.0)
#  rvcheck           0.1.3    2018-12-06 [1] CRAN (R 3.5.1)
#  S4Vectors         0.20.1   2018-11-09 [1] Bioconductor
#  scales            1.0.0    2018-08-09 [2] CRAN (R 3.5.1)
#  servr             0.13     2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo     * 1.1.1    2018-11-05 [1] CRAN (R 3.5.1)
#  stringi           1.4.3    2019-03-12 [2] CRAN (R 3.5.1)
#  stringr           1.4.0    2019-02-10 [1] CRAN (R 3.5.1)
#  tibble            2.0.1    2019-01-12 [1] CRAN (R 3.5.1)
#  tidyr             0.8.3    2019-03-01 [2] CRAN (R 3.5.1)
#  tidyselect        0.2.5    2018-10-11 [2] CRAN (R 3.5.1)
#  triebeard         0.3.0    2016-08-04 [1] CRAN (R 3.5.0)
#  tweenr            1.0.1    2018-12-14 [1] CRAN (R 3.5.1)
#  UpSetR            1.3.3    2017-03-21 [1] CRAN (R 3.5.0)
#  urltools          1.7.2    2019-02-04 [1] CRAN (R 3.5.1)
#  viridis           0.5.1    2018-03-29 [2] CRAN (R 3.5.0)
#  viridisLite       0.3.0    2018-02-01 [2] CRAN (R 3.5.0)
#  withr             2.1.2    2018-03-15 [2] CRAN (R 3.5.0)
#  xfun              0.5      2019-02-20 [1] CRAN (R 3.5.1)
#  xml2              1.2.0    2018-01-24 [2] CRAN (R 3.5.0)
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
