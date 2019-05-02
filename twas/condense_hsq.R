library('SummarizedExperiment')
library('sessioninfo')

files <- c(
    'DLPFC' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/DLPFC/gene/heritability/hsq_info.Rdata',
    'HIPPO' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene/heritability/hsq_info.Rdata'
)

# rse_files <- c(
#     'DLPFC' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/DLPFC/gene/subsetted_rse.Rdata',
#     'HIPPO' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene/subsetted_rse.Rdata'
# )
#
# rse <- lapply(rse_files, function(f) {
#     load(f, verbose = TRUE)
#     return(rse)
# })

## Genes are the same Same for both regions
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/DLPFC/gene/subsetted_rse.Rdata', verbose = TRUE)
genes <- names(rowRanges(rse))

hsq <- lapply(files, function(f) {
    load(f, verbose = TRUE)
    
    missing <- genes[!genes %in% hsq_info$ID]
    if(length(missing) > 0) {
        add <- data.frame(
            HSQ = NA,
            ID = missing,
            CHR = gsub('chr', '', as.character(seqnames(rowRanges(rse)[missing]))),
            P0 = start(rowRanges(rse)[missing]),
            P1 = end(rowRanges(rse)[missing]),
            geneID = missing,
            hsq = NA,
            hsq.se = NA,
            hsq.pv = NA,
            region = hsq_info$region[1],
            feature = hsq_info$feature[1],
            failFilter = NA,
            stringsAsFactors = FALSE
        )
        hsq_info <- rbind(hsq_info, add)
    }
    
    hsq_info <- hsq_info[match(genes, hsq_info$ID), ]
    rownames(hsq_info) <- NULL
    
    return(hsq_info)
})
names(hsq) <- names(files)

## Genes that didn't converge
sapply(hsq, function(x) { sum(is.na(x$HSQ))})
# DLPFC HIPPO
#   568   600

## Genes that didn't converge across regions
addmargins(table(
    'DLPFC' = is.na(hsq[['DLPFC']]$HSQ),
    'HIPPO' = is.na(hsq[['HIPPO']]$HSQ)
))
#        HIPPO
# DLPFC   FALSE  TRUE   Sum
#   FALSE 22284   550 22834
#   TRUE    518    50   568
#   Sum   22802   600 23402

## Number of genes that pass the heritability filters
sapply(hsq, function(x) { table(!x$failFilter, useNA = 'ifany')})
#       DLPFC HIPPO
# FALSE 15147 16934
# TRUE   7687  5868
# <NA>    568   600

## genes passing filters across regions
addmargins(table(
    'DLPFC' = !hsq[['DLPFC']]$failFilter,
    'HIPPO' = !hsq[['HIPPO']]$failFilter,
    useNA = 'ifany'
))
#        HIPPO
# DLPFC   FALSE  TRUE  <NA>   Sum
#   FALSE 13622  1020   505 15147
#   TRUE   2806  4836    45  7687
#   <NA>    506    12    50   568
#   Sum   16934  5868   600 23402

## Explore a bit how the heritability across regions looks
plot(x = hsq[['DLPFC']]$hsq, y = hsq[['HIPPO']]$hsq, xlab = 'DLPFC', ylab = 'HIPPO')

save(hsq, file = 'rda/hsq_genes.Rdata')


## Resume
load('rda/hsq_genes.Rdata', verbose = TRUE)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
                     function(x)
                       rgb(x[1], x[2], x[3], alpha=alpha))
}

pdf('pdf/heritability_DLPFC_HIPPO_gene.pdf', useDingbats = FALSE)
plot(
    x = (hsq[['DLPFC']]$hsq + hsq[['HIPPO']]$hsq) / 2,
    y = hsq[['DLPFC']]$hsq - hsq[['HIPPO']]$hsq,
    xlab = '(DLPFC + HIPPO) / 2',
    ylab = 'DLPFC - HIPPO',
    pch = 16,
    col = ifelse(
        !(hsq[['DLPFC']]$failFilter | hsq[['HIPPO']]$failFilter),
        add.alpha('magenta', 1/5),
        ifelse(
            hsq[['DLPFC']]$hsq < 0 | hsq[['HIPPO']]$hsq < 0,
            add.alpha('royalblue4', 1/5),
            add.alpha('black', 1/5)    
        )
    ),
    main = 'Heritability (cis +- 500kb)',
    cex.lab = 1.4,
    cex.axis = 1.4,
    cex.main = 1.5
)
abline(h = 0, col = 'red', lwd = 1.2)
legend('bottomright',
    legend = c(
        '<0 in 1 region',
        'Fails filter in 1 region',
        'Passes filters in both regions'
    ),
    col = c('royalblue4', 'black', 'magenta'), lwd = 2
)

plot(
    x = hsq[['HIPPO']]$hsq,
    y = hsq[['DLPFC']]$hsq,
    xlab = 'HIPPO',
    ylab = 'DLPFC',
    pch = 16,
    col = ifelse(
        !(hsq[['DLPFC']]$failFilter | hsq[['HIPPO']]$failFilter),
        add.alpha('magenta', 1/5),
        ifelse(
            hsq[['DLPFC']]$hsq < 0 | hsq[['HIPPO']]$hsq < 0,
            add.alpha('royalblue4', 1/5),
            add.alpha('black', 1/5)    
        )
    ),
    main = 'Heritability (cis +- 500kb)',
    cex.lab = 1.4,
    cex.axis = 1.4,
    cex.main = 1.5
)
abline(a = 0, b = 1, col = 'red', lwd = 1.2)
legend('bottomright',
    legend = c(
        '<0 in 1 region',
        'Fails filter in 1 region',
        'Passes filters in both regions'
    ),
    col = c('royalblue4', 'black', 'magenta'), lwd = 2
)
dev.off()


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.5.3 Patched (2019-03-11 r76311)
#  os       Red Hat Enterprise Linux Server release 6.9 (Santiago)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-05-01
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.5.1)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.6    2019-02-10 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  cli                    1.1.0     2019-03-19 [1] CRAN (R 3.5.3)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                  0.8.0.1   2019-02-15 [1] CRAN (R 3.5.1)
#  GenomeInfoDb         * 1.18.2    2019-02-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2                3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.5.1)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.5.0     2019-03-15 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  later                  0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.3)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.5.1)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.3)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                  0.3.2     2019-03-15 [2] CRAN (R 3.5.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  Rcpp                   1.0.1     2019-03-17 [1] CRAN (R 3.5.3)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  rlang                  0.3.3     2019-03-29 [1] CRAN (R 3.5.3)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  servr                  0.13      2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 2.1.1     2019-03-16 [1] CRAN (R 3.5.3)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.6       2019-04-02 [1] CRAN (R 3.5.3)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
