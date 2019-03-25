library('sessioninfo')
library('SummarizedExperiment')
library('jaffelab')

## Load the original data
load('rda/rse_and_modQsva.Rdata', verbose = TRUE)
load('rda/expr_and_cleaned.Rdata', verbose = TRUE)

## Also load the original correlation on the cleaned data (keeping Dx)
load('rda/indv_corr.Rdata', verbose = TRUE)

## Re-clean the data to keep the Sex effect only
regions <- c('DLPFC', 'HIPPO')
mod_reorder <- function(mod) {
    i <- c(1, grep('Sex', colnames(mod)))
    j <- seq_len(ncol(mod))[-i]
    mod[, c(i, j)]
}
## Check that Sex is indeed the second column
lapply(modQsva, function(x) { colnames(mod_reorder(x))})
# $DLPFC
#  [1] "(Intercept)"       "SexM"              "DxSchizo"
#  [4] "Age"               "mitoRate"          "rRNA_rate"
#  [7] "totalAssignedGene" "RIN"               "snpPC1"
# [10] "snpPC2"            "snpPC3"            "snpPC4"
# [13] "snpPC5"            "PC1"               "PC2"
# [16] "PC3"               "PC4"               "PC5"
# [19] "PC6"               "PC7"               "PC8"
# [22] "PC9"               "PC10"              "PC11"
# [25] "PC12"              "PC13"              "PC14"
# [28] "PC15"
#
# $HIPPO
#  [1] "(Intercept)"       "SexM"              "DxSchizo"
#  [4] "Age"               "mitoRate"          "rRNA_rate"
#  [7] "totalAssignedGene" "RIN"               "snpPC1"
# [10] "snpPC2"            "snpPC3"            "snpPC4"
# [13] "snpPC5"            "PC1"               "PC2"
# [16] "PC3"               "PC4"               "PC5"
# [19] "PC6"               "PC7"               "PC8"
# [22] "PC9"               "PC10"              "PC11"
# [25] "PC12"              "PC13"              "PC14"
# [28] "PC15"              "PC16"

## keep only the fully annotated junctions
for(region in names(expr)) {
    index <- which(rowRanges(simple_rse[[region]]$jxn)$inGencode)
    simple_rse[[region]]$jxn <- simple_rse[[region]]$jxn[index, ]
    expr[[region]]$jxnRp10m <- expr[[region]]$jxnRp10m[index, ]
}
# > table(rowRanges(simple_rse[[region]]$jxn)$inGencode)
#
#  FALSE   TRUE
#  93582 203599
rm(region, index)


## Drop the chrX and chrY features
drop_xy <- function(exp) {
    
    mapply(function(elist, rselist) {
        mapply(function(e, rse) {
            not_xy <- !seqnames(rowRanges(rse)) %in% c('chrX', 'chrY')
            stopifnot(identical(rownames(e), names(rowRanges(rse))))
            message(paste(Sys.time(), 'dropping', sum(!not_xy), 'out of', length(not_xy), 'that is,', round(sum(!not_xy) / length(not_xy) * 100, 4), 'percent'))
            e[as.vector(not_xy), ]
        }, elist, rselist, SIMPLIFY = FALSE)
    }, exp, simple_rse, SIMPLIFY = FALSE)
    
}

expr <- drop_xy(expr)
# 2019-03-25 11:55:06 dropping 875 out of 24652 that is, 3.5494 percent
# 2019-03-25 11:55:06 dropping 11603 out of 396583 that is, 2.9257 percent
# 2019-03-25 11:55:07 dropping 6881 out of 203599 that is, 3.3797 percent
# 2019-03-25 11:55:07 dropping 2613 out of 92732 that is, 2.8178 percent
## HIPPO (same as DLPFC)
# 2019-03-25 11:55:07 dropping 875 out of 24652 that is, 3.5494 percent
# 2019-03-25 11:55:07 dropping 11603 out of 396583 that is, 2.9257 percent
# 2019-03-25 11:55:08 dropping 6881 out of 203599 that is, 3.3797 percent
# 2019-03-25 11:55:08 dropping 2613 out of 92732 that is, 2.8178 percent
## Not needed
# cleaned <- drop_xy(cleaned)


cleaned_sex <- lapply(regions, function(region) {
    lapply(expr[[region]], cleaningY, mod_reorder(modQsva[[region]]), P = 2)
})
names(cleaned_sex) <- regions

## Re-compute the correlation in the cleaned (keeping Sex) data
paircor <- function(x, y) {
    stopifnot(nrow(x) == nrow(y))
    stopifnot(ncol(x) == ncol(y))
    sapply(seq_len(ncol(x)), function(i) {
        cor(x[, i], y[, i])
    })
}

computecor <- function(exp) {
    sets <- names(exp[[1]])
    res <- lapply(sets, function(feature) {
        message(paste(Sys.time(), 'processing feature', feature))
        paircor(exp[['DLPFC']][[feature]], exp[['HIPPO']][[feature]])
    })
    names(res) <- sets
    return(res)
}
indv_cleaned_sex <- computecor(cleaned_sex)

## Save for later
save(cleaned_sex, file = 'rda/cleaned_sex.Rdata')
save(indv_cleaned_sex, file = 'rda/indv_cleaned_sex.Rdata')


pdf('pdf/indv_box_by_sex.pdf', useDingbats = FALSE)
mapply(function(x, y, set, type) {
    m <- !is.na(x)
    sex <- colData(y)$Sex[m]
    sex <- factor(ifelse(sex == 'M', 'Male', ifelse(sex == 'F', 'Female', 'hmm')))
    f <- lm(x[m] ~ sex)
    p <- summary(f)$coef[2, 4]
    ylim <- range(x[m])
    boxplot(x[m] ~ sex, main = paste(type, '-', set, '\n p-value:', signif(p, 3)),
        xlab = 'Sex', outline = FALSE, ylab = 'Correlation', ylim = ylim)
    points(x[m] ~ jitter(as.numeric(sex), amount = 0.15), cex = 1.5, pch = 21, bg = ifelse(sex == 'Male', '#E69F00', '#CC79A7'))
}, indv_cleaned_sex, simple_rse[[1]], names(indv_cleaned), 'cleaned expr (keeping Sex)')
dev.off()

## Plotting with Dx kept, yet Sex removed makes no sense
# mapply(function(x, y, set, type) {
#     m <- !is.na(x)
#     sex <- colData(y)$Sex[m]
#     sex <- factor(ifelse(sex == 'M', 'Male', ifelse(sex == 'F', 'Female', 'hmm')))
#     f <- lm(x[m] ~ sex)
#     p <- summary(f)$coef[2, 4]
#     ylim <- range(x[m])
#     boxplot(x[m] ~ sex, main = paste(type, '-', set, '\n p-value:', signif(p, 3)),
#         xlab = 'Sex', outline = FALSE, ylab = 'Correlation', ylim = ylim)
#     points(x[m] ~ jitter(as.numeric(sex), amount = 0.15), cex = 1.5, pch = 21, bg = ifelse(sex == 'Male', '#E69F00', '#CC79A7'))
# }, indv_cleaned, simple_rse[[1]], names(indv_cleaned), 'cleaned expr (keeping Dx)')



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
#  package              * version   date       lib source
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.6    2019-02-10 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                  0.8.0.1   2019-02-15 [1] CRAN (R 3.5.1)
#  GenomeInfoDb         * 1.18.2    2019-02-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2                3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab             * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  later                  0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma                  3.38.3    2018-12-02 [1] Bioconductor
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                  0.3.1     2019-03-03 [2] CRAN (R 3.5.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr                  0.13      2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.5       2019-02-20 [1] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library