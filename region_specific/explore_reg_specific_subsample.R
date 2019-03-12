library('sessioninfo')
library('jaffelab')
library('purrr')
library('ggplot2')


dir.create('rda', showWarnings = FALSE)
dir.create('pdf', showWarnings = FALSE)

## Load BrainSeq model subsampled results
f_sub <- list.files('subsample/rda', pattern = 'limma_region_specific_adult_gene_[0-9]*.Rdata', full.names = TRUE)

top_list <- map(f_sub, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f)
    top$age <- 'adult'
    top$type <- 'gene'
    top$P.Bonf <- p.adjust(top$P.Value, 'bonf')
    return(top)
})
names(top_list) <- gsub('.*gene_|.Rdata', '', f_sub)

## Some of the 100 subsampled runs failed
## like number 63 that failed with:
# Error in if (abs(correlation) >= 1) stop("correlation is 1 or -1, so the model is degenerate") :
#   missing value where TRUE/FALSE needed
# Calls: system.time -> lmFit -> gls.series
length(top_list)
# [1] 82

message(paste(Sys.time(), 'loading rda/de_genes.Rdata'))
load('rda/de_genes.Rdata', verbose = TRUE)

message(paste(Sys.time(), 'loading rda/pcheck_both.Rdata'))
load('rda/pcheck_both.Rdata', verbose = TRUE)

## Make it easier to work with
adult <- subset(pcheck_both, age == 'adult' & type == 'gene')
fetal <- subset(pcheck_both, age == 'fetal' & type == 'gene')

## Match the gene names
genes <- gsub('adult_gene.', '', rownames(adult))
stopifnot(identical(rownames(top_list[[1]]), rownames(top_list[[2]])))
m <- match(genes, rownames(top_list[[1]]))

## Extract the number of times a gene was P bonf < 1%
de_mat <- map_dfc(top_list, ~ .x$P.Bonf[m] < 0.01)
adult$sub_freq <- rowSums(de_mat)
adult$sub_prop <- adult$sub_freq / length(top_list)
de_mat_span <- map_dfc(top_list, ~ .x$P.Bonf[m] < 0.01 & adult$span_P.Value < 0.05 & sign(.x$t[m]) == sign(adult$span_t))
adult$sub_freq_span <- rowSums(de_mat_span)
adult$sub_prop_span <- adult$sub_freq_span / length(top_list)

get_de <- function(x) {
    sign(x$t) == sign(x$span_t) & x$span_P.Value < 0.05 & x$P.Bonf < 0.01
}
adult$de <- get_de(adult)
stopifnot(sum(adult$de) == 1612)
fetal$de <- get_de(fetal)
stopifnot(sum(fetal$de) == 32)

tab_pbonf <- with(adult, table('DE original' = de, 'P Bonf < 1% in at least one rep' = sub_freq > 0))
addmargins(tab_pbonf)
#            P Bonf < 1% in at least one rep
# DE original FALSE  TRUE   Sum
#       FALSE 20808  2232 23040
#       TRUE    497  1115  1612
#       Sum   21305  3347 24652
getOR(tab_pbonf)
# [1] 20.91484
chisq.test(tab_pbonf)
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  tab_pbonf
# X-squared = 4537.7, df = 1, p-value < 2.2e-16

tab_span <- with(adult, table('DE original' = de, 'P Bonf <1% rep in BrainSpan > 1 sub rep' = sub_freq_span > 0))
addmargins(tab_span)
#            P Bonf <1% rep in BrainSpan > 1 sub rep
# DE original FALSE  TRUE   Sum
#       FALSE 23037     3 23040
#       TRUE    497  1115  1612
#       Sum   23534  1118 24652
getOR(tab_span)
# [1] 17227.54
chisq.test(tab_span)
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  tab_span
# X-squared = 16627, df = 1, p-value < 2.2e-16




tab_pbonf_only <- with(adult, table('Original P.Bonf < 1%' = P.Bonf < 0.01, 'P Bonf < 1% in at least one rep' = sub_freq > 0))
addmargins(tab_pbonf_only)
#                     P Bonf < 1% in at least one rep
# Original P.Bonf < 1% FALSE  TRUE   Sum
#                FALSE 16718    39 16757
#                TRUE   4587  3308  7895
#                Sum   21305  3347 24652
getOR(tab_pbonf_only)
# [1] 309.1409
chisq.test(tab_pbonf_only)
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  tab_pbonf_only
# X-squared = 7937, df = 1, p-value < 2.2e-16

# > chisq.test(tab_span)$p.value
# [1] 0
# > chisq.test(tab_pbonf)$p.value
# [1] 0
# > chisq.test(tab_pbonf_only)$p.value
# [1] 0

## Number of genes with P Bonf < 1% in the original
sum(adult$P.Bonf < 0.01)
# [1] 7895
sum(fetal$P.Bonf < 0.01)
# [1] 70
## And then the distribution of that number across the 82 sub-sampled runs
n_de <- colSums(de_mat)
summary(n_de)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  4.0   123.5   324.5   353.9   409.5  1642.0

t.test(n_de, mu = 70)
#     One Sample t-test
#
# data:  n_de
# t = 8.5926, df = 81, p-value = 5.08e-13
# alternative hypothesis: true mean is not equal to 70
# 95 percent confidence interval:
#  288.1811 419.6725
# sample estimates:
# mean of x
#  353.9268

## Number of DE genes (P Bonf < 1% that replicate in BrainSpan)
sum(adult$de)
# [1] 1612
sum(fetal$de)
# [1] 32
## And then their distribution across the 82 replicates
n_de_span <- colSums(de_mat_span)
summary(n_de_span)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  3.0   101.2   217.5   217.9   281.8   727.0

t.test(n_de_span, mu = 32)
#     One Sample t-test
#
# data:  n_de_span
# t = 11.257, df = 81, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 32
# 95 percent confidence interval:
#  185.0343 250.7462
# sample estimates:
# mean of x
#  217.8902
t.test(n_de_span, mu = 32)$p.value
# [1] 3.06133e-18


pdf('subsample/pdf/number_de_genes_adult_subsampled_to_prenatal_numbers.pdf', useDingbats = FALSE)
ggplot(data.frame(var = n_de, group = ''), aes(y = var, x = group, fill = 'red')) + 
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_point(shape = 21, position = position_jitter(width = 0.2)) +
    # geom_point() +
    theme_bw(base_size = 30) + 
    ylab('# genes with P Bonf<1%') + 
    xlab('') +
    scale_fill_manual(values = c('red' = '#CC79A7')) +
    guides(fill=FALSE) +
    labs(caption = 'Does not use BrainSpan') +
    geom_hline(aes(yintercept = 70), color = '#E69F00', linetype = 'dashed', size = 1.5)

ggplot(data.frame(var = n_de_span, group = ''), aes(y = var, x = group, fill = 'red')) + 
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_point(shape = 21, position = position_jitter(width = 0.2)) +
    # geom_point() +
    theme_bw(base_size = 30) + 
    ylab('# genes with P Bonf<1%') + 
    xlab('') +
    scale_fill_manual(values = c('red' = '#0072B2')) +
    guides(fill=FALSE) +
    labs(caption = 'Requires BrainSpan replication:\nSame sign and P span < 5%')+
    geom_hline(aes(yintercept = 32), color = '#E69F00', linetype = 'dashed', size = 1.5)
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
#  date     2019-03-12
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package          * version   date       lib source
#  assertthat         0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  bindr              0.1.1     2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp           0.2.2     2018-03-29 [1] CRAN (R 3.5.0)
#  BiocGenerics       0.28.0    2018-10-30 [1] Bioconductor
#  bitops             1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  cli                1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout         * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace         1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon             1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  digest             0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr              0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
#  GenomeInfoDb       1.18.1    2018-11-12 [1] Bioconductor
#  GenomeInfoDbData   1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges      1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2          * 3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue               1.3.0     2018-07-17 [1] CRAN (R 3.5.1)
#  gtable             0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools          0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets        1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv             1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges            2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab         * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  later              0.7.5     2018-09-18 [2] CRAN (R 3.5.1)
#  lattice            0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval           0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma              3.38.3    2018-12-02 [1] Bioconductor
#  magrittr           1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  munsell            0.5.0     2018-06-12 [2] CRAN (R 3.5.0)
#  pillar             1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig          2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr               1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises           1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr            * 0.2.5     2018-05-29 [2] CRAN (R 3.5.0)
#  R6                 2.3.0     2018-10-04 [2] CRAN (R 3.5.1)
#  rafalib          * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer       1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp               1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl              1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  rlang              0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote            * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors          0.20.1    2018-11-09 [1] Bioconductor
#  scales             1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented          0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr              0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo      * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  tibble             2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect         0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr              2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun               0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  XVector            0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc           1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library

