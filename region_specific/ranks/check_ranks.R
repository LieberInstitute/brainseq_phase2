library('sessioninfo')
library('purrr')
library('dplyr')
library('ggplot2')
library('ffpe')

## Load model results including BrainSpan
message(paste(Sys.time(), 'loading ../rda/pcheck_both.Rdata'))
load('../rda/pcheck_both.Rdata', verbose = TRUE)

## Rename here so it'll be easier later
pcheck_both$age[pcheck_both$age == 'fetal'] <- 'prenatal'
pcheck_both$span_age[pcheck_both$span_age == 'fetal'] <- 'prenatal'

## Save the feature id, otherwise we lose it with purrr::map
pcheck_both$feature_id <- gsub('.*[gene|exon|jxn|tx]\\.', '', rownames(pcheck_both))

## Split by age group, then feature
features <- c('gene', 'exon', 'jxn', 'tx')
info <- map(
    split(pcheck_both, factor(pcheck_both$age, levels = c('prenatal', 'adult'))),
    ~ split(.x, factor(.x$type, levels = features))
)


## Now add ranks
info <- map(info, ~ map(.x, ~ mutate(.data = .x, rank = rank(P.Value), rank_span = rank(span_P.Value)) ))
# x <- info[[1]][[1]]
# > length(unique(x$rank_span))
# [1] 23250
# > length(unique(x$rank))
# [1] 24637
# > nrow(x)
# [1] 24652

## Make it faster by just looking at the top 5k
maxrank <- 5000
info <- map(info, ~ map(.x, function(x) {
    message(paste(Sys.time(), 'processing', x$age[1], x$type[1]))
    x <- x[order(x$P.Value), ]
    id_sorted <- x$feature_id
    id_sorted_span <- x$feature_id[order(x$span_P.Value)]
    x$i <- seq_len(nrow(x))
    ## The concordance here is calculated like in http://leekgroup.github.io/recount-analyses/example_de/recount_SRP019936.html
    x$concordance <- map_int(x$i, ~ if(.x > maxrank) NA else sum(id_sorted[seq_len(.x)] %in% id_sorted_span[seq_len(.x)]))
    return(x)
}))

info_df <- map_dfr(info, ~ map_dfr(.x, ~ .x))
info_df$type <- factor(info_df$type, levels = features)

pdf('concordance_by_age.pdf', useDingbats = FALSE, width = 25)
p <- ggplot(subset(info_df, !is.na(concordance)), aes(x = i, y = concordance, color = age)) +
    geom_line(size = 2) +
    facet_grid(~ type) +
    geom_abline(slope = 1, intercept = 0, color = 'black', linetype = 'dashed', size = 1.5) +
    xlab('Ordered features in BSP2') +
    ylab('Ordered features in BrainSpan') +
    theme_bw(base_size = 30) +
    scale_color_manual(values = c('prenatal' = 'orange', 'adult' = 'purple'))
p
p + scale_x_log10() + scale_y_log10()
dev.off()


## Use ffpe::CATplot() instead
info_cat <- map2_dfr(
    info,
    names(info),
    function(.x, agegrp) {
        map2_dfr(
            .x,
            names(.x),
            ~ mutate(CATplot(.x$rank, .x$rank_span, maxrank = maxrank, ylim = c(0, 1), make.plot = FALSE), age = agegrp, feature = .y)
        )
    }
)

dat_ablines <- data.frame(feature = features, intercept = 0, slope = 1 / map_int(info[[1]], nrow))
dat_ablines
#      feature intercept        slope
# gene    gene         0 4.056466e-05
# exon    exon         0 2.521566e-06
# jxn      jxn         0 3.364953e-06
# tx        tx         0 1.078376e-05

pdf('concordance_vs_chance_by_age.pdf', useDingbats = FALSE, width = 25)
ggplot(info_cat, aes(x = rank, y = concordance, color = age)) +
    geom_line(size = 2) +
    facet_grid(~ factor(feature, levels = features)) +
    geom_abline(aes(slope = slope, intercept = intercept), data = dat_ablines, color = 'black', linetype = 'dashed', size = 1.5) +
    xlab('Rank by p-value') +
    ylab('Concordance') +
    theme_bw(base_size = 30) +
    scale_color_manual(values = c('prenatal' = 'orange', 'adult' = 'purple')) +
    labs(caption = 'Dotted line denotes expection by chance')
ggplot(subset(info_cat, rank <= 1000), aes(x = rank, y = concordance, color = age)) +
    geom_line(size = 2) +
    facet_grid(~ factor(feature, levels = features)) +
    geom_abline(aes(slope = slope, intercept = intercept), data = dat_ablines, color = 'black', linetype = 'dashed', size = 1.5) +
    xlab('Rank by p-value') +
    ylab('Concordance') +
    theme_bw(base_size = 30) +
    scale_color_manual(values = c('prenatal' = 'orange', 'adult' = 'purple')) +
    labs(caption = 'Dotted line denotes expection by chance')
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
#  date     2019-03-13
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  affy                   1.60.0    2018-10-30 [2] Bioconductor
#  affyio                 1.52.0    2018-10-30 [2] Bioconductor
#  annotate               1.60.0    2018-10-30 [1] Bioconductor
#  AnnotationDbi          1.44.0    2018-10-30 [1] Bioconductor
#  askpass                1.1       2019-01-13 [1] CRAN (R 3.5.1)
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  base64                 2.0       2016-05-10 [2] CRAN (R 3.5.0)
#  beanplot               1.2       2014-09-19 [2] CRAN (R 3.5.0)
#  bibtex                 0.4.2     2017-06-30 [1] CRAN (R 3.5.0)
#  bindr                  0.1.1     2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp             * 0.2.2     2018-03-29 [1] CRAN (R 3.5.0)
#  Biobase                2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics           0.28.0    2018-10-30 [1] Bioconductor
#  BiocManager            1.30.4    2018-11-13 [1] CRAN (R 3.5.1)
#  BiocParallel           1.16.5    2019-01-04 [1] Bioconductor
#  biomaRt                2.38.0    2018-10-30 [1] Bioconductor
#  Biostrings             2.50.2    2019-01-03 [1] Bioconductor
#  bit                    1.1-14    2018-05-29 [2] CRAN (R 3.5.0)
#  bit64                  0.9-7     2017-05-08 [2] CRAN (R 3.5.0)
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  blob                   1.1.1     2018-03-25 [2] CRAN (R 3.5.0)
#  bumphunter             1.24.5    2018-12-01 [1] Bioconductor
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  codetools              0.2-15    2016-10-05 [3] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  curl                   3.3       2019-01-10 [1] CRAN (R 3.5.1)
#  data.table             1.12.0    2019-01-13 [1] CRAN (R 3.5.1)
#  DBI                    1.0.0     2018-05-02 [2] CRAN (R 3.5.0)
#  DelayedArray           0.8.0     2018-10-30 [2] Bioconductor
#  DelayedMatrixStats     1.4.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  doRNG                  1.7.1     2018-06-22 [2] CRAN (R 3.5.1)
#  dplyr                * 0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
#  ffpe                 * 1.26.0    2018-10-30 [1] Bioconductor
#  foreach                1.4.4     2017-12-12 [2] CRAN (R 3.5.0)
#  genefilter             1.64.0    2018-10-30 [1] Bioconductor
#  GenomeInfoDb           1.18.1    2018-11-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicAlignments      1.18.1    2019-01-04 [1] Bioconductor
#  GenomicFeatures        1.34.1    2018-11-03 [1] Bioconductor
#  GenomicRanges          1.34.0    2018-10-30 [1] Bioconductor
#  GEOquery               2.50.5    2018-12-22 [1] Bioconductor
#  ggplot2              * 3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.0     2018-07-17 [1] CRAN (R 3.5.1)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  HDF5Array              1.10.1    2018-12-05 [1] Bioconductor
#  hms                    0.4.2     2018-03-10 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  httr                   1.4.0     2018-12-11 [1] CRAN (R 3.5.1)
#  illuminaio             0.24.0    2018-10-30 [2] Bioconductor
#  IRanges                2.16.0    2018-10-30 [1] Bioconductor
#  iterators              1.0.10    2018-07-13 [2] CRAN (R 3.5.1)
#  KernSmooth             2.23-15   2015-06-29 [3] CRAN (R 3.5.1)
#  labeling               0.3       2014-08-23 [2] CRAN (R 3.5.0)
#  later                  0.7.5     2018-09-18 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma                  3.38.3    2018-12-02 [1] Bioconductor
#  locfit                 1.5-9.1   2013-04-20 [2] CRAN (R 3.5.0)
#  lumi                   2.34.0    2018-10-30 [1] Bioconductor
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  MASS                   7.3-51.1  2018-11-01 [3] CRAN (R 3.5.1)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats            0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  mclust                 5.4.2     2018-11-17 [1] CRAN (R 3.5.1)
#  memoise                1.1.0     2017-04-21 [2] CRAN (R 3.5.0)
#  methylumi              2.28.0    2018-10-30 [1] Bioconductor
#  mgcv                   1.8-26    2018-11-21 [3] CRAN (R 3.5.1)
#  mime                   0.6       2018-10-05 [1] CRAN (R 3.5.1)
#  minfi                  1.28.3    2019-01-05 [2] Bioconductor
#  multtest               2.38.0    2018-10-30 [2] Bioconductor
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.0)
#  nleqslv                3.3.2     2018-05-17 [2] CRAN (R 3.5.0)
#  nlme                   3.1-137   2018-04-07 [3] CRAN (R 3.5.1)
#  nor1mix                1.2-3     2017-08-30 [2] CRAN (R 3.5.0)
#  openssl                1.2.1     2019-01-17 [1] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  pkgmaker               0.27      2018-05-25 [2] CRAN (R 3.5.0)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  preprocessCore         1.44.0    2018-10-30 [2] Bioconductor
#  prettyunits            1.0.2     2015-07-13 [1] CRAN (R 3.5.0)
#  progress               1.2.0     2018-06-14 [1] CRAN (R 3.5.1)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                * 0.2.5     2018-05-29 [2] CRAN (R 3.5.0)
#  quadprog               1.5-5     2013-04-17 [2] CRAN (R 3.5.0)
#  R6                     2.3.0     2018-10-04 [2] CRAN (R 3.5.1)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  readr                  1.3.1     2018-12-21 [1] CRAN (R 3.5.1)
#  registry               0.5       2017-12-03 [2] CRAN (R 3.5.0)
#  reshape                0.8.8     2018-10-23 [2] CRAN (R 3.5.1)
#  reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.5.0)
#  rhdf5                  2.26.2    2019-01-02 [2] Bioconductor
#  Rhdf5lib               1.4.2     2018-12-03 [1] Bioconductor
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  rngtools               1.3.1     2018-05-15 [2] CRAN (R 3.5.0)
#  Rsamtools              1.34.0    2018-10-30 [1] Bioconductor
#  RSQLite                2.1.1     2018-05-06 [2] CRAN (R 3.5.0)
#  rtracklayer            1.42.1    2018-11-22 [1] Bioconductor
#  S4Vectors              0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  servr                  0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  sfsmisc                1.1-3     2018-11-27 [1] CRAN (R 3.5.1)
#  siggenes               1.56.0    2018-10-30 [2] Bioconductor
#  stringi                1.2.4     2018-07-20 [2] CRAN (R 3.5.1)
#  stringr                1.3.1     2018-05-10 [1] CRAN (R 3.5.0)
#  SummarizedExperiment   1.12.0    2018-10-30 [1] Bioconductor
#  survival               2.43-3    2018-11-26 [3] CRAN (R 3.5.1)
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyr                  0.8.2     2018-10-28 [2] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  TTR                  * 0.23-4    2018-09-20 [1] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  XML                    3.98-1.16 2018-08-19 [2] CRAN (R 3.5.1)
#  xml2                   1.2.0     2018-01-24 [2] CRAN (R 3.5.0)
#  xtable                 1.8-3     2018-08-29 [2] CRAN (R 3.5.1)
#  xts                    0.11-2    2018-11-05 [1] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#  zoo                    1.8-4     2018-09-19 [2] CRAN (R 3.5.1)
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
