library('SummarizedExperiment')
library('gplots')
library('sessioninfo')
library('purrr')
library('dplyr')
library('ggplot2')
library('RColorBrewer')

## Load unfiltered RSEs
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_gene_unfiltered.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_exon_unfiltered.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_jxn_unfiltered.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_tx_unfiltered.Rdata", verbose = TRUE)

rses <- list(
    'gene' = rse_gene,
    'exon' = rse_exon,
    'jxn' = rse_jxn,
    'tx' = rse_tx
)

## From https://github.com/LieberInstitute/brainseq_phase2/blob/master/development/explore_limma_dev.R#L702-L718
## Only for those passing the expression cutoffs (matches the list from the explore_limma_dev.R script)
gene_ens <- lapply(names(rses), function(feat) {
    if(feat %in% c('gene', 'exon', 'jxn')) {
        res <- rowRanges(rses[[feat]])$ensemblID[rowRanges(rses[[feat]])$passExprsCut]
    } else {
        res <- gsub('\\..*', '', rowRanges(rses[[feat]])$gene_id[rowRanges(rses[[feat]])$passExprsCut])
    }
    res <- unique(res[!is.na(res)])
    return(res)
})
names(gene_ens) <- names(rses)


## Group by ensembl gene id
venn_info <- venn(gene_ens, show.plot = FALSE)

## Find those only present in genes
gene_only <- attr(venn_info, 'intersections')$gene
stopifnot(length(gene_only) == 955)


## Get all the ensembl IDs
gene_ens_all <- map2(rses, names(rses), 
    ~ if(.y != 'tx') rowRanges(.x)$ensemblID else gsub('\\..*', '', rowRanges(.x)$gene_id)
)
map_int(gene_ens_all, length)
# gene   exon    jxn     tx
# 58037 571623 837077 198093

## For each of the 955 genes, find the features that match them
i_map <- map(gene_ens_all[-1], function(x) {
    i <- which(x %in% gene_only)
    m <- match(x, gene_only)
    m2 <- m[!is.na(m)]
    stopifnot(identical(length(i), length(m2)))
    res <- split(i, gene_only[m2])
    
    ## Add missing ones
    missing <- gene_only[!gene_only %in% names(res)]
    if(length(missing) > 0) {
        res2 <- vector('list', length(missing))
        names(res2) <- missing
        res <- c(res, res2)
    }
    
    ## Re-order
    res[gene_only]
})

## Find the total number of missing genes
map_int(i_map, ~ sum(map_int(.x, length) == 0))
# exon  jxn   tx
#    1  868    0
## Next, in percent
map_int(i_map, ~ sum(map_int(.x, length) == 0)) / length(gene_only) * 100
#     exon       jxn        tx
# 0.104712 90.890052  0.000000

## Get the mean expression for each feature from these genes
mean_expr <- map2(
    i_map,
    names(i_map),
    ~ data.frame(
        meanExprs = rowRanges(rses[[.y]])$meanExprs[ unlist(.x)],
        gene_id = rep(names(.x), map_int(.x, length)),
        passExprsCut = rowRanges(rses[[.y]])$passExprsCut[ unlist(.x)],
        feature = .y,
        i = unlist(.x),
        stringsAsFactors = FALSE
    )
)

## Also combine into a single table
mean_expr_tab <- map_dfr(mean_expr, function(x) { x$feature_id <- rownames(x); return(x) })

## Should all be 0
stopifnot(all(map_int(mean_expr, ~ sum(.x$passExprsCut)) == 0))


## Get the highest expressed feature for each gene
max_expr <- map_dfr(mean_expr, function(x) {
    y <- split(x, x$gene_id)  
    i_max <- map_int(y, ~ .x$i[which.max(.x$meanExprs)])
    i_m <- match(i_max, x$i)
    res <- x[i_m, ]
    res$n_feat <- map_int(y, length)
    res$feature_id <- rownames(res)
    return(res)
})
## Check numbers
table(max_expr$feature)
# exon  jxn   tx
#  954   87  955

length(gene_only) - map_int(i_map, ~ sum(map_int(.x, length) == 0))
# exon  jxn   tx
#  954   87  955

## Add back the 0s
missing <- map(i_map[c('exon', 'jxn')], ~ which(map_int(.x, length) == 0))
missing <- map2_dfr(missing, names(missing),
    ~ data.frame(
        meanExprs = 0,
        gene_id = names(.x),
        passExprsCut = FALSE,
        feature = .y,
        i = NA,
        n_feat = 0,
        feature_id = NA,
        stringsAsFactors = FALSE
    )
)
max_expr_all <- rbind(max_expr, missing)

## Get the original feature colors we've been using
colors <- brewer.pal('Set1', n = 4)[2:4]
names(colors) <- names(mean_expr)

## Add cutoff lines used
dat_hlines <- data.frame(feature=c('exon', 'jxn', 'tx'), hline=c(0.3, 0.46, 0.38))
dat_hlines

## Visualize the maximum expr
pdf('max_expr_feature_by_gene.pdf', useDingbats = FALSE)
ggplot(max_expr_all, aes(y = meanExprs, x = feature, fill = feature)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_point(aes(fill = feature), shape = 21, position = position_jitter(width = 0.2)) +
    theme_bw(base_size = 30) +
    guides(fill = FALSE) +
    scale_fill_manual(values = colors) +
    geom_errorbar(data=dat_hlines, aes(y=NULL, ymax=hline, ymin=hline), colour="#AA0000", size = 1.5) +
    labs(y = 'Max Mean Expr\nfor each gene', x = 'Feature', caption = 'Red lines denote expression cutoff')
dev.off()
    
## Check the gene type of these 955 genes
gtype <- with(rowRanges(rses$gene), sort(table(gene_type[ensemblID %in% gene_only]), decreasing = TRUE))
cbind(gtype, percent = gtype / length(gene_only) * 100)
#                                    gtype    percent
# processed_pseudogene                 415 43.4554974
# lincRNA                              160 16.7539267
# antisense                            110 11.5183246
# protein_coding                        68  7.1204188
# sense_intronic                        43  4.5026178
# unprocessed_pseudogene                37  3.8743455
# TEC                                   36  3.7696335
# miRNA                                 31  3.2460733
# transcribed_processed_pseudogene      22  2.3036649
# sense_overlapping                      8  0.8376963
# misc_RNA                               7  0.7329843
# snRNA                                  4  0.4188482
# processed_transcript                   3  0.3141361
# rRNA                                   3  0.3141361
# snoRNA                                 2  0.2094241
# IG_C_pseudogene                        1  0.1047120
# IG_V_pseudogene                        1  0.1047120
# non_coding                             1  0.1047120
# transcribed_unitary_pseudogene         1  0.1047120
# transcribed_unprocessed_pseudogene     1  0.1047120
# unitary_pseudogene                     1  0.1047120

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
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  bindr                  0.1.1     2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 [1] CRAN (R 3.5.0)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.5    2019-01-04 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  caTools                1.17.1.1  2018-07-20 [2] CRAN (R 3.5.1)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                * 0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
#  gdata                  2.18.0    2017-06-06 [2] CRAN (R 3.5.0)
#  GenomeInfoDb         * 1.18.1    2018-11-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2              * 3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.0     2018-07-17 [1] CRAN (R 3.5.1)
#  gplots               * 3.0.1     2016-03-30 [1] CRAN (R 3.5.0)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  gtools                 3.8.1     2018-06-26 [2] CRAN (R 3.5.1)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  KernSmooth             2.23-15   2015-06-29 [3] CRAN (R 3.5.1)
#  labeling               0.3       2014-08-23 [2] CRAN (R 3.5.0)
#  later                  0.7.5     2018-09-18 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.0)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                * 0.2.5     2018-05-29 [2] CRAN (R 3.5.0)
#  R6                     2.3.0     2018-10-04 [2] CRAN (R 3.5.1)
#  RColorBrewer         * 1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.5.0)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  servr                  0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  stringi                1.2.4     2018-07-20 [2] CRAN (R 3.5.1)
#  stringr                1.3.1     2018-05-10 [1] CRAN (R 3.5.0)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
