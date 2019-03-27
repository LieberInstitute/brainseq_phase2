library('sessioninfo')
library('purrr')
library('jaffelab')
library('ggplot2')

outFeat <- lapply(
    c('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda'), function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    outTx$ensemblID <- gsub('\\..*', '', outTx$gene_id)
    return(list('gene' = outGene, 'exon' = outExon, 'jxn' = outJxn, 'tx' = outTx))
})
names(outFeat) <- c('DLPFC', 'HIPPO')


outFeat_cell <- lapply(
    c('/dcl01/lieber/ajaffe/lab/brainseq_phase2/cellComp/check_sczd/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda',
    '/dcl01/lieber/ajaffe/lab/brainseq_phase2/cellComp/check_sczd/rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda'), function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    outTx$ensemblID <- gsub('\\..*', '', outTx$gene_id)
    return(list('gene' = outGene, 'exon' = outExon, 'jxn' = outJxn, 'tx' = outTx))
})
names(outFeat_cell) <- c('DLPFC', 'HIPPO')

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
                     function(x)
                       rgb(x[1], x[2], x[3], alpha=alpha))
}

comp_log <- function(x, y, xlab, ylab, var = 'logFC', de = FALSE, n = 150, onlyx = FALSE) {
    if(de) {
        if(!onlyx) {
            common <- unique(c(
                head(x$ensemblID[order(x$adj.P.Val, decreasing = FALSE)], n),
                head(y$ensemblID[order(y$adj.P.Val, decreasing = FALSE)], n)
            ))
        } else {
            common <- unique(c(
                head(x$ensemblID[order(x$adj.P.Val, decreasing = FALSE)], n)
            ))
        }
        
    } else {
        common <- intersect(x$ensemblID, y$ensemblID)
    }
    x <- x[match(common, x$ensemblID), ]
    y <- y[match(common, y$ensemblID), ]
    corr = signif(cor(x[, var], y[, var], use = 'pairwise.complete.obs'), 3)
    
    ## Use colors by FDR < 0.05
    colgroup <- rep('none', nrow(x))
    colgroup[x$adj.P.Val < 0.05] <- 'x'
    colgroup[y$adj.P.Val < 0.05] <- 'y'
    colgroup[x$adj.P.Val < 0.05 & y$adj.P.Val < 0.05] <- 'xy'
    cols <- dplyr::case_when(
        colgroup == 'none' ~ add.alpha('black', ifelse(de, 1/2, 1/10)),
        colgroup == 'x' ~ add.alpha('magenta', ifelse(de, 1/2, 1/4)),
        colgroup == 'y' ~ add.alpha('magenta', ifelse(de, 1/2, 1/4)),
        colgroup == 'xy' ~ add.alpha('royalblue4', 1/2)
    )

    par(cex.axis = 2, cex.lab = 2, mar = c(5, 5, 4, 2) + 0.1)
    plot(x = x[, var], y = y[, var],
         xlab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), xlab),
         ylab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), ylab),
         col = cols, pch = 16)
    legend('topleft', legend = paste('r =', corr), cex = 1.8)
    # lines(loess.smooth(y = y[, var], x = x[, var]), col = 'red')
    # abline(lm(y[, var] ~ x[, var]), col = 'blue')
    abline(h = 0, col = 'grey20')
    abline(v = 0, col = 'grey20')
    abline(a = 0, b = 1, col = 'red')
}



pdf('pdf/scatter_models.pdf', useDingbats = FALSE)
comp_log(outFeat$DLPFC$gene, outFeat_cell$DLPFC$gene, 'DLPFC', 'DLPFC - RNA fraction adj.', var = 't')
comp_log(outFeat$HIPPO$gene, outFeat_cell$HIPPO$gene, 'HIPPO', 'HIPPO - RNA fraction adj.', var = 't')
comp_log(outFeat_cell$DLPFC$gene, outFeat_cell$HIPPO$gene, 'DLPFC - RNA fraction adj.', 'HIPPO - RNA fraction adj.', var = 't')

comp_log(outFeat$DLPFC$gene, outFeat_cell$DLPFC$gene, 'DLPFC', 'DLPFC - RNA fraction adj.')
comp_log(outFeat$HIPPO$gene, outFeat_cell$HIPPO$gene, 'HIPPO', 'HIPPO - RNA fraction adj.')
comp_log(outFeat_cell$DLPFC$gene, outFeat_cell$HIPPO$gene, 'DLPFC - RNA fraction adj.', 'HIPPO - RNA fraction adj.')
dev.off()


fact_vect <- function(x) {
    factor(x, levels = c('FALSE', 'TRUE'))
}

de_tabs <- map2(outFeat, outFeat_cell, ~ map2(.x, .y, 
    ~ table('DE original' = fact_vect(.x$adj.P.Val < 0.05), 'DE adj CellProp' = fact_vect(.y$adj.P.Val < 0.05))
))
de_tabs <- do.call(c, de_tabs)

map(de_tabs, addmargins)
# $DLPFC.gene
#            DE adj CellProp
# DE original FALSE  TRUE   Sum
#       FALSE 24407     0 24407
#       TRUE    217    28   245
#       Sum   24624    28 24652
#
# $DLPFC.exon
#            DE adj CellProp
# DE original  FALSE   TRUE    Sum
#       FALSE 396142      1 396143
#       TRUE     411     29    440
#       Sum   396553     30 396583
#
# $DLPFC.jxn
#            DE adj CellProp
# DE original  FALSE   TRUE    Sum
#       FALSE 297143      1 297144
#       TRUE      27     10     37
#       Sum   297170     11 297181
#
# $DLPFC.tx
#            DE adj CellProp
# DE original FALSE  TRUE   Sum
#       FALSE 92726     0 92726
#       TRUE      6     0     6
#       Sum   92732     0 92732
#
# $HIPPO.gene
#            DE adj CellProp
# DE original FALSE  TRUE   Sum
#       FALSE 24604     0 24604
#       TRUE     46     2    48
#       Sum   24650     2 24652
#
# $HIPPO.exon
#            DE adj CellProp
# DE original  FALSE   TRUE    Sum
#       FALSE 396379      7 396386
#       TRUE     192      5    197
#       Sum   396571     12 396583
#
# $HIPPO.jxn
#            DE adj CellProp
# DE original  FALSE   TRUE    Sum
#       FALSE 297123     17 297140
#       TRUE      27     14     41
#       Sum   297150     31 297181
#
# $HIPPO.tx
#            DE adj CellProp
# DE original FALSE  TRUE   Sum
#       FALSE 92732     0 92732
#       TRUE      0     0     0
#       Sum   92732     0 92732

options(width = 120)
map_dbl(de_tabs, getOR)
# DLPFC.gene DLPFC.exon  DLPFC.jxn   DLPFC.tx HIPPO.gene HIPPO.exon  HIPPO.jxn   HIPPO.tx
#        Inf  27951.625 110052.963        NaN        Inf   1474.624   9062.575        NaN

set.seed(20190326)
de_chi <- map(de_tabs, chisq.test, simulate.p.value = TRUE, B = 10000)
map_dbl(de_chi, ~ .x$p.value)
# DLPFC.gene DLPFC.exon  DLPFC.jxn   DLPFC.tx HIPPO.gene HIPPO.exon  HIPPO.jxn   HIPPO.tx
# 0.00009999 0.00009999 0.00009999        NaN 0.00019998 0.00009999 0.00009999        NaN



de_bias <- map2_dfr(outFeat, outFeat_cell, ~ map2_dfr(.x, .y,
        ~ data.frame(
            ratio_logfc = .y$logFC / .x$logFC,
            bias_logfc = abs(.y$logFC - .x$logFC) / abs(.x$logFC) * 100,
            ratio_t = .y$t / .x$t,
            bias_t = abs(.y$t - .x$t) / abs(.x$t) * 100,
            de_original = .x$adj.P.Val < 0.05,
            de_cell = .y$adj.P.Val < 0.05,
            stringsAsFactors = FALSE
        )
    )
)
de_bias$region <- rep(names(outFeat), each = nrow(de_bias) / 2)
de_bias$feature <- rep(rep(names(outFeat$DLPFC), map_int(outFeat$DLPFC, nrow)), 2)


# table(de_bias$de_original, de_bias$de_cell, de_bias$region, de_bias$feature)

sts_logfc <- boxplot.stats(de_bias$bias_logfc)$stats
sts_t <- boxplot.stats(de_bias$bias_t)$stats
features <- c('gene', 'exon', 'jxn', 'tx')

pdf('pdf/bias.pdf', useDingbats = TRUE, width = 16)
ggplot(de_bias, aes(x = region, y = bias_logfc, fill = region)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(~ factor(feature, levels = features)) +
    scale_fill_manual(values = c('DLPFC' = 'darkgoldenrod2', 'HIPPO' = 'steelblue1')) +
    xlab('Region') +
    ylab('Percent absolute bias on log2 FC') +
    coord_cartesian(ylim = c(0, max(sts_logfc) * 1.2)) +
    theme_bw(base_size = 30)
ggplot(de_bias, aes(x = region, y = bias_t, fill = region)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(~ factor(feature, levels = features)) +
    scale_fill_manual(values = c('DLPFC' = 'darkgoldenrod2', 'HIPPO' = 'steelblue1')) +
    xlab('Region') +
    ylab('Percent absolute bias on t-stats') +
    coord_cartesian(ylim = c(0, max(sts_t) * 1.2)) +
    theme_bw(base_size = 30)
dev.off()

pdf('pdf/bias_de_original.pdf', useDingbats = TRUE, width = 16)
ggplot(subset(de_bias, de_original), aes(x = region, y = bias_logfc, fill = region)) +
    geom_boxplot() +
    facet_grid(~ factor(feature, levels = features)) +
    scale_fill_manual(values = c('DLPFC' = 'darkgoldenrod2', 'HIPPO' = 'steelblue1')) +
    xlab('Region') +
    ylab('Percent absolute bias on log2 FC') +
    theme_bw(base_size = 30)
ggplot(subset(de_bias, de_original), aes(x = region, y = bias_t, fill = region)) +
    geom_boxplot() +
    facet_grid(~ factor(feature, levels = features)) +
    scale_fill_manual(values = c('DLPFC' = 'darkgoldenrod2', 'HIPPO' = 'steelblue1')) +
    xlab('Region') +
    ylab('Percent absolute bias on t-stats') +
    theme_bw(base_size = 30)
dev.off()


## Explore ratio of adj/original t-stats and log2 FC
sts_ratio_logfc <- boxplot.stats(abs(de_bias$ratio_logfc))$stats
sts_ratio_t <- boxplot.stats(abs(de_bias$ratio_t))$stats

pdf('pdf/ratio.pdf', useDingbats = TRUE, width = 16)
ggplot(de_bias, aes(x = region, y = abs(ratio_logfc), fill = region)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(~ factor(feature, levels = features)) +
    scale_fill_manual(values = c('DLPFC' = 'darkgoldenrod2', 'HIPPO' = 'steelblue1')) +
    xlab('Region') +
    ylab('Abs. ratio (adj / original) on log2 FC') +
    coord_cartesian(ylim = c(0, max(sts_ratio_logfc) * 1.2)) +
    theme_bw(base_size = 30) +
    geom_hline(aes(yintercept = 1), linetype = 'dashed', color = 'red')
ggplot(de_bias, aes(x = region, y = abs(ratio_t), fill = region)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(~ factor(feature, levels = features)) +
    scale_fill_manual(values = c('DLPFC' = 'darkgoldenrod2', 'HIPPO' = 'steelblue1')) +
    xlab('Region') +
    ylab('Abs. ratio (adj / original) on t-stats') +
    coord_cartesian(ylim = c(0, max(sts_ratio_t) * 1.2)) +
    theme_bw(base_size = 30) +
    geom_hline(aes(yintercept = 1), linetype = 'dashed', color = 'red')
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
#  date     2019-03-26
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package          * version   date       lib source
#  assertthat         0.2.1     2019-03-21 [2] CRAN (R 3.5.1)
#  BiocGenerics       0.28.0    2018-10-30 [1] Bioconductor
#  bitops             1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  cli                1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout         * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace         1.4-1     2019-03-18 [2] CRAN (R 3.5.1)
#  crayon             1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  digest             0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr              0.8.0.1   2019-02-15 [1] CRAN (R 3.5.1)
#  GenomeInfoDb       1.18.2    2019-02-12 [1] Bioconductor
#  GenomeInfoDbData   1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges      1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2          * 3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue               1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
#  gtable             0.3.0     2019-03-25 [2] CRAN (R 3.5.1)
#  htmltools          0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets        1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv             1.5.0     2019-03-15 [2] CRAN (R 3.5.1)
#  IRanges            2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab         * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  jsonlite           1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  labeling           0.3       2014-08-23 [2] CRAN (R 3.5.0)
#  later              0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice            0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval           0.2.2     2019-03-15 [2] CRAN (R 3.5.1)
#  limma              3.38.3    2018-12-02 [1] Bioconductor
#  magrittr           1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  munsell            0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  pillar             1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig          2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr               1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises           1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr            * 0.3.2     2019-03-15 [2] CRAN (R 3.5.1)
#  R6                 2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  rafalib          * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer       1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp               1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl              1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  reshape2           1.4.3     2017-12-11 [2] CRAN (R 3.5.0)
#  rlang              0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote            * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors          0.20.1    2018-11-09 [1] Bioconductor
#  scales             1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented          0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr              0.13      2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo      * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  stringi            1.4.3     2019-03-12 [2] CRAN (R 3.5.1)
#  stringr            1.4.0     2019-02-10 [1] CRAN (R 3.5.1)
#  tibble             2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect         0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr              2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun               0.5       2019-02-20 [1] CRAN (R 3.5.1)
#  XVector            0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc           1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
