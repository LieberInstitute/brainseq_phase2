library('sessioninfo')
library('SummarizedExperiment')
library('purrr')
library('jaffelab')

## Original code:
## https://github.com/LieberInstitute/brain-epigenomics/blob/2018d6a01856f758b592839486419203f4121962/meth_vs_expr/venn_and_summary_using_near.R#L296-L303
load('/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/sorted_nuclear_RNA/DE_limma_results_objects.rda', verbose = TRUE)
co <- grep('CellType', colnames(fit_gene_combined$p.value))
print(colnames(fit_gene_combined$p.value)[co])
top <- topTable(fit_gene_combined, coef = co, sort.by = 'none',
    number = nrow(fit_gene_combined$p.value))

## Add gene annotation info to the top results
top <- cbind(top, geneMap[rownames(top), ])
dim(top)
# [1] 46367    20
length(unique(top$ensemblID))
# [1] 46367
dim(geneMap)
# [1] 60252    14

## So, not all genes were tested in this analysis. Likely due to low expression.

## Get the FDR <5% for each cell type
split_top <- function(x) {
    top_sig <- map(split(x, sign(x$t)), ~ .x[.x$adj.P.Val < 0.05, ])
    names(top_sig) <- c('-1' = 'Neuron', '1' = 'Glia')[names(top_sig)]
    return(top_sig)
}
top_sig <- split_top(top)

## Number of DE genes (FDR<5%) by cell type
map_int(top_sig, nrow)
# Neuron   Glia
#   3473   6521

## Load our interaction eQTLs
load('../eqtl_tables/mergedEqtl_output_interaction_4features.rda', verbose = TRUE)

## Remove \\..* part from the EnsemblGeneIDs
allEqtl$EnsemblGeneID <- gsub('\\..*', '', allEqtl$EnsemblGeneID)

## Number of unique ENSEMBL gene ids we have
ens_unique <- unique(allEqtl$EnsemblGeneID)
length(ens_unique)
# [1] 20522

## How many of them are present in the Neuron vs Glia DE analysis
table(ens_unique %in% top$ensemblID)
# FALSE  TRUE
#   875 19647

## How many of them would map in any case (due to different versions of annotation)
table(ens_unique %in% geneMap$ensemblID)
# FALSE  TRUE
#   283 20239

## Code that made me realize that we needed to remove the \\..* part
## when creating ens_unique
head(ens_unique[!ens_unique %in% top$ensemblID])

## Ok, load our data to filter by ids
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_gene_unfiltered.Rdata', verbose = TRUE)

top_filt <- top[top$ensemblID %in% rowRanges(rse_gene)$ensemblID, ]
dim(top_filt)
# [1] 44338    20

nrow(top_filt) / nrow(top) * 100
# [1] 95.62404

mean(top_filt$adj.P.Val < 0.05) * 100
# [1] 22.2112

## It was still better to use the EnsemblID over the Gencode ID
table(gsub('_.*', '', top$gencodeID) %in% rowRanges(rse_gene)$gencodeID)
# FALSE  TRUE
#  2773 43594

top_sig_filt <- split_top(top_filt)

## Number of DE genes (FDR<5%) by cell type
## that match ensembl ids to our data
map_int(top_sig_filt, nrow)
# Neuron   Glia
#   3439   6409
sum(map_int(top_sig_filt, nrow))
# [1] 9848

## In percent of the original list (before filtering by matching ids)
map_int(top_sig_filt, nrow) / map_int(top_sig, nrow) * 100
#   Neuron     Glia
# 99.02102 98.28247

sum(map_int(top_sig_filt, nrow)) / sum(map_int(top_sig, nrow)) * 100
# [1] 98.53912

## Number of DE genes by cell type in our eQTLs
map_int(top_sig_filt, ~ sum(.x$ensemblID %in% ens_unique))
# Neuron   Glia
#   2462   4333
sum(map_int(top_sig_filt, ~ sum(.x$ensemblID %in% ens_unique)))
# [1] 6795

## in percent of the DE genes
sum(map_int(top_sig_filt, ~ sum(.x$ensemblID %in% ens_unique))) / sum(map_int(top_sig_filt, nrow)) * 100
# [1] 68.99878
map_int(top_sig_filt, ~ sum(.x$ensemblID %in% ens_unique)) / map_int(top_sig_filt, nrow) * 100
#   Neuron     Glia
# 71.59058 67.60805




## Now look at only the genes that have a eQTL at FDR <1%
ens_unique_sig <- unique(allEqtl$EnsemblGeneID[allEqtl$FDR < 0.01])
length(ens_unique_sig)
# [1] 1798

map_int(top_sig_filt, ~ sum(.x$ensemblID %in% ens_unique_sig))
# Neuron   Glia
#    150    596
sum(map_int(top_sig_filt, ~ sum(.x$ensemblID %in% ens_unique_sig)))
# [1] 746
sum(map_int(top_sig_filt, ~ sum(.x$ensemblID %in% ens_unique_sig))) / sum(map_int(top_sig_filt, nrow)) * 100
# [1] 7.575142
map_int(top_sig_filt, ~ sum(.x$ensemblID %in% ens_unique_sig)) / map_int(top_sig_filt, nrow) * 100
#   Neuron     Glia
# 4.361733 9.299423

## Percent of unique eQTL (FDR<1%) genes (across any of the 4 expression features)
## that is in the list of significant DE by Neun vs Glia
stopifnot(identical(
    sum(ens_unique_sig %in% top_filt$ensemblID[top_filt$adj.P.Val < 0.05]),
    sum(map_int(top_sig_filt, ~ sum(.x$ensemblID %in% ens_unique_sig)))
))
sum(ens_unique_sig %in% top_filt$ensemblID[top_filt$adj.P.Val < 0.05]) / length(ens_unique_sig) * 100
# [1] 41.49055

prop.test(
    sum(ens_unique_sig %in% top_filt$ensemblID[top_filt$adj.P.Val < 0.05]),
    n = length(ens_unique_sig),
    p = mean(top_filt$adj.P.Val < 0.05),
    alt = 'greater'
)
#     1-sample proportions test with continuity correction
#
# X-squared = 385.68, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is greater than 0.222112
# 95 percent confidence interval:
#  0.3956588 1.0000000
# sample estimates:
#         p
# 0.4149055
prop.test(
    sum(ens_unique_sig %in% top_filt$ensemblID[top_filt$adj.P.Val < 0.05]),
    n = length(ens_unique_sig),
    p = mean(top_filt$adj.P.Val < 0.05),
    alt = 'greater'
)$p.value
# 3.601288e-86

## Exclude transcript eQTLs
ens_unique_sig_notx <- unique(allEqtl$EnsemblGeneID[allEqtl$FDR < 0.01 & allEqtl$Type != 'Tx'])
stopifnot(sum(allEqtl$FDR < 0.01 & allEqtl$Type != 'Tx') == 205618)

## This number is more than the 1,484 from Figure 4C
## which is based on Gencode IDs instead of Ensembl IDs
## See https://github.com/LieberInstitute/brainseq_phase2/blob/62fa5770071f07552bae5cde62ac0b862daefc68/development/explore_limma_dev.R#L363-L409
length(ens_unique_sig_notx)
# [1] 1632
sum(ens_unique_sig_notx %in% top_filt$ensemblID[top_filt$adj.P.Val < 0.05])
# [1] 684
sum(ens_unique_sig_notx %in% top_filt$ensemblID[top_filt$adj.P.Val < 0.05]) / length(ens_unique_sig_notx) * 100
# [1] 41.91176

map_int(top_sig_filt, ~ sum(.x$ensemblID %in% ens_unique_sig_notx))
# Neuron   Glia
#    137    547
sum(map_int(top_sig_filt, ~ sum(.x$ensemblID %in% ens_unique_sig_notx))) / sum(map_int(top_sig_filt, nrow)) * 100
# [1] 6.945573
map_int(top_sig_filt, ~ sum(.x$ensemblID %in% ens_unique_sig_notx)) / map_int(top_sig_filt, nrow) * 100
#   Neuron     Glia
# 3.983716 8.534873

prop.test(
    sum(ens_unique_sig_notx %in% top_filt$ensemblID[top_filt$adj.P.Val < 0.05]),
    n = length(ens_unique_sig_notx),
    p = mean(top_filt$adj.P.Val < 0.05),
    alt = 'greater'
)
#     1-sample proportions test with continuity correction
#
# X-squared = 365.46, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is greater than 0.222112
# 95 percent confidence interval:
#  0.3988739 1.0000000
# sample estimates:
#         p
# 0.4191176
prop.test(
    sum(ens_unique_sig_notx %in% top_filt$ensemblID[top_filt$adj.P.Val < 0.05]),
    n = length(ens_unique_sig_notx),
    p = mean(top_filt$adj.P.Val < 0.05),
    alt = 'greater'
)$p.value
# [1] 9.125033e-82


eqtl_vs_neunglia <- function(x) {
    with(x,
        table(
            'Neuron/Glia DE' = EnsemblGeneID %in% top_filt$ensemblID[top_filt$adj.P.Val < 0.05],
            'eQTL FDR <1%' = FDR < 0.01
        )
    )
}

## From https://github.com/LieberInstitute/brainseq_phase2/blob/master/eQTL_full/check_cauc/cauc_eqtl_replication.R#L112-L119
or_chisq <- function(x) {
    x <- x[1:2, 1:2]
    list(
        'OR' = jaffelab::getOR(x),
        'chisq.test' = chisq.test(x),
        'p.value' = chisq.test(x)$p.value
    )
}

overall <- eqtl_vs_neunglia(allEqtl)
addmargins(overall)
#               eQTL FDR <1%
# Neuron/Glia DE   FALSE    TRUE     Sum
#          FALSE 1073734  154699 1228433
#          TRUE   475564   72814  548378
#          Sum   1549298  227513 1776811
addmargins(overall / sum(overall) * 100)
#               eQTL FDR <1%
# Neuron/Glia DE      FALSE       TRUE        Sum
#          FALSE  60.430400   8.706553  69.136954
#          TRUE   26.765030   4.098016  30.863046
#          Sum    87.195430  12.804570 100.000000
         
or_chisq(overall)
# $OR
# [1] 1.062711
#
# $chisq.test
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  x
# X-squared = 159.21, df = 1, p-value < 2.2e-16
#
#
# $p.value
# [1] 1.680764e-36

## Exclude transcript eQTLs
overall_notx <- eqtl_vs_neunglia(subset(allEqtl, Type != 'Tx'))
addmargins(overall_notx)
#               eQTL FDR <1%
# Neuron/Glia DE   FALSE    TRUE     Sum
#          FALSE  771326  123700  895026
#          TRUE   603839   81918  685757
#          Sum   1375165  205618 1580783
addmargins(overall_notx / sum(overall_notx) * 100)
#               eQTL FDR <1%
# Neuron/Glia DE      FALSE       TRUE        Sum
#          FALSE  48.793920   7.825236  56.619156
#          TRUE   38.198728   5.182115  43.380844
#          Sum    86.992649  13.007351 100.000000
or_chisq(overall_notx)
# $OR
# [1] 0.8459145
#
# $chisq.test
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  x
# X-squared = 1206.4, df = 1, p-value < 2.2e-16
#
#
# $p.value
# [1] 2.462958e-264

## Next, broken down by feature
features <- c('Gene', 'Exon', 'Jxn', 'Tx')
eqtls_by_type <- lapply(split(allEqtl, factor(allEqtl$Type, levels = features)), as.data.frame)

eqtl_tabs <- map(eqtls_by_type, eqtl_vs_neunglia)
map(eqtl_tabs, addmargins)
# $Gene
#               eQTL FDR <1%
# Neuron/Glia DE  FALSE   TRUE    Sum
#          FALSE  59326  25644  84970
#          TRUE   36508  14448  50956
#          Sum    95834  40092 135926
#
# $Exon
#               eQTL FDR <1%
# Neuron/Glia DE  FALSE   TRUE    Sum
#          FALSE 394719  41701 436420
#          TRUE  365588  48222 413810
#          Sum   760307  89923 850230
#
# $Jxn
#               eQTL FDR <1%
# Neuron/Glia DE  FALSE   TRUE    Sum
#          FALSE 317281  56355 373636
#          TRUE  201743  19248 220991
#          Sum   519024  75603 594627
#
# $Tx
#               eQTL FDR <1%
# Neuron/Glia DE  FALSE   TRUE    Sum
#          FALSE 100665  11751 112416
#          TRUE   73468  10144  83612
#          Sum   174133  21895 196028

## In percent
map(eqtl_tabs, ~ addmargins (.x / sum(.x) * 100))
# $Gene
#               eQTL FDR <1%
# Neuron/Glia DE     FALSE      TRUE       Sum
#          FALSE  43.64581  18.86615  62.51196
#          TRUE   26.85873  10.62931  37.48804
#          Sum    70.50454  29.49546 100.00000
#
# $Exon
#               eQTL FDR <1%
# Neuron/Glia DE      FALSE       TRUE        Sum
#          FALSE  46.424967   4.904673  51.329640
#          TRUE   42.998718   5.671642  48.670360
#          Sum    89.423685  10.576315 100.000000
#
# $Jxn
#               eQTL FDR <1%
# Neuron/Glia DE      FALSE       TRUE        Sum
#          FALSE  53.357987   9.477370  62.835357
#          TRUE   33.927655   3.236987  37.164643
#          Sum    87.285643  12.714357 100.000000
#
# $Tx
#               eQTL FDR <1%
# Neuron/Glia DE      FALSE       TRUE        Sum
#          FALSE  51.352358   5.994552  57.346910
#          TRUE   37.478319   5.174771  42.653090
#          Sum    88.830677  11.169323 100.000000

map(eqtl_tabs, or_chisq)
# $Gene
# $Gene$OR
# [1] 0.9155435
#
# $Gene$chisq.test
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  x
# X-squared = 50.995, df = 1, p-value = 9.258e-13
#
#
# $Gene$p.value
# [1] 9.25828e-13
#
#
# $Exon
# $Exon$OR
# [1] 1.248518
#
# $Exon$chisq.test
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  x
# X-squared = 988.25, df = 1, p-value < 2.2e-16
#
#
# $Exon$p.value
# [1] 6.434004e-217
#
#
# $Jxn
# $Jxn$OR
# [1] 0.5371539
#
# $Jxn$chisq.test
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  x
# X-squared = 5081.4, df = 1, p-value < 2.2e-16
#
#
# $Jxn$p.value
# [1] 0
#
#
# $Tx
# $Tx$OR
# [1] 1.182809
#
# $Tx$chisq.test
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  x
# X-squared = 136.08, df = 1, p-value < 2.2e-16
#
#
# $Tx$p.value
# [1] 1.915411e-31

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
#  date     2019-03-19
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
#  jaffelab               0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  later                  0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma                * 3.38.3    2018-12-02 [1] Bioconductor
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                * 0.3.1     2019-03-03 [2] CRAN (R 3.5.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  rafalib                1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
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
