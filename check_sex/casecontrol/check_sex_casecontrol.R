## Based on https://github.com/LieberInstitute/brainseq_phase2/blob/master/development/sz_effect_devel_overlap.R
library('GenomicRanges')
library('SummarizedExperiment')
library('purrr')
library('jaffelab')
library('ggplot2')
library('dplyr')
library('sessioninfo')

## For the chr info
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata", verbose = TRUE)


## List the files for each brain region
## either the original ones (DE by Dx) or the
## sex ones (DE by sex)
f_sex <- c(
    'DLPFC' = 'rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda',
    'HIPPO' = 'rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda'
)
stopifnot(all(file.exists(f_sex)))
f_ori <- paste0('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/', f_sex)
names(f_ori) <- names(f_sex)
stopifnot(all(file.exists(f_ori)))

## Create a similar object to:
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/development/sz_effect_devel_overlap.R#L5-L10
load_stats <- function(stats_files) {
    final <- do.call(c, map2(stats_files, names(stats_files), function(f, region) {
        message(paste(Sys.time(), 'loading', f))
        load(f, verbose = TRUE)
        
        stopifnot(identical(rownames(outGene), names(rowRanges(rse_gene))))
        stopifnot(identical(rownames(outExon), names(rowRanges(rse_exon))))
        stopifnot(identical(rownames(outJxn), names(rowRanges(rse_jxn))))
        stopifnot(identical(rownames(outTx), names(rowRanges(rse_tx))))
        
        res <- list(
            'Gene' = cbind(outGene, chr = seqnames(rowRanges(rse_gene))),
            'Exon' = cbind(outExon, chr = seqnames(rowRanges(rse_exon))),
            'Jxn' = cbind(outJxn, chr = seqnames(rowRanges(rse_jxn))),
            'Tx' = cbind(outTx, chr = seqnames(rowRanges(rse_tx)))
        )
        names(res) <- paste0(region, '_', names(res))
        
        ## Add Bonferroni
        res <- map(res, function(x) {
            x$p_bonf = p.adjust(x$P.Value, method = "bonf")
        	return(x)
        })
        return(res)
    }))
    names(final) <- gsub('.*\\.', '', names(final))
    return(final)
}

sz_stats <- load_stats(f_ori)
sex_stats <- load_stats(f_sex)

map_int(sz_stats, nrow)
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#      24652     396583     297181      92732      24652     396583     297181
#   HIPPO_Tx
#      92732

stopifnot(identical(map_int(sz_stats, nrow), map_int(sex_stats, nrow)))

## Line up
sz_stats <- map2(sz_stats, sex_stats, ~ .x[rownames(.y), ])

## Add sex stats
sz_stats <- map2(sz_stats, sex_stats, function(sz, sex) {
    sz$sexReg <- sex$adj.P.Val < 0.05
    sz$sexM <- sex$adj.P.Val < 0.05 & sex$t > 0
    sz$sexF <- sex$adj.P.Val < 0.05 & sex$t < 0
    sz$sext <- sex$t
    return(sz)
})

## Without chrX and Y
sz_stats_nosex <- map(sz_stats, ~ .x[!.x$chr %in% c('chrX', 'chrY'), ])
map_int(sz_stats_nosex, nrow)
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#      23777     384980     287599      90119      23777     384980     287599
#   HIPPO_Tx
#      90119
map_int(sz_stats_nosex, nrow) / map_int(sz_stats, nrow) * 100
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#   96.45059   97.07426   96.77570   97.18220   96.45059   97.07426   96.77570
#   HIPPO_Tx
#   97.18220


map_int(sz_stats, ~ sum(.x$sexReg))
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#        282       1576       1517        340        116       1261       1109
#   HIPPO_Tx
#        356
map_int(sz_stats_nosex, ~ sum(.x$sexReg))
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#        202        537        610        106         45        274        301
#   HIPPO_Tx
#        128

map_int(sz_stats_nosex, ~ sum(.x$sexReg)) / map_int(sz_stats, ~ sum(.x$sexReg)) * 100
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#   71.63121   34.07360   40.21094   31.17647   38.79310   21.72879   27.14157
#   HIPPO_Tx
#   35.95506

save(sz_stats, sz_stats_nosex, file = 'rdas/sz_stats.Rdata')


compute_overlap <- function(stats) {
    table(
        'SCZD DE' = factor(stats$adj.P.Val < 0.05, levels = c('FALSE', 'TRUE')),
        'Sex DE' = factor(stats$sexReg, levels = c('FALSE', 'TRUE'))
    )
}

## Overlap of SCZD DE features and Sex DE features
map(sz_stats, ~ addmargins(compute_overlap(.x)))
# $DLPFC_Gene
#        Sex DE
# SCZD DE FALSE  TRUE   Sum
#   FALSE 24130   277 24407
#   TRUE    240     5   245
#   Sum   24370   282 24652
#
# $DLPFC_Exon
#        Sex DE
# SCZD DE  FALSE   TRUE    Sum
#   FALSE 394567   1576 396143
#   TRUE     440      0    440
#   Sum   395007   1576 396583
#
# $DLPFC_Jxn
#        Sex DE
# SCZD DE  FALSE   TRUE    Sum
#   FALSE 295628   1516 297144
#   TRUE      36      1     37
#   Sum   295664   1517 297181
#
# $DLPFC_Tx
#        Sex DE
# SCZD DE FALSE  TRUE   Sum
#   FALSE 92386   340 92726
#   TRUE      6     0     6
#   Sum   92392   340 92732
#
# $HIPPO_Gene
#        Sex DE
# SCZD DE FALSE  TRUE   Sum
#   FALSE 24488   116 24604
#   TRUE     48     0    48
#   Sum   24536   116 24652
#
# $HIPPO_Exon
#        Sex DE
# SCZD DE  FALSE   TRUE    Sum
#   FALSE 395125   1261 396386
#   TRUE     197      0    197
#   Sum   395322   1261 396583
#
# $HIPPO_Jxn
#        Sex DE
# SCZD DE  FALSE   TRUE    Sum
#   FALSE 296033   1107 297140
#   TRUE      39      2     41
#   Sum   296072   1109 297181
#
# $HIPPO_Tx
#        Sex DE
# SCZD DE FALSE  TRUE   Sum
#   FALSE 92376   356 92732
#   TRUE      0     0     0
#   Sum   92376   356 92732


## Corresponding chisq p-values
map_dbl(sz_stats, ~ chisq.test(compute_overlap(.x), simulate.p.value = TRUE, B = 1e4)$p.value)
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
# 0.20927907 0.27497250 0.16818318 1.00000000 1.00000000 0.65663434 0.01189881
#   HIPPO_Tx
#        NaN

## Odds ratio
map_dbl(sz_stats, ~ getOR(compute_overlap(.x)))


## Now exclude autosome features
map(sz_stats_nosex, ~ addmargins(compute_overlap(.x)))
# $DLPFC_Gene
#        Sex DE
# SCZD DE FALSE  TRUE   Sum
#   FALSE 23345   197 23542
#   TRUE    230     5   235
#   Sum   23575   202 23777
#
# $DLPFC_Exon
#        Sex DE
# SCZD DE  FALSE   TRUE    Sum
#   FALSE 384030    537 384567
#   TRUE     413      0    413
#   Sum   384443    537 384980
#
# $DLPFC_Jxn
#        Sex DE
# SCZD DE  FALSE   TRUE    Sum
#   FALSE 286955    609 287564
#   TRUE      34      1     35
#   Sum   286989    610 287599
#
# $DLPFC_Tx
#        Sex DE
# SCZD DE FALSE  TRUE   Sum
#   FALSE 90007   106 90113
#   TRUE      6     0     6
#   Sum   90013   106 90119
#
# $HIPPO_Gene
#        Sex DE
# SCZD DE FALSE  TRUE   Sum
#   FALSE 23684    45 23729
#   TRUE     48     0    48
#   Sum   23732    45 23777
#
# $HIPPO_Exon
#        Sex DE
# SCZD DE  FALSE   TRUE    Sum
#   FALSE 384509    274 384783
#   TRUE     197      0    197
#   Sum   384706    274 384980
#
# $HIPPO_Jxn
#        Sex DE
# SCZD DE  FALSE   TRUE    Sum
#   FALSE 287259    299 287558
#   TRUE      39      2     41
#   Sum   287298    301 287599
#
# $HIPPO_Tx
#        Sex DE
# SCZD DE FALSE  TRUE   Sum
#   FALSE 89991   128 90119
#   TRUE      0     0     0
#   Sum   89991   128 90119

map_dbl(sz_stats_nosex, ~ chisq.test(foo(.x), simulate.p.value = TRUE, B = 1e4)$p.value)
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
# 0.04919508 0.67993201 0.06859314 1.00000000 1.00000000 1.00000000 0.00089991
#   HIPPO_Tx
#        NaN

## Odds ratio
map_dbl(sz_stats, ~ getOR(compute_overlap(.x)))
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#   1.814832   0.000000   5.416813   0.000000   0.000000   0.000000  13.713803
#   HIPPO_Tx
#        NaN
map_dbl(sz_stats_nosex, ~ getOR(compute_overlap(.x)))
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#   2.576142   0.000000  13.858543   0.000000   0.000000   0.000000  49.268330
#   HIPPO_Tx
#        NaN

## Correlation across t-statistics
map_dbl(sz_stats, ~ cor(.x$t, .x$sext))
map_dbl(sz_stats_nosex, ~ cor(.x$t, .x$sext))


plot_t_stats <- function(set, filename, stats) {
    i <- c(1, 5) + case_when(
        set == 'gene' ~ 0,
        set == 'exon' ~ 1,
        set == 'jxn' ~ 2,
        set == 'tx' ~ 3
    )
    pdf(file = filename, useDingbats = FALSE, width = 14)
    map2(stats[i], names(stats)[i], function(df, title) {
        df <- as.data.frame(df)
        df$FDR <- df$adj.P.Val
        print(
            ggplot(df, aes(x = t, y = sext, color = FDR < 0.05)) +
                geom_point() +
                facet_grid(~ sexReg) +
                theme_bw(base_size = 30) +
                ggtitle(title) +
                ylab('Sex t-statistic') + 
                xlab('SCZD case vs control t-statistic') +
                labs(caption = 'Separated by Sex FDR < 0.05 status')
                
        )
    })
    dev.off()
}


features <- c('gene', 'exon', 'jxn', 'tx')

map(features, ~
    plot_t_stats(.x, paste0('pdf/sczd_vs_sex_t_', .x, '.pdf'), sz_stats)
)
map(features, ~
    plot_t_stats(.x, paste0('pdf/sczd_vs_sex_t_', .x, '_autosomal.pdf'), sz_stats_nosex)
)


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
#  date     2019-03-07
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
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                * 0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
#  fansi                  0.4.0     2018-10-05 [1] CRAN (R 3.5.1)
#  GenomeInfoDb         * 1.18.1    2018-11-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2              * 3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.0     2018-07-17 [1] CRAN (R 3.5.1)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab             * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  labeling               0.3       2014-08-23 [2] CRAN (R 3.5.0)
#  later                  0.7.5     2018-09-18 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma                  3.38.3    2018-12-02 [1] Bioconductor
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
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.5.0)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr                  0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  stringi                1.2.4     2018-07-20 [2] CRAN (R 3.5.1)
#  stringr                1.3.1     2018-05-10 [1] CRAN (R 3.5.0)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  utf8                   1.1.4     2018-05-24 [1] CRAN (R 3.5.0)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
