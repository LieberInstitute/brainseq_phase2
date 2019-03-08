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
overlaps <- map(sz_stats, ~ addmargins(compute_overlap(.x)))
overlaps
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
set.seed(20190308)
p_vals <- map_dbl(sz_stats, ~ chisq.test(compute_overlap(.x), simulate.p.value = TRUE, B = 1e4)$p.value)
p_vals
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
# 0.21457854 0.26997300 0.16938306 1.00000000 1.00000000 0.66543346 0.00909909
#   HIPPO_Tx
#        NaN

## Odds ratio
map_dbl(sz_stats, ~ getOR(compute_overlap(.x)))


## Now exclude autosome features
overlaps_auto <- map(sz_stats_nosex, ~ addmargins(compute_overlap(.x)))
overlaps_auto
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

set.seed(20190308)
p_vals_auto <- map_dbl(sz_stats_nosex, ~ chisq.test(compute_overlap(.x), simulate.p.value = TRUE, B = 1e4)$p.value)
p_vals_auto
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
# 0.04969503 0.67283272 0.07049295 1.00000000 1.00000000 1.00000000 0.00079992
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
#  DLPFC_Gene   DLPFC_Exon    DLPFC_Jxn     DLPFC_Tx   HIPPO_Gene   HIPPO_Exon
# 0.047719208  0.036298540  0.032824981  0.016283529 -0.009889695  0.010819731
#   HIPPO_Jxn     HIPPO_Tx
# 0.000168062 -0.015039394
map_dbl(sz_stats_nosex, ~ cor(.x$t, .x$sext))
#   DLPFC_Gene    DLPFC_Exon     DLPFC_Jxn      DLPFC_Tx    HIPPO_Gene
# 0.0943334650  0.0760722723  0.0627606668  0.0391245717  0.0004921848
#   HIPPO_Exon     HIPPO_Jxn      HIPPO_Tx
# 0.0249700703 -0.0055431699  0.0244796302

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


## Make a table
make_table <- function(stats) {
    ov <- map(stats, ~ compute_overlap(.x))
    
    res <- map_dfr(ov,
        ~ as.data.frame(matrix(as.vector(.x), nrow = 1, dimnames = list(1, c('Null_both', 'DE_SCZD', 'DE_Sex', 'DE_both'))))
    )
    res$region <- ss(names(ov), '_', 1)
    res$feature <- tolower(ss(names(ov), '_', 2))
    res$OR <- map_dbl(ov, getOR)
    set.seed(20180308)
    res$pval <- map_dbl(ov, ~ chisq.test(.x, simulate.p.value = TRUE, B = 1e5)$p.value)
    res$pval_bonf <- p.adjust(res$pval, 'bonf')
    res$cor_t <- map_dbl(stats, ~ cor(.x$t, .x$sext))
    
    ## Re-order by feature
    res <- res[c(1,5,2,6,3,7,4,8), ]
    return(res)
}

sz_table <- make_table(sz_stats)
options(width = 120)
sz_table
#   Null_both DE_SCZD DE_Sex DE_both region feature        OR      pval  pval_bonf        cor_t
# 1     24130     240    277       5  DLPFC    gene  1.814832 0.2101379 1.00000000  0.047719208
# 5     24488      48    116       0  HIPPO    gene  0.000000 1.0000000 1.00000000 -0.009889695
# 2    394567     440   1576       0  DLPFC    exon  0.000000 0.2749373 1.00000000  0.036298540
# 6    395125     197   1261       0  HIPPO    exon  0.000000 0.6642134 1.00000000  0.010819731
# 3    295628      36   1516       1  DLPFC     jxn  5.416813 0.1715583 1.00000000  0.032824981
# 7    296033      39   1107       2  HIPPO     jxn 13.713803 0.0104399 0.07307927  0.000168062
# 4     92386       6    340       0  DLPFC      tx  0.000000 1.0000000 1.00000000  0.016283529
# 8     92376       0    356       0  HIPPO      tx       NaN       NaN        NaN -0.015039394

sz_table_nosex <- make_table(sz_stats_nosex)
sz_table_nosex
#   Null_both DE_SCZD DE_Sex DE_both region feature        OR         pval   pval_bonf         cor_t
# 1     23345     230    197       5  DLPFC    gene  2.576142 0.0506594934 0.354616454  0.0943334650
# 5     23684      48     45       0  HIPPO    gene  0.000000 1.0000000000 1.000000000  0.0004921848
# 2    384030     413    537       0  DLPFC    exon  0.000000 0.6751532485 1.000000000  0.0760722723
# 6    384509     197    274       0  HIPPO    exon  0.000000 1.0000000000 1.000000000  0.0249700703
# 3    286955      34    609       1  DLPFC     jxn 13.858543 0.0720592794 0.504414956  0.0627606668
# 7    287259      39    299       2  HIPPO     jxn 49.268330 0.0007699923 0.005389946 -0.0055431699
# 4     90007       6    106       0  DLPFC      tx  0.000000 1.0000000000 1.000000000  0.0391245717
# 8     89991       0    128       0  HIPPO      tx       NaN          NaN         NaN  0.0244796302

sz_all <- matrix(colSums(sz_table[, 1:4]), ncol = 2)
sz_all
#         [,1] [,2]
# [1,] 1614733 6549
# [2,]    1006    8
set.seed(20180308)
chisq.test(sz_all, simulate.p.value = TRUE, B = 1e5)$p.value
# [1] 0.07333927

sz_all_nosex <- matrix(colSums(sz_table_nosex[, 1:4]), ncol = 2)
sz_all_nosex
set.seed(20180308)
#         [,1] [,2]
# [1,] 1569780 2195
# [2,]     967    8
chisq.test(sz_all_nosex, simulate.p.value = TRUE, B = 1e5)$p.value
# [1] 0.0001399986

save(sz_table, sz_table_nosex, file = 'rdas/sz_table.Rdata')


export_table <- function(tab, filename) {
    write.csv(tab, file = paste0('rdas/', filename, '.csv'), row.names = FALSE, quote = FALSE)
}
export_table(sz_table, 'sz_vs_sex_summmary')
export_table(sz_table_nosex, 'sz_vs_sex_summmary_autosomal')

## Show those that overlap
options(width = 300)
feat_ov <- map(sz_stats_nosex, ~ .x[.x$adj.P.Val < 0.05 & .x$sexReg, ])
feat_ov[map_lgl(feat_ov, ~ nrow(.x) > 0)]
# $DLPFC_Gene
#                    Length          gencodeID       ensemblID      gene_type Symbol EntrezID Class  meanExprs NumTx    gencodeTx passExprsCut       logFC  AveExpr         t      P.Value  adj.P.Val          B   chr p_bonf sexReg  sexM sexF      sext
# ENSG00000132854.18   5980 ENSG00000132854.18 ENSG00000132854 protein_coding  KANK4   163782 InGen  0.3680731     4 ENST0000....         TRUE -0.26024883 1.056869 -4.029944 6.824644e-05 0.02209946  1.4449003  chr1      1   TRUE FALSE TRUE -3.892680
# ENSG00000177301.13  12344 ENSG00000177301.13 ENSG00000177301 protein_coding  KCNA2     3737 InGen 22.1524555     5 ENST0000....         TRUE -0.06685712 8.736897 -3.719575 2.317602e-04 0.03548667 -0.2766284  chr1      1   TRUE FALSE TRUE -3.485173
# ENSG00000151690.14   6031 ENSG00000151690.14 ENSG00000151690 protein_coding  MFSD6    54842 InGen 18.9680190    10 ENST0000....         TRUE -0.05382237 7.240094 -3.640553 3.123624e-04 0.04011692 -0.4524693  chr2      1   TRUE FALSE TRUE -3.813113
# ENSG00000105971.14   6926 ENSG00000105971.14 ENSG00000105971 protein_coding   CAV2      858 InGen  1.2640041    15 ENST0000....         TRUE  0.11617901 3.173690  3.670080 2.795711e-04 0.03828881  0.1033548  chr7      1   TRUE FALSE TRUE -3.732529
# ENSG00000177599.12   3032 ENSG00000177599.12 ENSG00000177599 protein_coding ZNF491   126069 InGen  1.8280137     4 ENST0000....         TRUE -0.08623048 2.379584 -3.677721 2.716280e-04 0.03783149  0.2173570 chr19      1   TRUE FALSE TRUE -5.094236
#
# $DLPFC_Jxn
#                            inGencode inGencodeStart inGencodeEnd gencodeGeneID ensemblID Symbol gencodeStrand gencodeTx numTx Class startExon endExon newGeneID newGeneSymbol isFusion meanExprs Length passExprsCut      logFC  AveExpr        t      P.Value  adj.P.Val        B   chr    p_bonf sexReg
# chr14:49853669-49862578(-)     FALSE          FALSE        FALSE          <NA>      <NA>   <NA>          <NA>               0 Novel        NA      NA      <NA>          <NA>    FALSE  12239.69    100         TRUE -0.5996511 9.082305 -4.98405 9.733285e-07 0.02066105 -1.05056 chr14 0.2892547   TRUE
#                             sexM sexF      sext
# chr14:49853669-49862578(-) FALSE TRUE -4.094468
#
# $HIPPO_Jxn
#                   inGencode inGencodeStart inGencodeEnd gencodeGeneID ensemblID Symbol gencodeStrand gencodeTx numTx Class startExon endExon newGeneID newGeneSymbol isFusion meanExprs Length passExprsCut    logFC AveExpr         t      P.Value    adj.P.Val         B  chr       p_bonf sexReg  sexM
# chrM:2514-2761(+)     FALSE          FALSE        FALSE          <NA>      <NA>   <NA>          <NA>               0 Novel        NA      NA      <NA>          <NA>    FALSE  4732.402    100         TRUE 4.437348 5.98525 10.608046 1.194987e-22 1.775637e-17 37.555741 chrM 3.551275e-17   TRUE FALSE
# chrM:2579-2761(+)     FALSE          FALSE        FALSE          <NA>      <NA>   <NA>          <NA>               0 Novel        NA      NA      <NA>          <NA>    FALSE  6053.336    100         TRUE 1.414241 8.04210  5.877925 1.068118e-08 5.290406e-04  5.730943 chrM 3.174243e-03   TRUE  TRUE
#                    sexF      sext
# chrM:2514-2761(+)  TRUE -5.246492
# chrM:2579-2761(+) FALSE  4.664882

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
