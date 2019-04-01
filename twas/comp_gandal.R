library('sessioninfo')
library('readxl')
library('ggplot2')
library('patchwork') # devtools::install_github("thomasp85/patchwork")
library('purrr')

## Read Gandal et al results
## http://science.sciencemag.org/content/362/6420/eaat8127


gandal <- read_xlsx('aat8127_Table_S4.xlsx', sheet = 3)

## Fix some entries
gandal[5551 - 1, 20] <- 7.71236473158186e-321
gandal[5875 - 1, 20] <- 1.77863632502849e-322
gandal[9464 - 1, 20] <- 1.76999017622627e-319
gandal[12786 - 1, 20] <- 8.80999999191014e-316
gandal[13232 - 1, 20] <- 6.8399999999901e-313

dim(gandal)
# [1] 14159    31
length(unique(gandal$GeneID))
# [1] 14159

## Load our data (PGC2 + CLOZUK)
load('rda/tt_objects.Rdata', verbose = TRUE)


twas_z_gandal <- subset(region_twas_z, feature == 'gene')
m <- match(gsub('\\..*', '', twas_z_gandal$geneid), gandal$GeneID)
table(is.na(m))
# FALSE  TRUE
#  5257  1410
twas_z_gandal$GANDAL <- gandal$TWAS.Z[m]
twas_z_gandal$GANDAL[is.na(m)] <- 0
twas_z_gandal$in_gandal <- !is.na(m)
twas_z_gandal$TWAS.P_GANDAL <- gandal$TWAS.P[m]
twas_z_gandal$TWAS.FDR_GANDAL <- p.adjust(gandal$TWAS.P, 'fdr')[m]
twas_z_gandal$TWAS.Bonf_GANDAL <- gandal$TWAS.Bonferroni[m]

p_comp <- function(df, xvar, yvar, x, y) {
    result <- rep('None', nrow(df))
    result[df[, xvar] < 0.05] <- x
    result[df[, yvar] < 0.05] <- y
    result[df[, xvar] < 0.05 & df[, yvar] < 0.05] <- 'Both'
    
    result <- factor(result, levels = c('None', x, y, 'Both'))
    print(table(result, useNA = 'ifany'))
    return(result)
}

## Check it works
x <- p_comp(twas_z_gandal, 'TWAS.FDR_DLPFC', 'TWAS.FDR_HIPPO', 'DLPFC', 'HIPPO')
stopifnot(identical(twas_z_gandal$FDR.5perc, x))

twas_z_gandal$FDR_GANDAL_DLPFC <- p_comp(twas_z_gandal, 'TWAS.FDR_GANDAL', 'TWAS.FDR_DLPFC', 'GANDAL', 'DLPFC')
# None GANDAL  DLPFC   Both
# 6055    206    183    223
twas_z_gandal$FDR_GANDAL_HIPPO <- p_comp(twas_z_gandal, 'TWAS.FDR_GANDAL', 'TWAS.FDR_HIPPO', 'GANDAL', 'HIPPO')
# None GANDAL  HIPPO   Both
# 6106    291    132    138

twas_z_gandal$Bonf_GANDAL_DLPFC <- p_comp(twas_z_gandal, 'TWAS.Bonf_GANDAL', 'TWAS.Bonf_DLPFC', 'GANDAL', 'DLPFC')
# None GANDAL  DLPFC   Both
# 6548     38     48     33
twas_z_gandal$Bonf_GANDAL_HIPPO <- p_comp(twas_z_gandal, 'TWAS.Bonf_GANDAL', 'TWAS.Bonf_HIPPO', 'GANDAL', 'HIPPO')
# None GANDAL  HIPPO   Both
# 6558     53     38     18

cols <- c('None' = 'grey80', 'DLPFC' = 'dark orange', 'HIPPO' = 'skyblue3', 'Both' = 'purple', 'GANDAL' = 'springgreen4')
save(twas_z_gandal, cols, file = 'rda/twas_z_gandal.Rdata')

pdf('pdf/twas_z_gandal.pdf', useDingbats = FALSE, width = 24, height = 14)
ggplot(twas_z_gandal,
    aes(x = DLPFC, y = GANDAL, color = FDR_GANDAL_DLPFC, shape = in_gandal)) +
    geom_point() +
    facet_grid(BEST.GWAS.status ~ feature) +
    coord_fixed() +
    theme_bw(base_size = 30) +
    ggtitle('TWAS Z: DLPFC vs Gandal et al') +
    scale_color_manual(values = cols) +
ggplot(twas_z_gandal,
    aes(x = HIPPO, y = GANDAL, color = FDR_GANDAL_HIPPO, shape = in_gandal)) +
    geom_point() +
    facet_grid(BEST.GWAS.status ~ feature) +
    coord_fixed() +
    theme_bw(base_size = 30) +
    ggtitle('TWAS Z: HIPPO vs Gandal et al') +
    scale_color_manual(values = cols)
    
ggplot(twas_z_gandal,
    aes(x = DLPFC, y = GANDAL, color = Bonf_GANDAL_DLPFC, shape = in_gandal)) +
    geom_point() +
    facet_grid(BEST.GWAS.status ~ feature) +
    coord_fixed() +
    theme_bw(base_size = 30) +
    ggtitle('TWAS Z: DLPFC vs Gandal et al') +
    scale_color_manual(values = cols) +
ggplot(twas_z_gandal,
    aes(x = HIPPO, y = GANDAL, color = Bonf_GANDAL_HIPPO, shape = in_gandal)) +
    geom_point() +
    facet_grid(BEST.GWAS.status ~ feature) +
    coord_fixed() +
    theme_bw(base_size = 30) +
    ggtitle('TWAS Z: HIPPO vs Gandal et al') +
    scale_color_manual(values = cols)
dev.off()


tabs <- list(
    'FDR' = with(twas_z_gandal, table(FDR_GANDAL_DLPFC, FDR_GANDAL_HIPPO)),
    'Bonf' = with(twas_z_gandal, table(Bonf_GANDAL_DLPFC, Bonf_GANDAL_HIPPO))
)

lapply(tabs, addmargins)
# $FDR
#                 FDR_GANDAL_HIPPO
# FDR_GANDAL_DLPFC None GANDAL HIPPO Both  Sum
#           None   5963      0    92    0 6055
#           GANDAL    0    166     0   40  206
#           DLPFC   143      0    40    0  183
#           Both      0    125     0   98  223
#           Sum    6106    291   132  138 6667
#
# $Bonf
#                  Bonf_GANDAL_HIPPO
# Bonf_GANDAL_DLPFC None GANDAL HIPPO Both  Sum
#            None   6525      0    23    0 6548
#            GANDAL    0     32     0    6   38
#            DLPFC    33      0    15    0   48
#            Both      0     21     0   12   33
#            Sum    6558     53    38   18 6667

32+12+21+6
# [1] 71

tabs2 <- list(
    'DLPFC' = with(twas_z_gandal, table('DLPFC' = TWAS.FDR_DLPFC < 0.05, 'GANDAL' = TWAS.P_GANDAL < 0.05, useNA = 'ifany')),
    'HIPPO' = with(twas_z_gandal, table('HIPPO' = TWAS.FDR_HIPPO < 0.05, 'GANDAL'= TWAS.P_GANDAL < 0.05, useNA = 'ifany'))
)

lapply(tabs2, addmargins)
# $DLPFC
#        GANDAL
# DLPFC   FALSE TRUE <NA>  Sum
#   FALSE  3449  661  966 5076
#   TRUE     52  284   70  406
#   <NA>    659  152  374 1185
#   Sum    4160 1097 1410 6667
#
# $HIPPO
#        GANDAL
# HIPPO   FALSE TRUE <NA>  Sum
#   FALSE  2479  437  791 3707
#   TRUE     40  177   53  270
#   <NA>   1641  483  566 2690
#   Sum    4160 1097 1410 6667

lapply(tabs2, addmargins)



## Adapted from https://github.com/LieberInstitute/brainseq_phase2/blob/master/twas/explore_twas_pgc2.R#L1655-L1699
twas_ov <- function(sig, prefix) {
    genes <- gandal$GeneID[p.adjust(gandal$TWAS.P, 'fdr') < 0.05]
    res <- map2_dfc(sig, names(sig), function(tsub, region) {
        res2 <- map_dfc(split(tsub, factor(tsub$feature, levels = features)), ~ genes %in% gsub('\\..*', '', .x$geneid))
        res2$any_feature <- pmap_lgl(res2, any)
        colnames(res2) <- paste0(region, '_', colnames(res2))
        return(res2)
    })
    colnames(res) <- paste0(prefix, colnames(res))
    return(res)
}

## Load the PGC2 TWAS results
get_pgc2 <- function() {
    e <- new.env()
    load('rda/pgc2_tt_objects.Rdata', verbose = TRUE, envir = e)
    e
}

pgc2 <- get_pgc2()

features <- c('gene', 'exon', 'jxn', 'tx')

table('Bonf' = gandal$TWAS.Bonferroni < 0.05, 'FDR' = p.adjust(gandal$TWAS.P, 'fdr') < 0.05)
#        FDR
# Bonf    FALSE  TRUE
#   FALSE 13054   912
#   TRUE      0   193

## Subset to only those with FDR < 5% (which includes the Bonf <5% ones)

gandal_gene <- cbind(
    'GeneID' = gandal$GeneID[p.adjust(gandal$TWAS.P, 'fdr') < 0.05],
    'bonferroni' = gandal$TWAS.Bonferroni[p.adjust(gandal$TWAS.P, 'fdr') < 0.05],
    twas_ov(pgc2$ttSig, 'pgc2_FDR_'),
    twas_ov(pgc2$ttSig_bonf, 'pgc2_Bonf_'),
    twas_ov(ttSig, 'psycm_FDR_'),
    twas_ov(ttSig_bonf, 'psycm_Bonf_')
)

## Is the gene in any of our results?
gandal_gene$in_any <- pmap_lgl(gandal_gene[, -(1:2)], any)
gandal_gene$in_any_FDR <- pmap_lgl(gandal_gene[, grep('FDR', colnames(gandal_gene))], any)
gandal_gene$in_any_Bonf <- pmap_lgl(gandal_gene[, grep('Bonf', colnames(gandal_gene))], any)

## Add gene info
gandal_gene <- cbind(gandal_gene, gandal[match(gandal_gene$GeneID, gandal$GeneID), c('gene_name', 'CHR', 'start', 'stop', 'strand', 'gene_type')])

dim(gandal_gene)
# [1] 1105   51

## Run a quick check
stopifnot(all(with(gandal_gene,
    pgc2_FDR_DLPFC_gene | pgc2_FDR_DLPFC_exon | pgc2_FDR_DLPFC_jxn | pgc2_FDR_DLPFC_tx == pgc2_FDR_DLPFC_any_feature
)))
save(gandal_gene, file = 'rda/gandal_gene.Rdata')


with(gandal_gene, addmargins(table(
    'CLOZUK+PGC2 DLPFC gene' = psycm_FDR_DLPFC_gene,
    'CLOZUK+PGC2 HIPPO gene' = psycm_FDR_HIPPO_gene
)))
#                       CLOZUK+PGC2 HIPPO gene
# CLOZUK+PGC2 DLPFC gene FALSE TRUE  Sum
#                  FALSE   842   40  882
#                  TRUE    125   98  223
#                  Sum     967  138 1105

with(subset(gandal_gene, bonferroni < 0.05), addmargins(table(
    'CLOZUK+PGC2 DLPFC gene' = psycm_Bonf_DLPFC_gene,
    'CLOZUK+PGC2 HIPPO gene' = psycm_Bonf_HIPPO_gene
)))
#                       CLOZUK+PGC2 HIPPO gene
# CLOZUK+PGC2 DLPFC gene FALSE TRUE Sum
#                  FALSE   154    6 160
#                  TRUE     21   12  33
#                  Sum     175   18 193

table(subset(gandal_gene, bonferroni < 0.05)$GeneID %in% gsub('\\..*', '', tt$geneid))
# FALSE  TRUE
#    75   118
table(subset(gandal_gene, bonferroni < 0.05)$GeneID %in% gsub('\\..*', '', tt$geneid[tt$feature == 'gene']))
# FALSE  TRUE
#   122    71

(21 + 12 + 6) / 71 * 100
# [1] 54.92958

(21 + 12 + 6) / 193 * 100
# [1] 20.20725

## Compare by PGC2 and CLOZUK+PGC2 TWAS by merging both brain regions
## TWAS FDR <5%
with(gandal_gene, addmargins(table(
    'PGC2' = pgc2_FDR_DLPFC_any_feature | pgc2_FDR_HIPPO_any_feature,
    'CLOZUK+PGC2' = psycm_FDR_DLPFC_any_feature | psycm_FDR_HIPPO_any_feature
)))
#        CLOZUK+PGC2
# PGC2    FALSE TRUE  Sum
#   FALSE   572  141  713
#   TRUE     11  381  392
#   Sum     583  522 1105

with(subset(gandal_gene, bonferroni < 0.05), addmargins(table(
    'PGC2' = pgc2_FDR_DLPFC_any_feature | pgc2_FDR_HIPPO_any_feature,
    'CLOZUK+PGC2' = psycm_FDR_DLPFC_any_feature | psycm_FDR_HIPPO_any_feature
)))
#        CLOZUK+PGC2
# PGC2    FALSE TRUE Sum
#   FALSE    84    9  93
#   TRUE      0  100 100
#   Sum      84  109 193

## TWAS Bonf <5%
with(subset(gandal_gene, bonferroni < 0.05), addmargins(table(
    'PGC2' = pgc2_Bonf_DLPFC_any_feature | pgc2_Bonf_HIPPO_any_feature,
    'CLOZUK+PGC2' = psycm_Bonf_DLPFC_any_feature | psycm_Bonf_HIPPO_any_feature
)))
#        CLOZUK+PGC2
# PGC2    FALSE TRUE Sum
#   FALSE   106   25 131
#   TRUE      1   61  62
#   Sum     107   86 193

with(subset(gandal_gene, bonferroni < 0.05), addmargins(table(
    'PGC2' = pgc2_Bonf_DLPFC_any_feature | pgc2_Bonf_HIPPO_any_feature,
    'CLOZUK+PGC2' = psycm_Bonf_DLPFC_any_feature | psycm_Bonf_HIPPO_any_feature
))) / sum(gandal_gene$bonferroni < 0.05) * 100

#        CLOZUK+PGC2
# PGC2          FALSE        TRUE         Sum
#   FALSE  54.9222798  12.9533679  67.8756477
#   TRUE    0.5181347  31.6062176  32.1243523
#   Sum    55.4404145  44.5595855 100.0000000

86 / 118 * 100
# [1] 72.88136
           
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
#  date     2019-04-01
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date       lib source
#  assertthat    0.2.1   2019-03-21 [2] CRAN (R 3.5.1)
#  cellranger    1.1.0   2016-07-27 [1] CRAN (R 3.5.0)
#  cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.1)
#  colorout    * 1.2-0   2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace    1.4-1   2019-03-18 [2] CRAN (R 3.5.1)
#  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)
#  digest        0.6.18  2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr         0.8.0.1 2019-02-15 [1] CRAN (R 3.5.1)
#  ggplot2     * 3.1.0   2018-10-25 [1] CRAN (R 3.5.1)
#  glue          1.3.1   2019-03-12 [1] CRAN (R 3.5.1)
#  gtable        0.3.0   2019-03-25 [2] CRAN (R 3.5.1)
#  htmltools     0.3.6   2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets   1.3     2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv        1.5.0   2019-03-15 [2] CRAN (R 3.5.1)
#  jsonlite      1.6     2018-12-07 [2] CRAN (R 3.5.1)
#  later         0.8.0   2019-02-11 [2] CRAN (R 3.5.1)
#  lattice       0.20-38 2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval      0.2.2   2019-03-15 [2] CRAN (R 3.5.1)
#  magrittr      1.5     2014-11-22 [1] CRAN (R 3.5.0)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 3.5.1)
#  patchwork   * 0.0.1   2019-04-01 [1] Github (thomasp85/patchwork@fd7958b)
#  pillar        1.3.1   2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig     2.0.2   2018-08-16 [1] CRAN (R 3.5.1)
#  plyr          1.8.4   2016-06-08 [2] CRAN (R 3.5.0)
#  png           0.1-7   2013-12-03 [2] CRAN (R 3.5.0)
#  promises      1.0.1   2018-04-13 [2] CRAN (R 3.5.0)
#  purrr       * 0.3.2   2019-03-15 [2] CRAN (R 3.5.1)
#  R6            2.4.0   2019-02-14 [2] CRAN (R 3.5.1)
#  Rcpp          1.0.0   2018-11-07 [1] CRAN (R 3.5.1)
#  readxl      * 1.3.1   2019-03-13 [2] CRAN (R 3.5.1)
#  rlang         0.3.1   2019-01-08 [1] CRAN (R 3.5.1)
#  rmote       * 0.3.4   2018-05-02 [1] deltarho (R 3.5.0)
#  scales        1.0.0   2018-08-09 [2] CRAN (R 3.5.1)
#  servr         0.13    2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo * 1.1.1   2018-11-05 [1] CRAN (R 3.5.1)
#  tibble        2.0.1   2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect    0.2.5   2018-10-11 [2] CRAN (R 3.5.1)
#  withr         2.1.2   2018-03-15 [2] CRAN (R 3.5.0)
#  xfun          0.5     2019-02-20 [1] CRAN (R 3.5.1)
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
