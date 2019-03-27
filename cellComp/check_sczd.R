library('sessioninfo')
library('purrr')
library('SummarizedExperiment')
library('broom')
library('jaffelab')
library('ggplot2')
library('tidyr')


files <- c('DLPFC' = '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_DLPFC.Rdata', 'HIPPO' = '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata')
file.exists(files)

load('methprop_pd.Rdata', verbose = TRUE)


sczd_lm <- map2_dfr(files, names(files), function(f, region) {
    load(f, verbose = TRUE)
    
    ## 
    # with qSVs:
    # modQsva <- modQsva
    # without qSVs:
    # modQsva <- mod
    # just dx:
    modQsva <- modQsva[, 1:2]
    
    
    res <- map_dfr(
        colnames(pd)[57:64],
        ~ cbind(tidy(lm(pd[keepIndex, .x] ~  modQsva - 1)), celltype = .x)
    )
    
    res$region <- region
    
    return(res)
})


dx_lm <- subset(sczd_lm, term == 'modQsvaDxSchizo')
dx_lm$p.bonf <- p.adjust(dx_lm$p.value, 'bonf')
dx_lm$p.fdr <- p.adjust(dx_lm$p.value, 'fdr')

options(width = 150)


## With qSVs:
dx_lm
#                term      estimate    std.error   statistic      p.value          celltype region      p.bonf        p.fdr
# 2   modQsvaDxSchizo -7.495378e-20 1.109544e-18 -0.06755366 0.9461794122 Fetal_replicating  DLPFC 1.000000000 0.9461794122
# 30  modQsvaDxSchizo  3.671069e-04 8.963282e-04  0.40956747 0.6823730944   Fetal_quiescent  DLPFC 1.000000000 0.7798549651
# 58  modQsvaDxSchizo -1.763644e-05 3.448260e-05 -0.51145894 0.6093513257               OPC  DLPFC 1.000000000 0.7499708624
# 86  modQsvaDxSchizo  4.364594e-03 2.810589e-03  1.55291077 0.1213456158           Neurons  DLPFC 1.000000000 0.3235883087
# 114 modQsvaDxSchizo  5.618424e-03 1.786842e-03  3.14433225 0.0018066469        Astrocytes  DLPFC 0.028906350 0.0072265876
# 142 modQsvaDxSchizo -5.402385e-03 4.676528e-03 -1.15521279 0.2487891032  Oligodendrocytes  DLPFC 1.000000000 0.4422917391
# 170 modQsvaDxSchizo -1.060408e-02 2.718298e-03 -3.90099859 0.0001148141         Microglia  DLPFC 0.001837025 0.0009185126
# 198 modQsvaDxSchizo  5.673974e-03 2.555703e-03  2.22012235 0.0270487446       Endothelial  DLPFC 0.432779914 0.0865559827
# 226 modQsvaDxSchizo -1.235050e-18 1.319839e-18 -0.93575813 0.3501404728 Fetal_replicating  HIPPO 1.000000000 0.5602247564
# 255 modQsvaDxSchizo  9.841033e-04 1.387236e-03  0.70939843 0.4786209780   Fetal_quiescent  HIPPO 1.000000000 0.6451113496
# 284 modQsvaDxSchizo -4.231540e-04 6.036354e-04 -0.70100918 0.4838335122               OPC  HIPPO 1.000000000 0.6451113496
# 313 modQsvaDxSchizo -3.398318e-03 2.644343e-03 -1.28512766 0.1997259139           Neurons  HIPPO 1.000000000 0.4087991349
# 342 modQsvaDxSchizo  9.362032e-03 2.331702e-03  4.01510647 0.0000749156        Astrocytes  HIPPO 0.001198650 0.0009185126
# 371 modQsvaDxSchizo  1.050929e-03 4.461831e-03  0.23553764 0.8139500575  Oligodendrocytes  HIPPO 1.000000000 0.8682133947
# 400 modQsvaDxSchizo -1.193057e-02 3.518362e-03 -3.39094564 0.0007885913         Microglia  HIPPO 0.012617461 0.0042058204
# 429 modQsvaDxSchizo  4.354981e-03 3.424139e-03  1.27184703 0.2043995675       Endothelial  HIPPO 1.000000000 0.4087991349

subset(dx_lm, p.bonf < 0.05)
#                term     estimate   std.error statistic      p.value   celltype region      p.bonf        p.fdr
# 114 modQsvaDxSchizo  0.005618424 0.001786842  3.144332 0.0018066469 Astrocytes  DLPFC 0.028906350 0.0072265876
# 170 modQsvaDxSchizo -0.010604078 0.002718298 -3.900999 0.0001148141  Microglia  DLPFC 0.001837025 0.0009185126
# 342 modQsvaDxSchizo  0.009362032 0.002331702  4.015106 0.0000749156 Astrocytes  HIPPO 0.001198650 0.0009185126
# 400 modQsvaDxSchizo -0.011930574 0.003518362 -3.390946 0.0007885913  Microglia  HIPPO 0.012617461 0.0042058204

## This is identical
subset(dx_lm, p.fdr < 0.05)


## Without qSVs:
dx_lm
#                term      estimate    std.error   statistic     p.value          celltype region     p.bonf      p.fdr
# 2   modQsvaDxSchizo  7.674786e-20 1.021891e-18  0.07510378 0.940173171 Fetal_replicating  DLPFC 1.00000000 0.94017317
# 15  modQsvaDxSchizo  1.042692e-04 8.362909e-04  0.12468057 0.900844890   Fetal_quiescent  DLPFC 1.00000000 0.94017317
# 28  modQsvaDxSchizo -1.621328e-05 3.449382e-05 -0.47003445 0.638610349               OPC  DLPFC 1.00000000 0.82506346
# 41  modQsvaDxSchizo -3.722487e-03 8.738396e-03 -0.42599203 0.670364058           Neurons  DLPFC 1.00000000 0.82506346
# 54  modQsvaDxSchizo  1.115355e-02 4.472193e-03  2.49397724 0.013072634        Astrocytes  DLPFC 0.20916214 0.06972071
# 67  modQsvaDxSchizo -1.343743e-02 8.195962e-03 -1.63951841 0.101964539  Oligodendrocytes  DLPFC 1.00000000 0.32628653
# 80  modQsvaDxSchizo -6.259620e-03 3.422548e-03 -1.82893588 0.068222771         Microglia  DLPFC 1.00000000 0.27289108
# 93  modQsvaDxSchizo  1.217794e-02 3.660912e-03  3.32647628 0.000968613       Endothelial  DLPFC 0.01549781 0.01549781
# 106 modQsvaDxSchizo -1.427219e-18 1.254737e-18 -1.13746478 0.256195014 Fetal_replicating  HIPPO 1.00000000 0.52484504
# 119 modQsvaDxSchizo  1.039186e-03 1.544726e-03  0.67273151 0.501603611   Fetal_quiescent  HIPPO 1.00000000 0.72960525
# 132 modQsvaDxSchizo  5.181142e-05 5.715428e-04  0.09065187 0.927825935               OPC  HIPPO 1.00000000 0.94017317
# 145 modQsvaDxSchizo -7.988527e-03 1.033172e-02 -0.77320432 0.439972240           Neurons  HIPPO 1.00000000 0.70395558
# 158 modQsvaDxSchizo  4.511292e-03 4.902272e-03  0.92024500 0.358137935        Astrocytes  HIPPO 1.00000000 0.63668966
# 171 modQsvaDxSchizo -1.177537e-02 8.187090e-03 -1.43828469 0.151330554  Oligodendrocytes  HIPPO 1.00000000 0.40354814
# 184 modQsvaDxSchizo -5.466201e-03 4.868962e-03 -1.12266239 0.262422519         Microglia  HIPPO 1.00000000 0.52484504
# 197 modQsvaDxSchizo  1.962781e-02 6.811114e-03  2.88173225 0.004222229       Endothelial  HIPPO 0.06755567 0.03377783
subset(dx_lm, p.bonf < 0.05)
#               term   estimate   std.error statistic     p.value    celltype region     p.bonf      p.fdr
# 93 modQsvaDxSchizo 0.01217794 0.003660912  3.326476 0.000968613 Endothelial  DLPFC 0.01549781 0.01549781
subset(dx_lm, p.fdr < 0.05)
#                term   estimate   std.error statistic     p.value    celltype region     p.bonf      p.fdr
# 93  modQsvaDxSchizo 0.01217794 0.003660912  3.326476 0.000968613 Endothelial  DLPFC 0.01549781 0.01549781
# 197 modQsvaDxSchizo 0.01962781 0.006811114  2.881732 0.004222229 Endothelial  HIPPO 0.06755567 0.03377783

## Just dx:
dx_lm
#               term      estimate    std.error    statistic      p.value          celltype region       p.bonf        p.fdr
# 2  modQsvaDxSchizo  8.460423e-21 9.584524e-19  0.008827171 9.929617e-01 Fetal_replicating  DLPFC 1.0000000000 0.9929616972
# 4  modQsvaDxSchizo -8.344877e-05 7.973923e-04 -0.104652081 9.167075e-01   Fetal_quiescent  DLPFC 1.0000000000 0.9825230827
# 6  modQsvaDxSchizo -3.609359e-06 3.642317e-05 -0.099095153 9.211154e-01               OPC  DLPFC 1.0000000000 0.9825230827
# 8  modQsvaDxSchizo -1.372423e-02 8.837534e-03 -1.552948380 1.212745e-01           Neurons  DLPFC 1.0000000000 0.2771988785
# 10 modQsvaDxSchizo  1.756957e-02 4.899042e-03  3.586328100 3.795030e-04        Astrocytes  DLPFC 0.0060720474 0.0020240158
# 12 modQsvaDxSchizo -1.118253e-02 7.903716e-03 -1.414844939 1.579394e-01  Oligodendrocytes  DLPFC 1.0000000000 0.3158787616
# 14 modQsvaDxSchizo -6.695610e-03 3.232257e-03 -2.071497030 3.899208e-02         Microglia  DLPFC 0.6238732297 0.1559683074
# 16 modQsvaDxSchizo  1.411986e-02 3.902716e-03  3.617957417 3.373895e-04       Endothelial  DLPFC 0.0053982313 0.0020240158
# 18 modQsvaDxSchizo -1.398931e-18 1.104318e-18 -1.266782755 2.061236e-01 Fetal_replicating  HIPPO 1.0000000000 0.3235014200
# 20 modQsvaDxSchizo -1.550131e-03 1.474590e-03 -1.051227996 2.939207e-01   Fetal_quiescent  HIPPO 1.0000000000 0.3918943260
# 22 modQsvaDxSchizo  6.222622e-04 5.090287e-04  1.222450054 2.224072e-01               OPC  HIPPO 1.0000000000 0.3235014200
# 24 modQsvaDxSchizo -1.603955e-02 9.939538e-03 -1.613712084 1.075428e-01           Neurons  HIPPO 1.0000000000 0.2771988785
# 26 modQsvaDxSchizo  4.435471e-03 4.569658e-03  0.970635199 3.324387e-01        Astrocytes  HIPPO 1.0000000000 0.4091553311
# 28 modQsvaDxSchizo -1.193916e-02 7.236380e-03 -1.649880355 9.991600e-02  Oligodendrocytes  HIPPO 1.0000000000 0.2771988785
# 30 modQsvaDxSchizo -5.717071e-03 4.618348e-03 -1.237903900 2.166290e-01         Microglia  HIPPO 1.0000000000 0.3235014200
# 32 modQsvaDxSchizo  3.018818e-02 7.053807e-03  4.279701025 2.453508e-05       Endothelial  HIPPO 0.0003925612 0.0003925612
subset(dx_lm, p.bonf < 0.05)
#               term   estimate   std.error statistic      p.value    celltype region       p.bonf        p.fdr
# 10 modQsvaDxSchizo 0.01756957 0.004899042  3.586328 3.795030e-04  Astrocytes  DLPFC 0.0060720474 0.0020240158
# 16 modQsvaDxSchizo 0.01411986 0.003902716  3.617957 3.373895e-04 Endothelial  DLPFC 0.0053982313 0.0020240158
# 32 modQsvaDxSchizo 0.03018818 0.007053807  4.279701 2.453508e-05 Endothelial  HIPPO 0.0003925612 0.0003925612
## Identical
# subset(dx_lm, p.fdr < 0.05)

sczd_cell <- map2_dfr(files, names(files), function(f, region) {
    load(f, verbose = TRUE)
    
    res <- map_dfr(
        colnames(pd)[57:64],
        function(cell) {
            data.frame(
                # cell = as.vector(cleaningY(matrix(pd[keepIndex, cell], nrow = 1), modQsva, P = 2)),
                cell = pd[keepIndex, cell],
                dx = pd$Dx[keepIndex],
                celltype = cell,
                stringsAsFactors = FALSE
            )
        }
    )
    
    res$region <- region
    
    return(res)
})

shorten_cell <- function(df, dx = TRUE) {
    if(dx) df$dx[df$dx == 'Schizo'] <- 'SCZD'
    df$celltype[df$celltype == 'Fetal_replicating'] <- 'FRN'
    df$celltype[df$celltype == 'Fetal_quiescent'] <- 'FQN'
    return(df)
}
sczd_cell <- shorten_cell(sczd_cell)

dx_lm$dx <- 'Schizo'
dx_lm <- shorten_cell(dx_lm)

pdf('sczd_cell.pdf', useDingbats = FALSE, width = 14, height = 25)
ggplot(sczd_cell, aes(x = dx, y = cell, fill = dx)) +
    geom_boxplot() +
    # geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    # geom_point(position = position_jitter(width = 0.2)) +
    facet_grid(celltype ~ region, scales = 'free_y') +
    # ylab('Cell RNA fraction - covariates removed') +
    ylab('Cell RNA fraction') +
    xlab('Diagnosis') +
    theme_bw(base_size = 30) +
    scale_fill_manual(values = c('Control' = 'orchid4', 'SCZD' = 'aquamarine4'), name = 'Diagnosis') +
    geom_text(
        data = dx_lm,
        color = 'black',
        size = 7,
        mapping = aes(x = 1.5, y = Inf, label = formatC(p.bonf, format = 'e', digits = 2), vjust = 1.5)
    ) +
    labs(caption = 'Bonferroni-adjusted p-value shown in each panel') +
    expand_limits(y = c(-0.0001, 0.0001))
dev.off()



sczd_qsv <- map2_dfr(files, names(files), function(f, region) {
    load(f, verbose = TRUE)
    
    stopifnot(all(sapply(pd$SAMPLE_ID[keepIndex], '[', 1) == rownames(qSVs)))
    
    res <- map_dfr(
        colnames(pd)[57:64],
        ~ cbind(tidy(lm(pd[keepIndex, .x] ~  qSVs)), celltype = .x)
    )
    res$region <- region
    
    return(res)
})
sczd_qsv$p.bonf <- p.adjust(sczd_qsv$p.value, 'bonf')
sczd_qsv$p.fdr <- p.adjust(sczd_qsv$p.value, 'fdr')
sczd_qsv <- shorten_cell(sczd_qsv, dx = FALSE)
dim(sczd_qsv)
# [1] 264   9
dim(subset(sczd_qsv, term != '(Intercept)'))
# [1] 248   9
head(sczd_qsv)
dim(subset(sczd_qsv, p.bonf < 0.05 & term != '(Intercept)'))
# [1] 98  9
dim(subset(sczd_qsv, p.fdr < 0.05 & term != '(Intercept)'))
# [1] 138   9
with(subset(sczd_qsv, p.bonf < 0.05 & term != '(Intercept)'), table(celltype, region))
#                   region
# celltype           DLPFC HIPPO
#   Astrocytes           9    11
#   Endothelial          9    10
#   FQN                  0     4
#   Microglia            8     6
#   Neurons              9    11
#   Oligodendrocytes     9     9
#   OPC                  2     1
with(subset(sczd_qsv, p.fdr < 0.05 & term != '(Intercept)'), table(celltype, region))
#                   region
# celltype           DLPFC HIPPO
#   Astrocytes          13    13
#   Endothelial         11    15
#   FQN                  2     6
#   FRN                  1     1
#   Microglia           11     9
#   Neurons             12    13
#   Oligodendrocytes    14    10
#   OPC                  6     1

sczd_qsv_tab <- map2_dfr(files, names(files), function(f, region) {
    load(f, verbose = TRUE)
    
    res <- map_dfr(
        colnames(pd)[57:64],
        function(cell) {
            cbind(data.frame(
                cell = pd[keepIndex, cell],
                dx = pd$Dx[keepIndex],
                celltype = cell,
                stringsAsFactors = FALSE
            ), qSVs)
        }
    )
    
    ## Make into long format
    res <- gather(res, qsv, qsv_value, grep('PC', colnames(res)))
    
    res$region <- region
    
    return(res)
})

## Simplify for plotting
sczd_qsv_tab <- shorten_cell(sczd_qsv_tab)
sczd_qsv_tab$qsv <- factor(sczd_qsv_tab$qsv, levels = unique(sczd_qsv_tab$qsv))

## Add info from the other table that I can then use for colors
m <- match(
    with(sczd_qsv_tab, paste('qSVs', qsv, '-', region, '-', celltype, sep = '')),
    with(sczd_qsv, paste(term, '-', region, '-', celltype, sep = ''))
)
sczd_qsv_tab <- cbind(sczd_qsv_tab, sczd_qsv[m, c('p.bonf', 'p.fdr')])

## Visualize all the qSVs versus the cell type proportions
pdf('sczd_cell_and_qsv.pdf', useDingbats = FALSE, width = 45, height = 55)
ggplot(sczd_qsv_tab, aes(x = qsv_value, y = cell, color = p.bonf < 0.05)) + 
    geom_point() +
    facet_grid(region * celltype ~ qsv, scales = 'free') +
    ylab('Cell RNA fraction') +
    xlab('Quality Surrogate Variable') +
    theme_bw(base_size = 30) +
    expand_limits(y = c(-0.01, 0.01)) +
    scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red'))
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
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.5.1)
#  backports              1.1.3     2018-12-14 [2] CRAN (R 3.5.1)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.6    2019-02-10 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  broom                * 0.5.1     2018-12-05 [1] CRAN (R 3.5.1)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                  0.8.0.1   2019-02-15 [1] CRAN (R 3.5.1)
#  generics               0.0.2     2018-11-29 [1] CRAN (R 3.5.1)
#  GenomeInfoDb         * 1.18.2    2019-02-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2              * 3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.5.1)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.5.0     2019-03-15 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab             * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  labeling               0.3       2014-08-23 [2] CRAN (R 3.5.0)
#  later                  0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.5.1)
#  limma                  3.38.3    2018-12-02 [1] Bioconductor
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  nlme                   3.1-137   2018-04-07 [3] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                * 0.3.2     2019-03-15 [2] CRAN (R 3.5.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.5.0)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr                  0.13      2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  stringi                1.4.3     2019-03-12 [2] CRAN (R 3.5.1)
#  stringr                1.4.0     2019-02-10 [1] CRAN (R 3.5.1)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyr                * 0.8.3     2019-03-01 [2] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.5       2019-02-20 [1] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
