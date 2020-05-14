library("SummarizedExperiment")
library("sessioninfo")

## Load the DE results
load(here::here("supp_tabs", "deres.Rdata"), verbose = TRUE)

## How big are they?
pryr::object_size(deres)
# 1.59 GB

lapply(deres, names)
# $development
# [1] "gene" "exon" "jxn"  "tx"
#
# $region
# [1] "gene" "exon" "jxn"  "tx"
#
# $sczd
# [1] "gene" "exon" "jxn"  "tx"

## Change commas to semi-colons before writing as a csv file
fix_csv <- function(df) {
    for (i in seq_len(ncol(df))) {
        if (any(grepl(",", df[, i]))) {
            message(paste(Sys.time(), "fixing column", colnames(df)[i]))
            df[, i] <- gsub(",", ";", df[, i])
        }
    }
    return(df)
}

for (i in 1:3) {
    message(paste(Sys.time(), "processing", names(deres)[i]))
    deres[[i]] <- lapply(deres[[i]], fix_csv)
}
# 2020-05-04 21:52:18 processing development
# 2020-05-04 21:52:21 fixing column gencodeTx
# 2020-05-04 21:53:00 fixing column gencodeTx
# 2020-05-04 21:53:40 fixing column gencodeTx
# 2020-05-04 21:54:02 processing region
# 2020-05-04 21:54:07 fixing column gencodeTx
# 2020-05-04 21:55:05 fixing column gencodeTx
# 2020-05-04 21:56:02 fixing column gencodeTx
# 2020-05-04 21:56:35 processing sczd
# 2020-05-04 21:56:38 fixing column gencodeTx
# 2020-05-04 21:57:04 fixing column gencodeTx
# 2020-05-04 21:57:49 fixing column gencodeTx

## Load the exon name map
exon_name_map <- read.delim("rda/BrainSeqPhaseII_exon_name_map.txt", row.names = 1)
head(exon_name_map)
#   libd_bsp2           gencode libd_gtex
# 1       e10 ENSE00001890219.1       e10
# 2       e11 ENSE00003507205.1       e11
# 3       e12 ENSE00003477500.1       e12
# 4       e13 ENSE00003565697.1       e13
# 5       e14 ENSE00003475637.1       e14
# 6       e15 ENSE00003502542.1       e15

## Fix the sczd rownames
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda", verbose = TRUE)
outFeat_hippo <- list(
    "gene" = outGene,
    "exon" = outExon,
    "jxn" = outJxn,
    "tx" = outTx
)

stopifnot(identical(deres$sczd$exon$t[-seq_len(nrow(outExon))], outExon$t))

for (i in names(outFeat_hippo)) {
    ## DLPFC first, then HIPPO based on https://github.com/LieberInstitute/brainseq_phase2/blob/master/supp_tabs/create_supp_tables.R#L128
    rownames(deres[[3]][[i]]) <-
        paste0(
            rep(paste0(c("DLPFC_", "HIPPO_"), i), each = nrow(outFeat_hippo[[i]])),
            ".",
            rep(rownames(outFeat_hippo[[i]]), 2)
        )
}

## add the feature_id column
add_feature_id <- function(x) {
    is_exon <- grepl("exon", rownames(x)[1])
    x$feature_id <- gsub(".*[gene|exon|jxn|tx]\\.", "", rownames(x))
    if (is_exon) {
        x$feature_id <- exon_name_map$gencode[match(x$feature_id, exon_name_map$libd_bsp2)]
    }
    return(x)
}
for (i in 1:3) {
    message(paste(Sys.time(), "processing", names(deres)[i]))
    deres[[i]] <- lapply(deres[[i]], add_feature_id)
}
# 2020-05-14 10:48:56 processing development
# 2020-05-14 10:48:59 processing regionspecific
# 2020-05-14 10:49:06 processing sczd_casecontrol

## Match the rest of the names used for the other browser files
names(deres) <- c("development", "regionspecific", "sczd_casecontrol")

## Check the exon case
lapply(deres, function(x) table(is.na(x$exon$feature_id)))
## Check all
sapply(deres, function(x) {
    sapply(x, function(y) {
        any(is.na(y$feature_id))
    })
})



## Explore the columns for each model
lapply(deres, function(x) head(x$gene, n = 2))
# $development
#                        Age.RegionHIPPO RegionHIPPO.fetal RegionHIPPO.birth
# gene.ENSG00000227232.5       0.2697299         0.1558122         0.1322596
# gene.ENSG00000278267.1       1.9601751         1.2619333        -1.4895673
#                        RegionHIPPO.infant RegionHIPPO.child RegionHIPPO.teen
# gene.ENSG00000227232.5         -0.4566381        0.05835068      0.003025649
# gene.ENSG00000278267.1         -0.4741538       -0.04049151      0.047792118
#                        RegionHIPPO.adult   AveExpr        F       P.Value
# gene.ENSG00000227232.5      -0.016199796  1.028798 376.4546 1.292818e-213
# gene.ENSG00000278267.1      -0.003501927 -1.989715 131.4779 2.722490e-116
#                            adj.P.Val type        P.Bonf span_Age.RegionHIPPO
# gene.ENSG00000227232.5 5.855326e-213 gene 3.187054e-209            8.0021738
# gene.ENSG00000278267.1 8.087096e-116 gene 6.711481e-112           -0.3196381
#                        span_RegionHIPPO.fetal span_RegionHIPPO.birth
# gene.ENSG00000227232.5              4.8358656             -7.6366303
# gene.ENSG00000278267.1             -0.9955068             -0.9892802
#                        span_RegionHIPPO.infant span_RegionHIPPO.child
# gene.ENSG00000227232.5              -0.3427108              0.3053811
# gene.ENSG00000278267.1               1.5505467             -0.3984407
#                        span_RegionHIPPO.teen span_AveExpr   span_F span_P.Value
# gene.ENSG00000227232.5            -0.4897329    -3.059094 4.880797 3.885935e-04
# gene.ENSG00000278267.1             0.2114331    -4.532517 7.143530 8.320897e-06
#                        span_adj.P.Val span_type span_P.Bonf seqnames start
# gene.ENSG00000227232.5   1.397463e-03      gene   1.0000000     chr1 14404
# gene.ENSG00000278267.1   4.274365e-05      gene   0.2051268     chr1 17369
#                          end width strand Length         gencodeID
# gene.ENSG00000227232.5 29570 15167      -   1351 ENSG00000227232.5
# gene.ENSG00000278267.1 17436    68      -     68 ENSG00000278267.1
#                              ensemblID              gene_type    Symbol
# gene.ENSG00000227232.5 ENSG00000227232 unprocessed_pseudogene    WASH7P
# gene.ENSG00000278267.1 ENSG00000278267                  miRNA MIR6859-1
#                         EntrezID Class meanExprs NumTx         gencodeTx
# gene.ENSG00000227232.5        NA InGen  1.696980     1 ENST00000488147.1
# gene.ENSG00000278267.1 102466751 InGen  4.355235     1 ENST00000619216.1
#                        passExprsCut replicates_in_BrainSpan        feature_id
# gene.ENSG00000227232.5         TRUE                    TRUE ENSG00000227232.5
# gene.ENSG00000278267.1         TRUE                    TRUE ENSG00000278267.1
#
# $regionspecific
#                                    logFC    AveExpr          t     P.Value
# adult_gene.ENSG00000227232.5 -0.21909569  0.9292268 -2.7371981 0.006440476
# adult_gene.ENSG00000278267.1 -0.02520089 -2.0587727 -0.1611352 0.872058842
#                               adj.P.Val         B   age type P.Bonf span_logFC
# adult_gene.ENSG00000227232.5 0.01122054 -3.890969 adult gene      1  0.4870191
# adult_gene.ENSG00000278267.1 0.89795725 -7.093203 adult gene      1 -0.1678198
#                              span_AveExpr     span_t span_P.Value
# adult_gene.ENSG00000227232.5    -2.872806  0.3855199    0.7102749
# adult_gene.ENSG00000278267.1    -4.524791 -0.3009702    0.7713967
#                              span_adj.P.Val    span_B span_age span_type
# adult_gene.ENSG00000227232.5      0.9741138 -5.688121    adult      gene
# adult_gene.ENSG00000278267.1      0.9955806 -5.704816    adult      gene
#                              span_P.Bonf seqnames start   end width strand
# adult_gene.ENSG00000227232.5           1     chr1 14404 29570 15167      -
# adult_gene.ENSG00000278267.1           1     chr1 17369 17436    68      -
#                              Length         gencodeID       ensemblID
# adult_gene.ENSG00000227232.5   1351 ENSG00000227232.5 ENSG00000227232
# adult_gene.ENSG00000278267.1     68 ENSG00000278267.1 ENSG00000278267
#                                           gene_type    Symbol  EntrezID Class
# adult_gene.ENSG00000227232.5 unprocessed_pseudogene    WASH7P        NA InGen
# adult_gene.ENSG00000278267.1                  miRNA MIR6859-1 102466751 InGen
#                              meanExprs NumTx         gencodeTx passExprsCut
# adult_gene.ENSG00000227232.5  1.696980     1 ENST00000488147.1         TRUE
# adult_gene.ENSG00000278267.1  4.355235     1 ENST00000619216.1         TRUE
#                              replicates_in_BrainSpan        feature_id
# adult_gene.ENSG00000227232.5                   FALSE ENSG00000227232.5
# adult_gene.ENSG00000278267.1                   FALSE ENSG00000278267.1
#
# $sczd_casecontrol
#                              Length         gencodeID       ensemblID
# DLPFC_gene.ENSG00000227232.5   1351 ENSG00000227232.5 ENSG00000227232
# DLPFC_gene.ENSG00000278267.1     68 ENSG00000278267.1 ENSG00000278267
#                                           gene_type    Symbol  EntrezID Class
# DLPFC_gene.ENSG00000227232.5 unprocessed_pseudogene    WASH7P        NA InGen
# DLPFC_gene.ENSG00000278267.1                  miRNA MIR6859-1 102466751 InGen
#                              meanExprs NumTx         gencodeTx passExprsCut
# DLPFC_gene.ENSG00000227232.5  1.696980     1 ENST00000488147.1         TRUE
# DLPFC_gene.ENSG00000278267.1  4.355235     1 ENST00000619216.1         TRUE
#                                     logFC   AveExpr           t   P.Value
# DLPFC_gene.ENSG00000227232.5 -0.021942394  1.345579 -0.38564882 0.6999870
# DLPFC_gene.ENSG00000278267.1  0.007710115 -1.679423  0.08617349 0.9313769
#                              adj.P.Val         B region type        feature_id
# DLPFC_gene.ENSG00000227232.5 0.9100562 -5.820263  DLPFC gene ENSG00000227232.5
# DLPFC_gene.ENSG00000278267.1 0.9819371 -5.287980  DLPFC gene ENSG00000278267.1


## Export to csv files
for (i in 1:3) {
    for (j in 1:4) {
        message(paste(
            Sys.time(),
            "processing",
            names(deres)[i],
            "at the",
            names(deres[[i]])[j],
            "level"
        ))
        write.csv(deres[[i]][[j]],
            file = paste0(
                "BrainSeqPhaseII_stats_",
                names(deres)[i],
                "_",
                names(deres[[i]])[j],
                ".csv"
            )
        )
    }
}
# 2020-05-14 10:53:13 processing development at the gene level
# 2020-05-14 10:53:15 processing development at the exon level
# 2020-05-14 10:53:42 processing development at the jxn level
# 2020-05-14 10:54:04 processing development at the tx level
# 2020-05-14 10:54:11 processing regionspecific at the gene level
# 2020-05-14 10:54:13 processing regionspecific at the exon level
# 2020-05-14 10:54:51 processing regionspecific at the jxn level
# 2020-05-14 10:55:23 processing regionspecific at the tx level
# 2020-05-14 10:55:33 processing sczd_casecontrol at the gene level
# 2020-05-14 10:55:34 processing sczd_casecontrol at the exon level
# 2020-05-14 10:55:54 processing sczd_casecontrol at the jxn level
# 2020-05-14 10:56:13 processing sczd_casecontrol at the tx level

system("chmod 770 BrainSeqPhaseII_stats_*")
system("ls -lh BrainSeqPhaseII_stats_*")
# -rwxrwx--- 1 lcollado lieber_jaffe 247M May 14 10:53 BrainSeqPhaseII_stats_development_exon.csv
# -rwxrwx--- 1 lcollado lieber_jaffe  18M May 14 10:53 BrainSeqPhaseII_stats_development_gene.csv
# -rwxrwx--- 1 lcollado lieber_jaffe 206M May 14 10:54 BrainSeqPhaseII_stats_development_jxn.csv
# -rwxrwx--- 1 lcollado lieber_jaffe  65M May 14 10:54 BrainSeqPhaseII_stats_development_tx.csv
# -rwxrwx--- 1 lcollado lieber_jaffe 378M May 14 10:54 BrainSeqPhaseII_stats_regionspecific_exon.csv
# -rwxrwx--- 1 lcollado lieber_jaffe  28M May 14 10:54 BrainSeqPhaseII_stats_regionspecific_gene.csv
# -rwxrwx--- 1 lcollado lieber_jaffe 320M May 14 10:55 BrainSeqPhaseII_stats_regionspecific_jxn.csv
# -rwxrwx--- 1 lcollado lieber_jaffe 105M May 14 10:55 BrainSeqPhaseII_stats_regionspecific_tx.csv
# -rwxrwx--- 1 lcollado lieber_jaffe 246M May 14 10:55 BrainSeqPhaseII_stats_sczd_casecontrol_exon.csv
# -rwxrwx--- 1 lcollado lieber_jaffe  20M May 14 10:55 BrainSeqPhaseII_stats_sczd_casecontrol_gene.csv
# -rwxrwx--- 1 lcollado lieber_jaffe 223M May 14 10:56 BrainSeqPhaseII_stats_sczd_casecontrol_jxn.csv
# -rwxrwx--- 1 lcollado lieber_jaffe  76M May 14 10:56 BrainSeqPhaseII_stats_sczd_casecontrol_tx.csv
system("wc -l BrainSeqPhaseII_stats_*")
#  396956 BrainSeqPhaseII_stats_development_exon.csv
#   25240 BrainSeqPhaseII_stats_development_gene.csv
#  297528 BrainSeqPhaseII_stats_development_jxn.csv
#   92733 BrainSeqPhaseII_stats_development_tx.csv
#  793911 BrainSeqPhaseII_stats_regionspecific_exon.csv
#   50479 BrainSeqPhaseII_stats_regionspecific_gene.csv
#  595055 BrainSeqPhaseII_stats_regionspecific_jxn.csv
#  185465 BrainSeqPhaseII_stats_regionspecific_tx.csv
#  793919 BrainSeqPhaseII_stats_sczd_casecontrol_exon.csv
#   50479 BrainSeqPhaseII_stats_sczd_casecontrol_gene.csv
#  595055 BrainSeqPhaseII_stats_sczd_casecontrol_jxn.csv
#  185465 BrainSeqPhaseII_stats_sczd_casecontrol_tx.csv
# 4062285 total

styler::style_file("export_de_stats.R", transformers = biocthis::bioc_style())

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.0.0 Patched (2020-05-13 r78451)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2020-05-14
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version  date       lib source
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.0)
#  backports              1.1.6    2020-04-05 [1] CRAN (R 4.0.0)
#  Biobase              * 2.48.0   2020-04-27 [1] Bioconductor
#  BiocGenerics         * 0.34.0   2020-04-27 [1] Bioconductor
#  bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.0)
#  cli                    2.0.2    2020-02-28 [1] CRAN (R 4.0.0)
#  codetools              0.2-16   2018-12-24 [3] CRAN (R 4.0.0)
#  colorout             * 1.2-2    2020-05-08 [1] Github (jalvesaq/colorout@726d681)
#  colorspace             1.4-1    2019-03-18 [2] CRAN (R 4.0.0)
#  crayon                 1.3.4    2017-09-16 [1] CRAN (R 4.0.0)
#  DelayedArray         * 0.14.0   2020-04-27 [1] Bioconductor
#  digest                 0.6.25   2020-02-23 [1] CRAN (R 4.0.0)
#  dplyr                  0.8.5    2020-03-07 [1] CRAN (R 4.0.0)
#  ellipsis               0.3.0    2019-09-20 [1] CRAN (R 4.0.0)
#  fansi                  0.4.1    2020-01-08 [1] CRAN (R 4.0.0)
#  GenomeInfoDb         * 1.24.0   2020-04-27 [1] Bioconductor
#  GenomeInfoDbData       1.2.3    2020-04-21 [2] Bioconductor
#  GenomicRanges        * 1.40.0   2020-04-27 [1] Bioconductor
#  ggplot2                3.3.0    2020-03-05 [1] CRAN (R 4.0.0)
#  glue                   1.4.0    2020-04-03 [1] CRAN (R 4.0.0)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.0)
#  here                   0.1      2017-05-28 [1] CRAN (R 4.0.0)
#  htmltools              0.4.0    2019-10-04 [1] CRAN (R 4.0.0)
#  htmlwidgets            1.5.1    2019-10-08 [1] CRAN (R 4.0.0)
#  httpuv                 1.5.2    2019-09-11 [1] CRAN (R 4.0.0)
#  IRanges              * 2.22.1   2020-04-28 [1] Bioconductor
#  jsonlite               1.6.1    2020-02-02 [2] CRAN (R 4.0.0)
#  later                  1.0.0    2019-10-04 [1] CRAN (R 4.0.0)
#  lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.0)
#  lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.0)
#  magrittr               1.5      2014-11-22 [1] CRAN (R 4.0.0)
#  Matrix                 1.2-18   2019-11-27 [3] CRAN (R 4.0.0)
#  matrixStats          * 0.56.0   2020-03-13 [1] CRAN (R 4.0.0)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.0)
#  pillar                 1.4.4    2020-05-05 [1] CRAN (R 4.0.0)
#  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.0.0)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 4.0.0)
#  promises               1.1.0    2019-10-04 [1] CRAN (R 4.0.0)
#  pryr                   0.1.4    2018-02-18 [2] CRAN (R 4.0.0)
#  purrr                  0.3.4    2020-04-17 [1] CRAN (R 4.0.0)
#  R6                     2.4.1    2019-11-12 [2] CRAN (R 4.0.0)
#  Rcpp                   1.0.4.6  2020-04-09 [1] CRAN (R 4.0.0)
#  RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.0)
#  rlang                  0.4.6    2020-05-02 [1] CRAN (R 4.0.0)
#  rmote                * 0.3.4    2020-05-08 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              1.3-2    2018-01-03 [2] CRAN (R 4.0.0)
#  S4Vectors            * 0.26.0   2020-04-27 [1] Bioconductor
#  scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.0)
#  servr                  0.16     2020-03-02 [1] CRAN (R 4.0.0)
#  sessioninfo          * 1.1.1    2018-11-05 [1] CRAN (R 4.0.0)
#  stringi                1.4.6    2020-02-17 [2] CRAN (R 4.0.0)
#  stringr                1.4.0    2019-02-10 [1] CRAN (R 4.0.0)
#  SummarizedExperiment * 1.18.1   2020-04-30 [1] Bioconductor
#  tibble                 3.0.1    2020-04-20 [1] CRAN (R 4.0.0)
#  tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.0)
#  vctrs                  0.2.4    2020-03-10 [1] CRAN (R 4.0.0)
#  withr                  2.2.0    2020-04-20 [1] CRAN (R 4.0.0)
#  xfun                   0.13     2020-04-13 [1] CRAN (R 4.0.0)
#  XVector                0.28.0   2020-04-27 [1] Bioconductor
#  zlibbioc               1.34.0   2020-04-27 [1] Bioconductor
#
# [1] /users/lcollado/R/4.0
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/library
