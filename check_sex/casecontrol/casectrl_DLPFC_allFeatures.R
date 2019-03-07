# Adapted from /dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/casectrl_DLPFC_allFeatures.R 
library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library('devtools')

dir.create('rdas', showWarnings = FALSE)
dir.create('pdf', showWarnings = FALSE)

### Run with Gold samples using the qSVs made without the age>17 samples
## and without the HIPPO Gold samples (qsv are DLPFC specific)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata", verbose = TRUE)
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_DLPFC.Rdata', verbose = TRUE)

## Drop samples absent in mod and modQsVA
rse_gene <- rse_gene[, keepIndex]
rse_jxn <- rse_jxn[, keepIndex]
rse_exon <- rse_exon[, keepIndex]
rse_tx <- rse_tx[, keepIndex]

## Keep region-specific samples
keepIndex = which(rse_gene$Age>17 &
			rse_gene$Region == "DLPFC")
rse_gene <- rse_gene[, keepIndex]
rse_jxn <- rse_jxn[, keepIndex]
rse_exon <- rse_exon[, keepIndex]
rse_tx <- rse_tx[, keepIndex]
mod <- mod[keepIndex, ]
modQsva <- modQsva[keepIndex, ]

## Make sex the second coefficient
check_sex <- function(mat) {
    i_sex <- grep('sex', tolower(colnames(mat)))
    i_rest <- seq_len(ncol(mat))
    i_rest <- i_rest[-c(1, i_sex)]
    mat[, c(1, i_sex, i_rest)]
}
mod <- check_sex(mod)
modQsva <- check_sex(modQsva)
colnames(mod)
#  [1] "(Intercept)"       "SexM"              "DxSchizo"
#  [4] "Age"               "mitoRate"          "rRNA_rate"
#  [7] "totalAssignedGene" "RIN"               "snpPC1"
# [10] "snpPC2"            "snpPC3"            "snpPC4"
# [13] "snpPC5"
colnames(modQsva)
#  [1] "(Intercept)"       "SexM"              "DxSchizo"
#  [4] "Age"               "mitoRate"          "rRNA_rate"
#  [7] "totalAssignedGene" "RIN"               "snpPC1"
# [10] "snpPC2"            "snpPC3"            "snpPC4"
# [13] "snpPC5"            "PC1"               "PC2"
# [16] "PC3"               "PC4"               "PC5"
# [19] "PC6"               "PC7"               "PC8"
# [22] "PC9"               "PC10"              "PC11"
# [25] "PC12"              "PC13"              "PC14"
# [28] "PC15"

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)
pdf('pdf/dlpfc_voom_qsva_noHGoldQSV_matchDLPFC_gene.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]


##### Exon ######
dee = DGEList(counts = assays(rse_exon)$counts,
	genes = rowData(rse_exon))
dee = calcNormFactors(dee)
pdf('pdf/dlpfc_voom_qsva_noHGoldQSV_matchDLPFC_exon.pdf', useDingbats = FALSE)
vExon = voom(dee,modQsva, plot=TRUE)
dev.off()
fitExon = lmFit(vExon)
eBExon = eBayes(fitExon)
sigExon = topTable(eBExon,coef=2,
	p.value = 1,number=nrow(rse_exon))
outExon = sigExon[rownames(rse_exon),]

##### Junction ######
dje = DGEList(counts = assays(rse_jxn)$counts,
	genes = rowData(rse_jxn))
dje = calcNormFactors(dje)
pdf('pdf/dlpfc_voom_qsva_noHGoldQSV_matchDLPFC_jxn.pdf', useDingbats = FALSE)
vJxn = voom(dje,modQsva, plot=TRUE)
dev.off()
fitJxn = lmFit(vJxn)
eBJxn = eBayes(fitJxn)
sigJxn = topTable(eBJxn,coef=2,
	p.value = 1,number=nrow(rse_jxn))
outJxn = sigJxn[rownames(rse_jxn),]

##### Transcript ######
fitTx = lmFit(log2(assays(rse_tx)$tpm + 0.5), modQsva)
eBTx = eBayes(fitTx)
sigTx = topTable(eBTx,coef=2,
	p.value = 1,number=nrow(rse_tx))
outTx = sigTx[rownames(rse_tx),]
outTx <- cbind(outTx, rowData(rse_tx))

save(outGene, outExon, outJxn,outTx,
	file = "rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda")


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
#  backports              1.1.3     2018-12-14 [2] CRAN (R 3.5.1)
#  bindr                  0.1.1     2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 [1] CRAN (R 3.5.0)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.5    2019-01-04 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  callr                  3.1.1     2018-12-21 [2] CRAN (R 3.5.1)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  desc                   1.2.0     2018-05-01 [2] CRAN (R 3.5.1)
#  devtools             * 2.0.1     2018-10-26 [1] CRAN (R 3.5.1)
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                  0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
#  edgeR                * 3.24.3    2019-01-02 [1] Bioconductor
#  fs                     1.2.6     2018-08-23 [2] CRAN (R 3.5.1)
#  GenomeInfoDb         * 1.18.1    2018-11-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2                3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.0     2018-07-17 [1] CRAN (R 3.5.1)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab             * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  later                  0.7.5     2018-09-18 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma                * 3.38.3    2018-12-02 [1] Bioconductor
#  locfit                 1.5-9.1   2013-04-20 [2] CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  memoise                1.1.0     2017-04-21 [2] CRAN (R 3.5.0)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.0)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgbuild               1.0.2     2018-10-16 [2] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  pkgload                1.0.2     2018-10-29 [2] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  prettyunits            1.0.2     2015-07-13 [1] CRAN (R 3.5.0)
#  processx               3.2.1     2018-12-05 [1] CRAN (R 3.5.1)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  ps                     1.3.0     2018-12-21 [2] CRAN (R 3.5.1)
#  purrr                  0.2.5     2018-05-29 [2] CRAN (R 3.5.0)
#  R6                     2.3.0     2018-10-04 [2] CRAN (R 3.5.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  remotes                2.0.2     2018-10-30 [1] CRAN (R 3.5.1)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  rprojroot              1.3-2     2018-01-03 [2] CRAN (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr                  0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo            1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  testthat               2.0.1     2018-10-13 [1] CRAN (R 3.5.1)
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  usethis              * 1.4.0     2018-08-14 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
