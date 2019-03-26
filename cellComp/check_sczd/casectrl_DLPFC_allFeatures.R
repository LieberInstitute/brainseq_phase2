# Adapted from casectrl_DLPFC.R and
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/caseControl_analysis_DLPFC.R
library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library('devtools')


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

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=FALSE)
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

# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.4.3 Patched (2018-01-20 r74142)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       US/Eastern
#  date     2018-04-26
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  base                 * 3.4.3     2018-01-20 local
#  Biobase              * 2.38.0    2017-11-07 Bioconductor
#  BiocGenerics         * 0.24.0    2017-11-29 Bioconductor
#  bitops                 1.0-6     2013-08-17 CRAN (R 3.4.1)
#  colorout             * 1.2-0     2018-02-19 Github (jalvesaq/colorout@2f01173)
#  colorspace             1.3-2     2016-12-14 CRAN (R 3.4.1)
#  compiler               3.4.3     2018-01-20 local
#  datasets             * 3.4.3     2018-01-20 local
#  DelayedArray         * 0.4.1     2017-11-07 Bioconductor
#  devtools             * 1.13.5    2018-02-18 CRAN (R 3.4.3)
#  digest                 0.6.15    2018-01-28 cran (@0.6.15)
#  edgeR                * 3.20.9    2018-04-18 Bioconductor
#  GenomeInfoDb         * 1.14.0    2017-11-29 Bioconductor
#  GenomeInfoDbData       1.0.0     2018-01-09 Bioconductor
#  GenomicRanges        * 1.30.3    2018-04-18 Bioconductor
#  ggplot2                2.2.1     2016-12-30 CRAN (R 3.4.1)
#  graphics             * 3.4.3     2018-01-20 local
#  grDevices            * 3.4.3     2018-01-20 local
#  grid                   3.4.3     2018-01-20 local
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.4.1)
#  htmltools              0.3.6     2017-04-28 CRAN (R 3.4.1)
#  htmlwidgets            1.2       2018-04-19 CRAN (R 3.4.3)
#  httpuv                 1.3.6.2   2018-03-02 CRAN (R 3.4.3)
#  IRanges              * 2.12.0    2017-11-29 Bioconductor
#  jaffelab             * 0.99.20   2018-04-19 Github (LieberInstitute/jaffelab@04c470a)
#  later                  0.7.1     2018-03-07 CRAN (R 3.4.3)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.4.3)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.4.2)
#  limma                * 3.34.9    2018-04-18 Bioconductor
#  locfit                 1.5-9.1   2013-04-20 CRAN (R 3.4.1)
#  Matrix                 1.2-12    2017-11-30 CRAN (R 3.4.3)
#  matrixStats          * 0.53.1    2018-02-11 CRAN (R 3.4.3)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.4.1)
#  methods              * 3.4.3     2018-01-20 local
#  mime                   0.5       2016-07-07 CRAN (R 3.4.1)
#  munsell                0.4.3     2016-02-13 CRAN (R 3.4.1)
#  parallel             * 3.4.3     2018-01-20 local
#  pillar                 1.2.1     2018-02-27 CRAN (R 3.4.3)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.4.1)
#  png                    0.1-7     2013-12-03 CRAN (R 3.4.1)
#  rafalib              * 1.0.0     2015-08-09 CRAN (R 3.4.1)
#  RColorBrewer           1.1-2     2014-12-07 CRAN (R 3.4.1)
#  Rcpp                   0.12.16   2018-03-13 CRAN (R 3.4.3)
#  RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.4.2)
#  rlang                  0.2.0     2018-02-20 CRAN (R 3.4.3)
#  rmote                * 0.3.4     2018-02-16 deltarho (R 3.4.3)
#  S4Vectors            * 0.16.0    2017-11-29 Bioconductor
#  scales                 0.5.0     2017-08-24 CRAN (R 3.4.1)
#  segmented              0.5-3.0   2017-11-30 CRAN (R 3.4.2)
#  servr                  0.9       2018-03-25 CRAN (R 3.4.3)
#  stats                * 3.4.3     2018-01-20 local
#  stats4               * 3.4.3     2018-01-20 local
#  SummarizedExperiment * 1.8.1     2018-01-09 Bioconductor
#  tibble                 1.4.2     2018-01-22 CRAN (R 3.4.3)
#  tools                  3.4.3     2018-01-20 local
#  utils                * 3.4.3     2018-01-20 local
#  withr                  2.1.2     2018-03-15 CRAN (R 3.4.3)
#  xfun                   0.1       2018-01-22 CRAN (R 3.4.3)
#  XVector                0.18.0    2017-11-29 Bioconductor
#  zlibbioc               1.24.0    2017-11-07 Bioconductor
#
