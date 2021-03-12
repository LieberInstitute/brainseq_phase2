####

### libraries
library("SummarizedExperiment")
library("jaffelab")
library("MatrixEQTL")
library("sva")
library("here")
library("sessioninfo")

######################
### load data ####
######################

load(here::here("expr_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here::here("expr_cutoff", "rse_exon.Rdata"), verbose = TRUE)
load(here::here("expr_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
load(here::here("expr_cutoff", "rse_tx.Rdata"), verbose = TRUE)

## keep adult samples & correct region
keepInd = which(colData(rse_gene)$Age > 13 & colData(rse_gene)$Region == "DLPFC")
rse_gene = rse_gene[,keepInd]
rse_exon = rse_exon[,keepInd]
rse_jxn = rse_jxn[,keepInd]
rse_tx = rse_tx[,keepInd]

## extract pd and rpkms
pd = colData(rse_gene)
geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm


######################
### snp data ####
######################

## load SNP data
load(
    here::here(
        "eQTL_full",
        "subset_hancock",
        "BrainSeq_Phase2_RiboZero_Genotypes_n551_subset_hancock.Rdata"
    ),
    verbose = TRUE
)

### make mds and snp dimensions equal to N
###(repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum


## drop SNPs not mapping to hg38
table(is.na(snpMap$chr_hg38))
# FALSE
#    51
   
# keepIndex = which(!is.na(snpMap$chr_hg38))
# snpMap = snpMap[keepIndex,]
# snp = snp[keepIndex,]


######################
# statistical model ##
######################

pd$Dx = factor(pd$Dx,
	levels = c("Control", "Schizo"))

mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]),
	data = pd)
colnames(mod)[4:8] = colnames(mds)[1:5]


######################
# create SNP objects #
######################

theSnps = SlicedData$new(as.matrix(snp))
theSnps$ResliceCombined(sliceSize = 50000)

snpspos = snpMap[,c("SNP","chr_hg38","pos_hg38")]
colnames(snpspos) = c("name","chr","pos")


#######################
####### do PCA ########
#######################

load(
    here::here(
        "eQTL_full",
        "eqtl_tables",
        "rdas",
        "pcs_dlpfc_4features_filtered_over13.rda"
    ),
    verbose = TRUE
)

colnames(mod)[-1]
# [1] "DxSchizo" "SexM"     "snpPC1"   "snpPC2"   "snpPC3"   "snpPC4"   "snpPC5"

covsGene = SlicedData$new(t(cbind(mod[,-1],genePCs)))
covsExon = SlicedData$new(t(cbind(mod[,-1],exonPCs)))
covsJxn = SlicedData$new(t(cbind(mod[,-1],jxnPCs)))
covsTx = SlicedData$new(t(cbind(mod[,-1],txPCs)))

##########################
### feature annotation ###
##########################

###### gene level
posGene = as.data.frame(rowRanges(rse_gene))[,1:3]
posGene$name = rownames(posGene)
posGene = posGene[,c(4,1:3)]

##### exon level 
posExon = as.data.frame(rowRanges(rse_exon))[,1:3]
posExon$name = rownames(posExon)
posExon = posExon[,c(4,1:3)]

##### junction level 
posJxn = as.data.frame(rowRanges(rse_jxn))[,1:3]
posJxn$name = rownames(posJxn)
posJxn = posJxn[,c(4,1:3)]
names(posJxn)[2:4] = c("Chr", "Start","End")

##### transcript level 
posTx = as.data.frame(rowRanges(rse_tx))[,1:3]
posTx$name = rownames(posTx)
posTx = posTx[,c(4,1:3)]
names(posTx)[2:4] = c("Chr", "Start","End")


#############################
### sliced expression data ##
geneSlice = SlicedData$new(log2(geneRpkm+1))
exonSlice = SlicedData$new(log2(exonRpkm+1))
jxnSlice = SlicedData$new(log2(jxnRp10m+1))
txSlice = SlicedData$new(log2(txTpm+1))

geneSlice$ResliceCombined(sliceSize = 5000)
exonSlice$ResliceCombined(sliceSize = 5000)
jxnSlice$ResliceCombined(sliceSize = 5000)
txSlice$ResliceCombined(sliceSize = 5000)


##########################
### Run EQTLs ############
##########################
print("Begin eQTL analysis")

meGene = Matrix_eQTL_main(snps=theSnps, gene = geneSlice, 
	cvrt = covsGene, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posGene, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	

meExon = Matrix_eQTL_main(snps=theSnps, gene = exonSlice, 
	cvrt = covsExon, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posExon, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	

meJxn = Matrix_eQTL_main(snps=theSnps, gene = jxnSlice, 
	cvrt = covsJxn, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posJxn, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	
	
meTx = Matrix_eQTL_main(snps=theSnps, gene = txSlice, 
	cvrt = covsTx, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posTx, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	

save(meGene, meExon, meJxn, meTx,
	file=here::here(
        "eQTL_full",
        "subset_hancock",
        "matrixEqtl_output_dlpfc_4features_subset_hancock.Rdata")
)

	
######################
###### annotate ######

# extract
geneEqtl = meGene$cis$eqtls
geneEqtl$gene = as.character(geneEqtl$gene)
geneEqtl$snps = as.character(geneEqtl$snps)

exonEqtl = meExon$cis$eqtls
exonEqtl$gene = as.character(exonEqtl$gene)
exonEqtl$snps = as.character(exonEqtl$snps)

jxnEqtl = meJxn$cis$eqtls
jxnEqtl$gene = as.character(jxnEqtl$gene)
jxnEqtl$snps = as.character(jxnEqtl$snps)

txEqtl = meTx$cis$eqtls
txEqtl$gene = as.character(txEqtl$gene)
txEqtl$snps = as.character(txEqtl$snps)

################################
# add gene annotation info #####
################################

geneEqtl$Symbol = rowRanges(rse_gene)$Symbol[match(geneEqtl$gene, rownames(rse_gene))]
geneEqtl$EnsemblGeneID = rowRanges(rse_gene)$ensemblID[match(geneEqtl$gene, rownames(rse_gene))]
geneEqtl$Type = "Gene"
geneEqtl$Class = "InGen"
geneEqtl = DataFrame(geneEqtl)
# geneEqtl$gene_type = rowRanges(rse_gene)$gene_type[match(geneEqtl$gene, rownames(rse_gene))]

exonEqtl$Symbol = rowRanges(rse_exon)$Symbol[match(exonEqtl$gene, rownames(rse_exon))]
exonEqtl$EnsemblGeneID = rowRanges(rse_exon)$ensemblID[match(exonEqtl$gene, rownames(rse_exon))]
exonEqtl$Type = "Exon"
exonEqtl$Class = "InGen"
exonEqtl = DataFrame(exonEqtl)
# exonEqtl$gene_type = rowRanges(rse_exon)$gene_type[match(exonEqtl$gene, rownames(rse_exon))]

jxnEqtl$Symbol = rowRanges(rse_jxn)$newGeneSymbol[match(jxnEqtl$gene, rownames(rse_jxn))]
jxnEqtl$EnsemblGeneID = rowRanges(rse_jxn)$newGeneID[match(jxnEqtl$gene, rownames(rse_jxn))]
jxnEqtl$Type = "Jxn"
jxnEqtl$Class = rowRanges(rse_jxn)$Class[match(jxnEqtl$gene, rownames(rse_jxn))]
jxnEqtl = DataFrame(jxnEqtl)
# jxnEqtl$gene_type = rowRanges(rse_jxn)$gene_type[match(jxnEqtl$gene, rownames(rse_jxn))]

txEqtl$Symbol = rowRanges(rse_tx)$gene_name[match(txEqtl$gene, rownames(rse_tx))]
txEqtl$EnsemblGeneID = ss(rowRanges(rse_tx)$gene_id[match(txEqtl$gene, rownames(rse_tx))],"\\.",1)
txEqtl$Type = "Tx"
txEqtl$Class = "InGen"
txEqtl = DataFrame(txEqtl)
# txEqtl$gene_type = rowRanges(rse_tx)$gene_type[match(txEqtl$gene, rownames(rse_tx))]


# merge
allEqtl = rbind(geneEqtl, exonEqtl, jxnEqtl, txEqtl)
allEqtl$gencodeTx = CharacterList(c(as.list(rowRanges(rse_gene)$gencodeTx[match(geneEqtl$gene, 
	rownames(rse_gene))]),
	as.list(rowRanges(rse_exon)$gencodeTx[match(exonEqtl$gene, rownames(rse_exon))]),
	as.list(rowRanges(rse_jxn)$gencodeTx[match(jxnEqtl$gene, rownames(rse_jxn))]),
	as.list(txEqtl$gene)))
    
allEqtl$Region <- "DLPFC"
save(allEqtl, file=here::here(
        "eQTL_full",
        "subset_hancock",
        "mergedEqtl_output_dlpfc_4features_subset_hancock.Rdata")
    ,compress=TRUE)

dim(allEqtl)
# [1] 6834   12

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.0.4 RC (2021-02-08 r79975)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-03-11
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version  date       lib source
#  annotate               1.68.0   2020-10-27 [1] Bioconductor
#  AnnotationDbi          1.52.0   2020-10-27 [1] Bioconductor
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.3)
#  Biobase              * 2.50.0   2020-10-27 [1] Bioconductor
#  BiocGenerics         * 0.36.0   2020-10-27 [1] Bioconductor
#  BiocParallel         * 1.24.1   2020-11-06 [1] Bioconductor
#  bit                    4.0.4    2020-08-04 [2] CRAN (R 4.0.3)
#  bit64                  4.0.5    2020-08-30 [2] CRAN (R 4.0.3)
#  bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.3)
#  blob                   1.2.1    2020-01-20 [2] CRAN (R 4.0.3)
#  cachem                 1.0.4    2021-02-13 [2] CRAN (R 4.0.4)
#  cli                    2.2.0    2020-11-20 [1] CRAN (R 4.0.3)
#  colorout               1.2-2    2020-05-09 [1] Github (jalvesaq/colorout@726d681)
#  colorspace             2.0-0    2020-11-11 [2] CRAN (R 4.0.3)
#  crayon                 1.4.1    2021-02-08 [1] CRAN (R 4.0.4)
#  DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.0.3)
#  DelayedArray           0.16.0   2020-10-27 [1] Bioconductor
#  digest                 0.6.27   2020-10-24 [1] CRAN (R 4.0.3)
#  dplyr                  1.0.2    2020-08-18 [1] CRAN (R 4.0.3)
#  edgeR                  3.32.0   2020-10-27 [1] Bioconductor
#  ellipsis               0.3.1    2020-05-15 [1] CRAN (R 4.0.3)
#  fansi                  0.4.1    2020-01-08 [1] CRAN (R 4.0.0)
#  fastmap                1.0.1    2019-10-08 [1] CRAN (R 4.0.0)
#  genefilter           * 1.72.0   2020-10-27 [1] Bioconductor
#  generics               0.1.0    2020-10-31 [1] CRAN (R 4.0.3)
#  GenomeInfoDb         * 1.26.2   2020-12-08 [1] Bioconductor
#  GenomeInfoDbData       1.2.4    2020-11-30 [2] Bioconductor
#  GenomicRanges        * 1.42.0   2020-10-27 [1] Bioconductor
#  ggplot2                3.3.3    2020-12-30 [1] CRAN (R 4.0.3)
#  glue                   1.4.2    2020-08-27 [1] CRAN (R 4.0.3)
#  googledrive            1.0.1    2020-05-05 [1] CRAN (R 4.0.3)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.3)
#  here                   1.0.1    2020-12-13 [1] CRAN (R 4.0.3)
#  htmltools              0.5.0    2020-06-16 [1] CRAN (R 4.0.3)
#  htmlwidgets            1.5.3    2020-12-10 [1] CRAN (R 4.0.3)
#  httpuv                 1.5.4    2020-06-06 [1] CRAN (R 4.0.3)
#  httr                   1.4.2    2020-07-20 [1] CRAN (R 4.0.3)
#  IRanges              * 2.24.1   2020-12-12 [1] Bioconductor
#  jaffelab             * 0.99.30  2021-01-12 [1] Github (LieberInstitute/jaffelab@42637ff)
#  jsonlite               1.7.2    2020-12-09 [2] CRAN (R 4.0.3)
#  later                  1.1.0.1  2020-06-05 [1] CRAN (R 4.0.3)
#  lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.4)
#  lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.0)
#  limma                  3.46.0   2020-10-27 [1] Bioconductor
#  locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.0.3)
#  magrittr               2.0.1    2020-11-17 [1] CRAN (R 4.0.3)
#  Matrix                 1.3-2    2021-01-06 [3] CRAN (R 4.0.4)
#  MatrixEQTL           * 2.3      2019-12-22 [1] CRAN (R 4.0.0)
#  MatrixGenerics       * 1.2.1    2021-01-30 [2] Bioconductor
#  matrixStats          * 0.58.0   2021-01-29 [1] CRAN (R 4.0.3)
#  memoise                2.0.0    2021-01-26 [2] CRAN (R 4.0.3)
#  mgcv                 * 1.8-34   2021-02-16 [3] CRAN (R 4.0.4)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.3)
#  nlme                 * 3.1-152  2021-02-04 [3] CRAN (R 4.0.4)
#  pillar                 1.4.7    2020-11-20 [1] CRAN (R 4.0.3)
#  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 4.0.0)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 4.0.3)
#  promises               1.1.1    2020-06-09 [1] CRAN (R 4.0.3)
#  purrr                  0.3.4    2020-04-17 [1] CRAN (R 4.0.0)
#  R6                     2.5.0    2020-10-28 [2] CRAN (R 4.0.3)
#  rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.0.0)
#  RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.0.3)
#  Rcpp                   1.0.5    2020-07-06 [1] CRAN (R 4.0.3)
#  RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.3)
#  rlang                  0.4.10   2020-12-30 [1] CRAN (R 4.0.3)
#  rmote                  0.3.4    2020-05-09 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.0.3)
#  RSQLite                2.2.3    2021-01-24 [2] CRAN (R 4.0.3)
#  S4Vectors            * 0.28.1   2020-12-09 [1] Bioconductor
#  scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.3)
#  segmented              1.3-1    2020-12-10 [1] CRAN (R 4.0.3)
#  servr                  0.21     2020-12-14 [1] CRAN (R 4.0.3)
#  sessioninfo          * 1.1.1    2018-11-05 [1] CRAN (R 4.0.3)
#  SummarizedExperiment * 1.20.0   2020-10-27 [1] Bioconductor
#  survival               3.2-7    2020-09-28 [1] CRAN (R 4.0.3)
#  sva                  * 3.38.0   2020-10-27 [2] Bioconductor
#  tibble                 3.0.4    2020-10-12 [1] CRAN (R 4.0.3)
#  tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.3)
#  vctrs                  0.3.6    2020-12-17 [1] CRAN (R 4.0.3)
#  withr                  2.3.0    2020-09-22 [1] CRAN (R 4.0.3)
#  xfun                   0.20     2021-01-06 [1] CRAN (R 4.0.3)
#  XML                    3.99-0.5 2020-07-23 [2] CRAN (R 4.0.3)
#  xtable                 1.8-4    2019-04-21 [2] CRAN (R 4.0.3)
#  XVector                0.30.0   2020-10-27 [1] Bioconductor
#  zlibbioc               1.36.0   2020-10-27 [1] Bioconductor
#
# [1] /users/lcollado/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
#
