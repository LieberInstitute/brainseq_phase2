# based on make_boxplots_hippo.R
library('SummarizedExperiment')
library('jaffelab')
library('MatrixEQTL')
library('sva')
library('RColorBrewer')
library('devtools')

######################
### load data ####
######################

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata", verbose = TRUE)

## keep adult samples - keep both regions
keepInd = which(colData(rse_gene)$Age > 13)
rse_gene = rse_gene[,keepInd]
rse_exon = rse_exon[,keepInd]
rse_jxn = rse_jxn[,keepInd]
rse_tx = rse_tx[,keepInd]

## extract pd and rpkms
pd = colData(rse_gene)
pd$Dx <- factor(ifelse(pd$Dx == 'Control', 'Control', 'SCZD'))
geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm



######################
### snp data ####
######################

## load SNP data
load("../genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda", verbose = TRUE)

### make mds and snp dimensions equal to N
###(repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum

## drop SNPs not mapping to hg38
keepIndex = which(!is.na(snpMap$chr_hg38))
snpMap = snpMap[keepIndex,]
snp = snp[keepIndex,]



################
## load table
load("eqtl_tables/matrixEqtl_output_interaction_4features.rda", verbose = TRUE)
allEqtl <- rbind(
    cbind(meGene$cis$eqtls, Type = 'Gene', Symbol = rowRanges(rse_gene)$Symbol[match(meGene$cis$eqtls$gene, rownames(rse_gene))]),
    cbind(meExon$cis$eqtls, Type = 'Exon', Symbol = rowRanges(rse_exon)$Symbol[match(meExon$cis$eqtls$gene, rownames(rse_exon))]),
    cbind(meJxn$cis$eqtls, Type = 'Jxn', Symbol = rowRanges(rse_jxn)$newGeneSymbol[match(meJxn$cis$eqtls$gene, rownames(rse_jxn))]),
    cbind(meTx$cis$eqtls, Type = 'Tx', Symbol = rowRanges(rse_tx)$gene_name[match(meTx$cis$eqtls$gene, rownames(rse_tx))])
)
dim(allEqtl)
allEqtl$snps <- as.character(allEqtl$snps)
allEqtl$gene <- as.character(allEqtl$gene)
allEqtl$Symbol <- as.character(allEqtl$Symbol)
save(allEqtl, file="eqtl_tables/mergedEqtl_output_interaction_4features_missingSomeAnnotation.rda", compress=TRUE)
#save(allEqtl, file="~/mergedEqtl_output_interaction_4features_missingSomeAnnotation.rda", compress=TRUE)
Sys.time()

# interaction = allEqtl[which(allEqtl$Type=="Gene"),]
# interaction$Symbol = as.character(interaction$Symbol)

################
## load expression


mod = model.matrix(~ Region + Dx + Sex + as.matrix(mds[,1:5]),
	data = pd)
colnames(mod)[5:9] = colnames(mds)[1:5]

load('eqtl_tables/rdas/pcs_4features_combined_regions_filtered_over13.rda', verbose = TRUE)

## residualize expression		
gExprs = log2(geneRpkm+1)		
gExprs = cleaningY(gExprs, cbind(mod, genePCs), P=2)

eExprs = log2(exonRpkm+1)		
eExprs = cleaningY(eExprs, cbind(mod, exonPCs), P=2)

jExprs = log2(jxnRp10m+1)		
jExprs = cleaningY(jExprs, cbind(mod, jxnPCs), P=2)

tExprs = log2(txTpm+1)		
tExprs = cleaningY(tExprs, cbind(mod, txPCs), P=2)


exprsAdj = rbind(gExprs,eExprs,jExprs,tExprs)



# pdf("interaction_top_eqtl_gene.pdf", h=6, w=6)
# par(mar=c(5,5,5,2), cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
# palette(brewer.pal(8,"Spectral"))
# ## plot
# for (i in 1:25) {
#     symi = interaction[i,"Symbol"]
#     symi[is.na(symi)]=""
#     snpi = interaction[i,"snps"]
#     feati = interaction[i,"gene"]
#     p_i = signif(interaction[i,"pvalue"],3)
#
#     boxplot(gExprs[feati,] ~ snp[snpi,],
#             xlab=snpi, ylab="Residualized Expression",
#             main=paste0(symi,"\n",feati," (Gene)"),
#             ylim = c(range(gExprs[feati,])), outline=FALSE)
#     points(gExprs[feati,] ~ jitter(snp[snpi,]+1),
#                pch=21,
#                bg=as.numeric(snp[snpi,])+2,cex=1.5)
#     legend("top",paste0("p=",p_i))
# }
# dev.off()



#################
### PGC index snps

pgc = read.csv("../eQTL_GWAS_riskSNPs/41588_2018_59_MOESM3_ESM.csv", stringsAsFactors=FALSE)

index = pgc$"Index.SNP..dbSNP.b141."
interaction2 = allEqtl[which(allEqtl$snps %in% index & allEqtl$FDR < 0.1),]
interaction2 = interaction2[order(interaction2$FDR, decreasing=FALSE),]
dim(interaction2)
# [1] 2 8
options(width = 160)
interaction2
## This is matching by name. To match by position check create_eqtl_table.R code by Emily Burke
#             snps                       gene statistic       pvalue         FDR        beta Type Symbol
# 1059226 rs324015 chr12:57154371-57154478(+) -4.788587 2.021204e-06 0.008699934 -0.09129091  Jxn   LRP1
# 327344  rs324015                    e714506 -4.137357 3.907839e-05 0.088492024 -0.07329732 Exon   LRP1

## Subset the 3 snps listed in create_eqtl_table.R
interaction2 <- allEqtl[grep("rs4144797|rs12293670|rs324015", allEqtl$snps), ]
dim(interaction2)
# [1] 16  8
interaction2 <- subset(interaction2, FDR < 0.1)
interaction2
#                             snps                         gene  statistic       pvalue          FDR        beta Type Symbol
# 327344                  rs324015                      e714506  -4.137357 3.907839e-05 8.849202e-02 -0.07329732 Exon   LRP1
# 988426  rs12293670:124612932:A:G chr11:124742693-124745500(+) -10.950543 5.266493e-26 7.296935e-21 -0.84182172  Jxn   <NA>
# 997933   rs4144797:233562197:T:C  chr2:232890181-232890219(-)  -7.727053 3.503959e-14 9.357701e-10 -0.73216733  Jxn   <NA>
# 1059226                 rs324015   chr12:57154371-57154478(+)  -4.788587 2.021204e-06 8.699934e-03 -0.09129091  Jxn   LRP1

options(width = 400)
rowRanges(rse_jxn)[interaction2$gene[interaction2$Type == 'Jxn']]
# GRanges object with 3 ranges and 18 metadata columns:
#                                seqnames              ranges strand | inGencode inGencodeStart inGencodeEnd      gencodeGeneID       ensemblID      Symbol gencodeStrand                           gencodeTx     numTx       Class startExon   endExon          newGeneID newGeneSymbol  isFusion        meanExprs    Length passExprsCut
#                                   <Rle>           <IRanges>  <Rle> | <logical>      <logical>    <logical>        <character>     <character> <character>   <character>                     <CharacterList> <integer> <character> <integer> <integer>        <character>   <character> <logical>        <numeric> <numeric>    <logical>
#   chr11:124742693-124745500(+)    chr11 124742693-124745500      + |     FALSE          FALSE        FALSE               <NA>            <NA>        <NA>          <NA>                                             0       Novel      <NA>      <NA>               <NA>          <NA>     FALSE 1.72236371055124       100         TRUE
#    chr2:232890181-232890219(-)     chr2 232890181-232890219      - |     FALSE          FALSE        FALSE               <NA>            <NA>        <NA>          <NA>                                             0       Novel      <NA>      <NA>               <NA>          <NA>     FALSE 1.38160646207221       100         TRUE
#     chr12:57154371-57154478(+)    chr12   57154371-57154478      + |      TRUE           TRUE         TRUE ENSG00000123384.13 ENSG00000123384        LRP1             + ENST00000243077.7,ENST00000554174.1         2       InGen    342567    342568 ENSG00000123384.13          LRP1     FALSE 168.867477061955       100         TRUE
#   -------
#   seqinfo: 25 sequences from an unspecified genome; no seqlengths
width(rowRanges(rse_jxn)[interaction2$gene[interaction2$Type == 'Jxn']])
# [1] 2808   39  108

rowRanges(rse_exon)[interaction2$gene[interaction2$Type == 'Exon']]
# GRanges object with 1 range and 11 metadata columns:
#           seqnames            ranges strand |    Length          gencodeID       ensemblID      gene_type      Symbol  EntrezID       Class        meanExprs     NumTx         gencodeTx passExprsCut
#              <Rle>         <IRanges>  <Rle> | <integer>        <character>     <character>    <character> <character> <integer> <character>        <numeric> <integer>   <CharacterList>    <logical>
#   e714506    chr12 57184839-57184990      + |       152 ENSG00000123384.13 ENSG00000123384 protein_coding        LRP1      4035       InGen 26.1312217740406         1 ENST00000243077.7         TRUE
#   -------
#   seqinfo: 25 sequences from an unspecified genome; no seqlengths

pdf("~/interaction_top_eqtl_PGC_indexSNPs.pdf", h=6, w=6, useDingbats=FALSE)
par(mar=c(5,5,5,2), cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
#palette(brewer.pal(8,"Spectral"))
palette(c('skyblue3', 'dark orange'))	
## plot
for (i in seq_len(nrow(interaction2))) {
	symi = interaction2[i,"Symbol"]
	symi[is.na(symi)]=""
	snpi = interaction2[i,"snps"]
	feati = interaction2[i,"gene"]
	p_i = signif(interaction2[i,"pvalue"],3)
	typei = interaction2[i,"Type"]
	
	boxplot(exprsAdj[feati,] ~ snp[snpi,] + substr(pd$Region, 1, 1),
			xlab=snpi, ylab="Residualized Expression", 
			main=paste0(symi,"\n",feati," (",typei,")"), 
			ylim = c(range(exprsAdj[feati,])), outline=FALSE)
	points(exprsAdj[feati,] ~ jitter(snp[snpi,] + ifelse(pd$Region == 'HIPPO', 4, 1)),
			   pch=21, 
			   bg= ifelse(pd$Region == 'HIPPO', 'skyblue3', 'dark orange'), cex=1.5)
	legend("topright",paste0("p=",p_i))
}
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.5.0 Patched (2018-04-30 r74679)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       US/Eastern
#  date     2018-09-05
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  annotate               1.58.0    2018-05-03 Bioconductor
#  AnnotationDbi          1.42.1    2018-05-17 Bioconductor
#  assertthat             0.2.0     2017-04-11 CRAN (R 3.5.0)
#  base                 * 3.5.0     2018-05-02 local
#  bindr                  0.1.1     2018-03-13 CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 CRAN (R 3.5.0)
#  Biobase              * 2.40.0    2018-05-02 Bioconductor
#  BiocGenerics         * 0.26.0    2018-05-03 Bioconductor
#  BiocParallel         * 1.14.2    2018-07-08 Bioconductor
#  bit                    1.1-14    2018-05-29 CRAN (R 3.5.0)
#  bit64                  0.9-7     2017-05-08 CRAN (R 3.5.0)
#  bitops                 1.0-6     2013-08-17 CRAN (R 3.5.0)
#  blob                   1.1.1     2018-03-25 CRAN (R 3.5.0)
#  colorout             * 1.2-0     2018-05-02 Github (jalvesaq/colorout@c42088d)
#  colorspace             1.3-2     2016-12-14 CRAN (R 3.5.0)
#  compiler               3.5.0     2018-05-02 local
#  crayon                 1.3.4     2017-09-16 CRAN (R 3.5.0)
#  datasets             * 3.5.0     2018-05-02 local
#  DBI                    1.0.0     2018-05-02 CRAN (R 3.5.0)
#  DelayedArray         * 0.6.2     2018-07-23 Bioconductor
#  devtools             * 1.13.6    2018-06-27 CRAN (R 3.5.0)
#  digest                 0.6.15    2018-01-28 CRAN (R 3.5.0)
#  dplyr                  0.7.6     2018-06-29 CRAN (R 3.5.0)
#  genefilter           * 1.62.0    2018-05-03 Bioconductor
#  GenomeInfoDb         * 1.16.0    2018-05-03 Bioconductor
#  GenomeInfoDbData       1.1.0     2018-04-17 Bioconductor
#  GenomicRanges        * 1.32.6    2018-07-20 Bioconductor
#  ggplot2                3.0.0     2018-07-03 CRAN (R 3.5.0)
#  glue                   1.3.0     2018-07-17 CRAN (R 3.5.0)
#  graphics             * 3.5.0     2018-05-02 local
#  grDevices            * 3.5.0     2018-05-02 local
#  grid                   3.5.0     2018-05-02 local
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 CRAN (R 3.5.0)
#  htmlwidgets            1.2       2018-04-19 CRAN (R 3.5.0)
#  httpuv                 1.4.5     2018-07-19 CRAN (R 3.5.0)
#  IRanges              * 2.14.10   2018-05-17 Bioconductor
#  jaffelab             * 0.99.21   2018-05-03 Github (LieberInstitute/jaffelab@7ed0ab7)
#  later                  0.7.3     2018-06-08 CRAN (R 3.5.0)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.5.0)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.5.0)
#  limma                  3.36.2    2018-06-21 Bioconductor
#  magrittr               1.5       2014-11-22 CRAN (R 3.5.0)
#  Matrix                 1.2-14    2018-04-13 CRAN (R 3.5.0)
#  MatrixEQTL           * 2.2       2018-01-13 CRAN (R 3.5.0)
#  matrixStats          * 0.54.0    2018-07-23 CRAN (R 3.5.0)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.5.0)
#  methods              * 3.5.0     2018-05-02 local
#  mgcv                 * 1.8-23    2018-01-21 CRAN (R 3.5.0)
#  mime                   0.5       2016-07-07 CRAN (R 3.5.0)
#  munsell                0.5.0     2018-06-12 CRAN (R 3.5.0)
#  nlme                 * 3.1-137   2018-04-07 CRAN (R 3.5.0)
#  parallel             * 3.5.0     2018-05-02 local
#  pillar                 1.3.0     2018-07-14 CRAN (R 3.5.0)
#  pkgconfig              2.0.1     2017-03-21 CRAN (R 3.5.0)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 CRAN (R 3.5.0)
#  purrr                  0.2.5     2018-05-29 CRAN (R 3.5.0)
#  R6                     2.2.2     2017-06-17 CRAN (R 3.5.0)
#  rafalib              * 1.0.0     2015-08-09 CRAN (R 3.5.0)
#  RColorBrewer         * 1.1-2     2014-12-07 CRAN (R 3.5.0)
#  Rcpp                   0.12.18   2018-07-23 CRAN (R 3.5.0)
#  RCurl                  1.95-4.11 2018-07-15 CRAN (R 3.5.0)
#  rlang                  0.2.1     2018-05-30 cran (@0.2.1)
#  rmote                * 0.3.4     2018-05-02 deltarho (R 3.5.0)
#  RSQLite                2.1.1     2018-05-06 CRAN (R 3.5.0)
#  S4Vectors            * 0.18.3    2018-06-13 Bioconductor
#  scales                 0.5.0     2017-08-24 CRAN (R 3.5.0)
#  segmented              0.5-3.0   2017-11-30 CRAN (R 3.5.0)
#  servr                  0.10      2018-05-30 CRAN (R 3.5.0)
#  splines                3.5.0     2018-05-02 local
#  stats                * 3.5.0     2018-05-02 local
#  stats4               * 3.5.0     2018-05-02 local
#  SummarizedExperiment * 1.10.1    2018-05-17 Bioconductor
#  survival               2.42-3    2018-04-16 CRAN (R 3.5.0)
#  sva                  * 3.28.0    2018-05-02 Bioconductor
#  tibble                 1.4.2     2018-01-22 CRAN (R 3.5.0)
#  tidyselect             0.2.4     2018-02-26 CRAN (R 3.5.0)
#  tools                  3.5.0     2018-05-02 local
#  utils                * 3.5.0     2018-05-02 local
#  withr                  2.1.2     2018-03-15 CRAN (R 3.5.0)
#  xfun                   0.3       2018-07-06 CRAN (R 3.5.0)
#  XML                    3.98-1.12 2018-07-15 CRAN (R 3.5.0)
#  xtable                 1.8-2     2016-02-05 CRAN (R 3.5.0)
#  XVector                0.20.0    2018-05-03 Bioconductor
#  zlibbioc               1.26.0    2018-05-02 Bioconductor
 