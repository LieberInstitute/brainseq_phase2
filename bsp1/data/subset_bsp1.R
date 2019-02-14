# qrsh -l mem_free=150G,h_vmem=150G,h_fsize=100G
library('SummarizedExperiment')
library('recount')
library('sessioninfo')

dir.create('/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data', recursive = TRUE, showWarnings = FALSE)


## Load BSP1
load('/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseGene_merged_n732.rda', verbose = TRUE)
load('/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseExon_merged_n732.rda', verbose = TRUE)
load('/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseJxn_merged_n732.rda', verbose = TRUE)
load('/dcl01/ajaffe/data/lab/brainseq_phase1/count_data/dlpfc_polyA_brainseq_phase1_hg38_rseTx_merged_n732.rda', verbose = TRUE)
bsp1_gene <- rse_gene
bsp1_exon <- rse_exon
bsp1_jxn <- rse_jxn
bsp1_tx <- rse_tx

## Load BSP2
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_gene_unfiltered.Rdata', verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_exon_unfiltered.Rdata', verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_jxn_unfiltered.Rdata', verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_tx_unfiltered.Rdata', verbose = TRUE)

## Load genotype data
load('/dcl01/ajaffe/data/lab/brainseq_phase1/genotype_data/brainseq_phase1_Genotypes_n732.rda', verbose = TRUE)

fix_pd <- function(rse) {
    ## Only do this for the SCZD vs control analysis, not for all...
    # rse <- rse[, colData(rse)$Dx %in% c('Control', 'Schizo')]
    m <- match(rse$BrNum, rownames(mds))
    stopifnot(!any(is.na(m)))

    colData(rse) <- cbind(colData(rse), mds[m, 1:10])
    
    return(rse)
}

#### Filter Data
### First: gene
stopifnot(identical(names(rowRanges(rse_gene)), names(rowRanges(bsp1_gene))))

m_gene <- match(names(rowRanges(rse_gene)[rowRanges(rse_gene)$passExprsCut]), names(rowRanges(bsp1_gene)))
stopifnot(!any(is.na(m_gene)))
bsp1_gene <- bsp1_gene[m_gene, ]

## Not needed
# ## Check with the filtered data
# load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata', verbose = TRUE)
# stopifnot(identical(names(rowRanges(rse_gene)), names(rowRanges(bsp1_gene))))
stopifnot(identical(names(rowRanges(rse_gene)[rowRanges(rse_gene)$passExprsCut]), names(rowRanges(bsp1_gene))))

## Add RPKM and fix the pheno columns
bsp1_gene <- fix_pd(bsp1_gene)
assays(bsp1_gene)$rpkm <- recount::getRPKM(bsp1_gene, 'Length')
rowRanges(bsp1_gene)$meanExprs <- rowMeans(assays(bsp1_gene)$rpkm)

## How much is above the BSP2 expr cutoff?
table(rowRanges(bsp1_gene)$meanExprs > 0.25)
# FALSE  TRUE
#  5187 19465
table(rowRanges(bsp1_gene)$meanExprs > 0.25) / nrow(bsp1_gene) * 100
#    FALSE     TRUE
# 21.04089 78.95911

## Ok, save it now
dim(bsp1_gene)
# [1] 24652   732
save(bsp1_gene, file = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/bsp1_gene.Rdata')






### Second: exon
stopifnot(identical(nrow(rse_exon), nrow(bsp1_exon)))
identical(names(rowRanges(rse_exon)), names(rowRanges(bsp1_exon)))

## Ok, there are 11 exons that match the genomic coordinates but not the strands, will drop those
## well, only 2/11 passed the exprs cutoff in BSP2
ov <- findOverlaps(rse_exon, bsp1_exon, type = 'equal')
table(table(queryHits(ov)))
#      1
# 571612
missing <- which(!seq_len(nrow(rse_exon)) %in% queryHits(ov))

ov2 <- findOverlaps(bsp1_exon, rse_exon, type = 'equal')
table(table(queryHits(ov2)))
#      1
# 571612
missing2 <- which(!seq_len(nrow(bsp1_exon)) %in% queryHits(ov2))


options(width = 500)
rowRanges(rse_exon)[missing]
# GRanges object with 11 ranges and 11 metadata columns:
#            seqnames              ranges strand |    Length         gencodeID       ensemblID      gene_type      Symbol  EntrezID       Class            meanExprs     NumTx         gencodeTx passExprsCut
#               <Rle>           <IRanges>  <Rle> | <integer>       <character>     <character>    <character> <character> <integer> <character>            <numeric> <integer>   <CharacterList>    <logical>
#     e82220     chr1 170151378-170151462      + |        85 ENSG00000263390.1 ENSG00000263390          miRNA                  <NA>       InGen  0.00367306761698502         0              <NA>        FALSE
#    e181558     chr2 206783234-206783308      + |        75 ENSG00000263468.1 ENSG00000263468          miRNA                  <NA>       InGen    0.172298062345905         0              <NA>        FALSE
#    e434110     chr7   30289794-30289890      - |        97 ENSG00000283237.1 ENSG00000283237          miRNA                  <NA>       InGen                    0         1 ENST00000637270.1        FALSE
#    e435118     chr7   32732981-32733077      - |        97 ENSG00000283461.1 ENSG00000283461          miRNA                  <NA>       InGen  0.00905475475855007         1 ENST00000637473.1        FALSE
#    e514496     chr8 123348034-123348130      - |        97 ENSG00000283172.1 ENSG00000283172          miRNA                  <NA>       InGen    0.453936954830773         1 ENST00000636914.1         TRUE
#    e604240    chr10 121949412-121949487      - |        76 ENSG00000273767.1 ENSG00000273767       misc_RNA                  <NA>       InGen   0.0872228115005848         0              <NA>        FALSE
#    e683741    chr12     6666648-6668115      + |      1468 ENSG00000219410.5 ENSG00000219410      antisense                  <NA>       InGen  0.00188420449349526         1 ENST00000396799.3        FALSE
#    e718157    chr12   64622509-64622605      - |        97 ENSG00000266016.3 ENSG00000266016          miRNA                  <NA>       InGen                    0         1 ENST00000584743.3        FALSE
#    e862727    chr16     3384459-3384941      + |       483 ENSG00000262621.4 ENSG00000262621 protein_coding                  <NA>       InGen    0.301426795813603         1 ENST00000618352.1         TRUE
#    e862730    chr16     3397192-3397745      + |       554 ENSG00000262621.4 ENSG00000262621 protein_coding                  <NA>       InGen    0.298294154766214         1 ENST00000618352.1        FALSE
#   e1157108     chrX   55451495-55451582      + |        88 ENSG00000266328.3 ENSG00000266328          miRNA                  <NA>       InGen 0.000718119673219427         0              <NA>        FALSE
#   -------
#   seqinfo: 25 sequences from an unspecified genome; no seqlengths
  
rowRanges(bsp1_exon)[missing2]
# GRanges object with 11 ranges and 15 metadata columns:
#            seqnames              ranges strand |    Length          gencodeID       ensemblID                      gene_type      Symbol  EntrezID       Class            meanExprs                        coord   NumENSE                                        exon_gencodeID   NumLIBD                                             exon_libdID     NumTx                                                                                                                     gencodeTx
#               <Rle>           <IRanges>  <Rle> | <integer>        <character>     <character>                    <character> <character> <integer> <character>            <numeric>                  <character> <integer>                                           <character> <integer>                                             <character> <integer>                                                                                                                   <character>
#     e82219     chr1 170151378-170151462      - |        85  ENSG00000283340.1 ENSG00000283340                          miRNA   MIR3119-1 100422839       InGen  0.00123494697100833  chr1:170151378-170151462(-)         1                                     ENSE00003800987.1         0                                                                 1                                                                                                             ENST00000637673.1
#    e181557     chr2 206783234-206783308      - |        75  ENSG00000283469.1 ENSG00000283469                          miRNA   MIR3130-1 100422993       InGen  0.00810798823668647  chr2:206783234-206783308(-)         1                                     ENSE00003796367.1         0                                                                 1                                                                                                             ENST00000637816.1
#    e434111     chr7   30289794-30289890      + |        97  ENSG00000207771.1 ENSG00000207771                          miRNA    MIR550A1    693133       InGen   0.0161333761108822    chr7:30289794-30289890(+)         1                                     ENSE00001808817.1         1                                                 e434111         1                                                                                                             ENST00000385037.1
#    e435119     chr7   32732981-32733077      + |        97  ENSG00000207573.1 ENSG00000207573                          miRNA    MIR550A2    693134       InGen   0.0300062474200371    chr7:32732981-32733077(+)         1                                     ENSE00001808098.1         1                                                 e435119         1                                                                                                             ENST00000384841.1
#    e514497     chr8 123348034-123348130      + |        97  ENSG00000207704.2 ENSG00000207704                          miRNA   MIR548AA1    693130       InGen    0.137786975806104  chr8:123348034-123348130(+)         1                                     ENSE00001499978.2         1                                                 e514497         1                                                                                                             ENST00000384971.2
#    e604238    chr10 121949412-121949487      + |        76  ENSG00000226864.2 ENSG00000226864 transcribed_unitary_pseudogene    ATE1-AS1 100130887       InGen    0.175151417519368 chr10:121949412-121949487(+)         1                                     ENSE00003793137.1         0                                                                 1                                                                                                             ENST00000636460.1
#    e683765    chr12     6666648-6668115      - |      1468 ENSG00000126746.17 ENSG00000126746                 protein_coding      ZNF384    171017       InGen     3.91773936635725     chr12:6666648-6668115(-)         1                                     ENSE00002321416.1         1                                                 e683765         1                                                                                                             ENST00000396801.7
#    e718158    chr12   64622509-64622605      + |        97  ENSG00000207546.1 ENSG00000207546                          miRNA     MIR548C    693129       InGen  0.00601548503479784   chr12:64622509-64622605(+)         1                                     ENSE00001808994.1         1                                                 e718158         1                                                                                                             ENST00000384815.1
#    e862739    chr16     3384459-3384941      - |       483 ENSG00000140987.19 ENSG00000140987                 protein_coding     ZSCAN32     54925       InGen     1.34889934773721     chr16:3384459-3384941(-)         3 ENSE00003471480.1;ENSE00003634190.1;ENSE00003465353.1         7 e862739;e862743;e862750;e862757;e862761;e862764;e862771         7 ENST00000396846.7;ENST00000618425.4;ENST00000304926.7;ENST00000396852.8;ENST00000576500.5;ENST00000439568.2;ENST00000573327.5
#    e862773    chr16     3397192-3397745      - |       554 ENSG00000140987.19 ENSG00000140987                 protein_coding     ZSCAN32     54925       InGen     1.63683142496467     chr16:3397192-3397745(-)         1                                     ENSE00001730526.1         2                                         e862773;e862790         2                                                                                           ENST00000574940.5;ENST00000573719.1
#   e1157107     chrX   55451495-55451582      - |        88  ENSG00000283334.1 ENSG00000283334                          miRNA   MIR4536-2 100616155       InGen 0.000787436229558355    chrX:55451495-55451582(-)         1                                     ENSE00003793502.1         0                                                                 1                                                                                                             ENST00000636519.1
#   -------
#   seqinfo: 25 sequences from an unspecified genome; no seqlengths

## checked the first one that is on both (different strands) that passed the exprs cutoff
## https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A123348034-123348130&hgsid=710804097_cORrkdA10wn9gCJrihUsKeauXoDJ
## indeed, it exists on both strands 

## Match by chr position/strand
ov_exon <- findOverlaps(rse_exon[rowRanges(rse_exon)$passExprsCut], bsp1_exon, type = 'equal')
stopifnot(length(queryHits(ov_exon)) == sum(rowRanges(rse_exon)$passExprsCut) - 2)
stopifnot(length(queryHits(ov_exon)) == length(unique(queryHits(ov_exon))))
stopifnot(length(queryHits(ov_exon)) == length(unique(subjectHits(ov_exon))))

## Subset now
bsp1_exon <- bsp1_exon[subjectHits(ov_exon), ]

## Check that the LIBD exon names are the same ones we used in BSP2
ov_exon_rev <- findOverlaps(bsp1_exon, rse_exon[rowRanges(rse_exon)$passExprsCut], type = 'equal')
stopifnot(all(queryHits(ov_exon_rev) == seq_len(nrow(bsp1_exon))))

identical(names(bsp1_exon), names(rse_exon[rowRanges(rse_exon)$passExprsCut])[subjectHits(ov_exon_rev)])
# [1] FALSE
table(names(bsp1_exon) == names(rse_exon[rowRanges(rse_exon)$passExprsCut])[subjectHits(ov_exon_rev)])
# FALSE   TRUE
#  1341 395240

## The ranges look ok
stopifnot(identical(ranges(bsp1_exon, use.names = FALSE), ranges(rse_exon[rowRanges(rse_exon)$passExprsCut], use.names = FALSE)[subjectHits(ov_exon_rev)] ))

## So now we can force the exon LIBD names (IDs) to match
names(bsp1_exon) <- names(rse_exon[rowRanges(rse_exon)$passExprsCut])[subjectHits(ov_exon_rev)]

## Check that it all looks good
m_exon <- match(names(bsp1_exon), names(rse_exon))
stopifnot(identical(ranges(bsp1_exon), ranges(rse_exon[m_exon])))

## Add RPKM and fix the pheno columns
bsp1_exon <- fix_pd(bsp1_exon)
assays(bsp1_exon)$rpkm <- recount::getRPKM(bsp1_exon, 'Length')
rowRanges(bsp1_exon)$meanExprs <- rowMeans(assays(bsp1_exon)$rpkm)

## How much is above the BSP2 expr cutoff?
table(rowRanges(bsp1_exon)$meanExprs > 0.30)
# FALSE   TRUE
# 44343 352238
table(rowRanges(bsp1_exon)$meanExprs > 0.30) / nrow(bsp1_exon) * 100
#    FALSE     TRUE
# 11.18132 88.81868

## Save it finally
dim(bsp1_exon)
# [1] 396581    732
save(bsp1_exon, file = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/bsp1_exon.Rdata')





### Thid: jxn
dim(rse_jxn)
# [1] 837077    900
dim(bsp1_jxn)
# [1] 1832800     732

table(countOverlaps(rse_jxn, bsp1_jxn, type = 'equal'))
#      0      1
# 167602 669475
table(countOverlaps(rse_jxn[rowRanges(rse_jxn)$passExprsCut], bsp1_jxn, type = 'equal'))
#     0      1
# 23577 273604

## Match by chr position/strand
ov_jxn <- findOverlaps(rse_jxn[rowRanges(rse_jxn)$passExprsCut], bsp1_jxn, type = 'equal')
stopifnot(length(queryHits(ov_jxn)) == sum(rowRanges(rse_jxn)$passExprsCut) - 23577)
stopifnot(length(queryHits(ov_jxn)) == length(unique(queryHits(ov_jxn))))
stopifnot(length(queryHits(ov_jxn)) == length(unique(subjectHits(ov_jxn)))) ## this isn't true
## due to the fact that BSP1 is unstranded and BSP2 is stranded

## Subset
bsp1_jxn <- bsp1_jxn[unique(subjectHits(ov_jxn)), ]

## Check the jxn names
ov_jxn_rev <- findOverlaps(bsp1_jxn, rse_jxn[rowRanges(rse_jxn)$passExprsCut], type = 'equal')
table(countOverlaps(bsp1_jxn, rse_jxn[rowRanges(rse_jxn)$passExprsCut], type = 'equal'))
#      1      2
# 267872   2866
stopifnot(267872 + 2866 - nrow(bsp1_jxn) == 0)
stopifnot(267872 + 2866 * 2 - length(queryHits(ov_jxn_rev)) == 0)
stopifnot(all(queryHits(ov_jxn_rev) %in% seq_len(nrow(bsp1_jxn))))


## Remove the strand info part of the names when checking
stopifnot(identical(
    gsub('\\(.\\)', '', names(bsp1_jxn)[queryHits(ov_jxn_rev)]),
    gsub('\\(.\\)', '', names(rse_jxn[rowRanges(rse_jxn)$passExprsCut])[subjectHits(ov_jxn_rev)])
))

## Add RP10M and fix the pheno columns
bsp1_jxn <- fix_pd(bsp1_jxn)
rowRanges(bsp1_jxn)$Length <- 100
assays(bsp1_jxn)$rp10m <- recount::getRPKM(bsp1_jxn, 'Length')
rowRanges(bsp1_jxn)$meanExprs <- rowMeans(assays(bsp1_jxn)$rp10m)

## How much is above the BSP2 expr cutoff?
table(rowRanges(bsp1_jxn)$meanExprs > 0.46)
# FALSE   TRUE
# 41330 229408
table(rowRanges(bsp1_jxn)$meanExprs > 0.46) / nrow(bsp1_jxn) * 100
#    FALSE     TRUE
# 15.26568 84.73432

## Save for later use
dim(bsp1_jxn)
# [1] 270738    732
save(bsp1_jxn, file = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/bsp1_jxn.Rdata')




### Fourth: tx


## Ok, names are the same
stopifnot(identical(names(rowRanges(rse_tx)), names(rowRanges(bsp1_tx))))


## Code adapted from the gene case
m_tx <- match(names(rowRanges(rse_tx)[rowRanges(rse_tx)$passExprsCut]), names(rowRanges(bsp1_tx)))
stopifnot(!any(is.na(m_tx)))
bsp1_tx <- bsp1_tx[m_tx, ]

## Check again
stopifnot(identical(names(rowRanges(rse_tx)[rowRanges(rse_tx)$passExprsCut]), names(rowRanges(bsp1_tx))))

## Add mean TPM and fix the pheno columns
bsp1_tx <- fix_pd(bsp1_tx)
rowRanges(bsp1_tx)$meanExprs <- rowMeans(assay(bsp1_tx))

## How much is above the BSP2 expr cutoff?
table(rowRanges(bsp1_tx)$meanExprs > 0.38)
# FALSE  TRUE
# 24660 68072
table(rowRanges(bsp1_tx)$meanExprs > 0.38) / nrow(bsp1_tx) * 100
#    FALSE     TRUE
# 26.59276 73.40724


## Save for later use
dim(bsp1_tx)
# [1] 92732   732
save(bsp1_tx, file = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/bsp1_tx.Rdata')


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
#  date     2019-02-14
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  acepack                1.4.1     2016-10-29 [2] CRAN (R 3.5.0)
#  AnnotationDbi          1.44.0    2018-10-30 [1] Bioconductor
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  backports              1.1.3     2018-12-14 [2] CRAN (R 3.5.1)
#  base64enc              0.1-3     2015-07-28 [2] CRAN (R 3.5.0)
#  bibtex                 0.4.2     2017-06-30 [1] CRAN (R 3.5.0)
#  bindr                  0.1.1     2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 [1] CRAN (R 3.5.0)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.5    2019-01-04 [1] Bioconductor
#  biomaRt                2.38.0    2018-10-30 [1] Bioconductor
#  Biostrings             2.50.2    2019-01-03 [1] Bioconductor
#  bit                    1.1-14    2018-05-29 [2] CRAN (R 3.5.0)
#  bit64                  0.9-7     2017-05-08 [2] CRAN (R 3.5.0)
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  blob                   1.1.1     2018-03-25 [2] CRAN (R 3.5.0)
#  BSgenome               1.50.0    2018-10-30 [1] Bioconductor
#  bumphunter             1.24.5    2018-12-01 [1] Bioconductor
#  checkmate              1.9.1     2019-01-15 [1] CRAN (R 3.5.1)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  cluster                2.0.7-1   2018-04-13 [3] CRAN (R 3.5.1)
#  codetools              0.2-15    2016-10-05 [3] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  data.table             1.12.0    2019-01-13 [1] CRAN (R 3.5.1)
#  DBI                    1.0.0     2018-05-02 [2] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  derfinder              1.16.1    2018-12-03 [1] Bioconductor
#  derfinderHelper        1.16.1    2018-12-03 [1] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  doRNG                  1.7.1     2018-06-22 [2] CRAN (R 3.5.1)
#  downloader             0.4       2015-07-09 [1] CRAN (R 3.5.0)
#  dplyr                  0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
#  foreach                1.4.4     2017-12-12 [2] CRAN (R 3.5.0)
#  foreign                0.8-71    2018-07-20 [3] CRAN (R 3.5.1)
#  Formula                1.2-3     2018-05-03 [2] CRAN (R 3.5.0)
#  GenomeInfoDb         * 1.18.1    2018-11-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicAlignments      1.18.1    2019-01-04 [1] Bioconductor
#  GenomicFeatures        1.34.1    2018-11-03 [1] Bioconductor
#  GenomicFiles           1.18.0    2018-10-30 [1] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  GEOquery               2.50.5    2018-12-22 [1] Bioconductor
#  ggplot2                3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.0     2018-07-17 [1] CRAN (R 3.5.1)
#  gridExtra              2.3       2017-09-09 [2] CRAN (R 3.5.0)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  Hmisc                  4.1-1     2018-01-03 [1] CRAN (R 3.5.0)
#  hms                    0.4.2     2018-03-10 [2] CRAN (R 3.5.0)
#  htmlTable              1.13.1    2019-01-07 [2] CRAN (R 3.5.1)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  httr                   1.4.0     2018-12-11 [1] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  iterators              1.0.10    2018-07-13 [2] CRAN (R 3.5.1)
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  knitr                  1.21      2018-12-10 [2] CRAN (R 3.5.1)
#  later                  0.7.5     2018-09-18 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  latticeExtra           0.6-28    2016-02-09 [2] CRAN (R 3.5.0)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma                  3.38.3    2018-12-02 [1] Bioconductor
#  locfit                 1.5-9.1   2013-04-20 [2] CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  memoise                1.1.0     2017-04-21 [2] CRAN (R 3.5.0)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.0)
#  nnet                   7.3-12    2016-02-02 [3] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  pkgmaker               0.27      2018-05-25 [2] CRAN (R 3.5.0)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  prettyunits            1.0.2     2015-07-13 [1] CRAN (R 3.5.0)
#  progress               1.2.0     2018-06-14 [1] CRAN (R 3.5.1)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                  0.2.5     2018-05-29 [2] CRAN (R 3.5.0)
#  qvalue                 2.14.1    2019-01-10 [1] Bioconductor
#  R6                     2.3.0     2018-10-04 [2] CRAN (R 3.5.1)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  readr                  1.3.1     2018-12-21 [1] CRAN (R 3.5.1)
#  recount              * 1.8.1     2018-12-03 [1] Bioconductor
#  registry               0.5       2017-12-03 [2] CRAN (R 3.5.0)
#  rentrez                1.2.1     2018-03-05 [1] CRAN (R 3.5.0)
#  reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.5.0)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  rngtools               1.3.1     2018-05-15 [2] CRAN (R 3.5.0)
#  rpart                  4.1-13    2018-02-23 [3] CRAN (R 3.5.1)
#  Rsamtools              1.34.0    2018-10-30 [1] Bioconductor
#  RSQLite                2.1.1     2018-05-06 [2] CRAN (R 3.5.0)
#  rstudioapi             0.9.0     2019-01-09 [2] CRAN (R 3.5.1)
#  rtracklayer            1.42.1    2018-11-22 [1] Bioconductor
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  servr                  0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  stringi                1.2.4     2018-07-20 [2] CRAN (R 3.5.1)
#  stringr                1.3.1     2018-05-10 [1] CRAN (R 3.5.0)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  survival               2.43-3    2018-11-26 [3] CRAN (R 3.5.1)
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyr                  0.8.2     2018-10-28 [2] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  VariantAnnotation      1.28.10   2019-01-21 [1] Bioconductor
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  XML                    3.98-1.16 2018-08-19 [2] CRAN (R 3.5.1)
#  xml2                   1.2.0     2018-01-24 [2] CRAN (R 3.5.0)
#  xtable                 1.8-3     2018-08-29 [2] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
