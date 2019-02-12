# qrsh -l mem_free=150G,h_vmem=150G,h_fsize=100G
library('SummarizedExperiment')
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

#### Filter Data


### First: gene
stopifnot(identical(names(rowRanges(rse_gene)), names(rowRanges(bsp1_gene))))

m_gene <- match(names(rowRanges(rse_gene)[rowRanges(rse_gene)$passExprsCut]), names(rowRanges(bsp1_gene)))
stopifnot(!any(is.na(m_gene)))
bsp1_gene <- bsp1_gene[m_gene, ]

## Check with the filtered data
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata', verbose = TRUE)
stopifnot(identical(names(rowRanges(rse_gene)), names(rowRanges(bsp1_gene))))

## Add RPKM
assays(bsp1_gene)$rpkm <- recount::getRPKM(bsp1_gene, 'Length')
rowRanges(bsp1_gene)$meanExprs <- rowMeans(assays(bsp1_gene)$rpkm)
colData(bsp1_gene)$Region <- 'DLPFC'

## Ok, save it now
save(bsp1_gene, file = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/bsp1_gene.Rdata')



### Second: exon
stopifnot(identical(nrow(rse_exon), nrow(bsp1_exon)))
identical(names(rowRanges(rse_exon)), names(rowRanges(bsp1_exon)))

## Ok, there are 11 exons that match the genomic coordinates but not the strands, will drop those
## well, only 2/11 passed the exprs cutoff in BSP2
ov <- findOverlaps(rse_exon, bsp1_exon, type = 'equal')
table(table(queryHits(ov)))
missing <- which(!seq_len(nrow(rse_exon)) %in% queryHits(ov))

ov2 <- findOverlaps(bsp1_exon, rse_exon, type = 'equal')
table(table(queryHits(ov2)))
missing2 <- which(!seq_len(nrow(bsp1_exon)) %in% queryHits(ov2))


options(width = 500)
rowRanges(rse_exon_unfiltered)[missing]
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

## Add RPKM
assays(bsp1_exon)$rpkm <- recount::getRPKM(bsp1_exon, 'Length')
rowRanges(bsp1_exon)$meanExprs <- rowMeans(assays(bsp1_exon)$rpkm)
colData(bsp1_exon)$Region <- 'DLPFC'

## Save it finally
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

## Subjest
bsp1_jxn <- bsp1_jxn[unique(subjectHits(ov_jxn)), ]

## Check the jxn names
ov_jxn_rev <- findOverlaps(bsp1_jxn, rse_jxn[rowRanges(rse_jxn)$passExprsCut], type = 'equal')
stopifnot(all(queryHits(ov_jxn_rev) == seq_len(nrow(bsp1_jxn))))

identical(names(bsp1_jxn), names(rse_jxn[rowRanges(rse_jxn)$passExprsCut])[subjectHits(ov_jxn_rev)])

## Add RP10M
rowRanges(bsp1_jxn)$Length <- 100
assays(bsp1_jxn)$rp10m <- recount::getRPKM(bsp1_jxn, 'Length')
rowRanges(bsp1_jxn)$meanExprs <- rowMeans(assays(bsp1_jxn)$rp10m)
colData(bsp1_jxn)$Region <- 'DLPFC'

## Save for later use
save(bsp1_jxn, file = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/bsp1_jxn.Rdata')


### Fourth: tx



stopifnot(identical(names(rowRanges(rse_tx)), names(rowRanges(bsp1_tx))))




rowRanges(bsp1_tx)$meanExprs <- rowMeans(assay(bsp1_tx))
colData(bsp1_tx)$Region <- 'DLPFC'




# load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata', verbose = TRUE)
# load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata', verbose = TRUE)
# load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata', verbose = TRUE)
# load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata', verbose = TRUE)
#
#
# ov3 <- findOverlaps(rse_exon, bsp1_exon, type = 'equal')
# table(table(queryHits(ov3)))
# missing3 <- which(!seq_len(nrow(rse_exon)) %in% queryHits(ov3))
# rowRanges(rse_exon)[missing3]


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
