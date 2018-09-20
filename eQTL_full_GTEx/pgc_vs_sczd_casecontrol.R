library('data.table')
library('GenomicRanges')
library('devtools')

## Load DE results
message(paste(Sys.time(), 'loading SCZD case-control results'))
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/correlation/rda/out_info.Rdata', verbose = TRUE)

## Load sig eQTL results
message(paste(Sys.time(), 'loading significant eQTL results'))
load('rdas/merged_GTEx_BrainSeq_QTLs_dlpfc.Rdata', verbose = TRUE)
load('rdas/merged_GTEx_BrainSeq_QTLs_hippo.Rdata', verbose = TRUE)

## Load snp annotation
message(paste(Sys.time(), 'loading snp info'))
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda", verbose = TRUE)
rm(snp, mds)
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

## Load risk loci
message(paste(Sys.time(), 'loading risk loci info'))
riskLoci = read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/rAggr_results_179.csv", stringsAsFactors=FALSE)	# 10,981 snps
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## Subset to risk snps
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS)	# keep 9,736 snps
stopifnot(length(keepIndex) == 9736)
snpMap = snpMap[keepIndex,]
rm(keepIndex)

## Label if each SNP is a PGC (or in the neighborhood) snp
is_pgc <- function(DT) {
    DT$pgc <- DT$snps %in% snpMap$SNP
    print(table(DT$pgc))
    return(DT)
}

message(paste(Sys.time(), 'finding PGC snps in DLPFC'))
dlpfc <- lapply(dlpfc, is_pgc)
## Below each SNP might count more than once.
## The "TRUE" are the associations involving PGC snps.
## gene
#   FALSE    TRUE
# 1571579    6384
#
## exon
#   FALSE    TRUE
# 8879854   48455
#
## jxn
#   FALSE    TRUE
# 5236726   23400
#
## tx
#   FALSE    TRUE
# 2260036   11228
message(paste(Sys.time(), 'finding PGC snps in HIPPO'))
hippo <- lapply(hippo, is_pgc)
#   FALSE    TRUE
# 1071310    4423
#
#   FALSE    TRUE
# 6178899   27476
#
#   FALSE    TRUE
# 3939170   16079
#
#   FALSE    TRUE
# 1694175    8400

by_snp <- function(DT) {
#    setkey(DT, snps)
    summ <- DT[, sum(pgc), keyby = snps]
    print('Number of associations per SNP involving PGC SNPs')
    print(table(summ$V1))
    data.frame('no_pgc' = sum(summ$V1 == 0), 'pgc' = nrow(summ) - sum(summ$V1 == 0))
}

get_nums <- function(qtls, des) {
    res <- mapply(function(qtl, de) {
        snp_info <- by_snp(qtl)
        
        data.frame(
            'sczd_casecontrol_DE' = sum(de$adj.P.Val < 0.05),
            'sczd_casecontrol_notDE' = sum(de$adj.P.Val >= 0.05),
            'snp_assoc_pgc' = snp_info$pgc,
            'snp_assoc_no_pgc' = snp_info$no_pgc
        )
        
    }, qtls, des)
    t(res)
}

dlpfc_nums <- get_nums(dlpfc,
    list(
        'gene' = outGene$DLPFC_matchQSV,
        'exon' = outFeat$DLPFC$exon,
        'jxn' = outFeat$DLPFC$jxn,
        'tx' = outFeat$DLPFC$tx
    )
)
## gene
# [1] "Number of associations per SNP involving PGC SNPs"
#
#      0      1      2      3      4      5      6      7
# 965350   2704    755    277    238     53     18      2
#
## exon
# [1] "Number of associations per SNP involving PGC SNPs"
#
#       0       1       2       3       4       5       6       7       8       9
# 1182287     937     260     253     286     315     253     134     135     183
#      10      11      12      13      14      15      16      17      18      19
#     185     748     324      75      81      80      74      19      34      35
#      20      21      22      23      24      25      26      27      28      29
#      69      98      58      58      29      10      18      29      68      13
#      30      31      32      33      34      35      36      37      38      39
#       7      11      10       3       9       7       1       6       3       5
#      40      41      42      43      44      45      46      47      48      49
#      10       5       4       7      10       7      23      28       8       4
#      53      54      60
#       1       2      12
#
## jxn
# [1] "Number of associations per SNP involving PGC SNPs"
#
#       0       1       2       3       4       5       6       7       8       9
# 1373355    1593    1045     682     212     373     363     602      35      31
#      10      11      12      13      14      15      16      17      18      19
#      36      48      90      96      34      81      27      27      12      15
#      20      21      22      23      24      25      26      38      39      40
#      15      18       3      10       2       2       1       1       1       3
#      41      42      43
#       8       1       1
#
## tx
# [1] "Number of associations per SNP involving PGC SNPs"
#
#      0      1      2      3      4      5      6      7      8      9     10
# 958885   1588   1350    594    384    132    104     75    117     38     49
#     11     12
#      3      1
dlpfc_nums
#      sczd_casecontrol_DE sczd_casecontrol_notDE snp_assoc_pgc snp_assoc_no_pgc
# gene 245                 24407                  4047          965350
# exon 440                 396143                 5044          1182287
# jxn  37                  297144                 5468          1373355
# tx   6                   92726                  4435          958885

hippo_nums <- get_nums(hippo,
    list(
        'gene' = outGene$HIPPO_matchQSV,
        'exon' = outFeat$HIPPO$exon,
        'jxn' = outFeat$HIPPO$jxn,
        'tx' = outFeat$HIPPO$tx
    )
)

## gene
# [1] "Number of associations per SNP involving PGC SNPs"
#
#      0      1      2      3      4      5      6
# 679518   1536    762    256    116     25      1
#
## exon
# [1] "Number of associations per SNP involving PGC SNPs"
#
#      0      1      2      3      4      5      6      7      8      9     10
# 891293    760    356    362    559    298    149    130     68     92     43
#     11     12     13     14     15     16     17     18     19     20     21
#    138     83     57     77     57     41     21     42     60     91     81
#     22     23     24     25     26     27     28     29     30     31     32
#     76     37     23     54     12     22      7      3      2      3      1
#     33     34     35
#      2      2      1
#
## jxn
# [1] "Number of associations per SNP involving PGC SNPs"
#
#       0       1       2       3       4       5       6       7       8       9
# 1138472    1330    1114     674     348      78     201      71      29      33
#      10      11      12      13      14      15      16      17      18      19
#      53     116      94     108      24       8       2       2      33       5
#      20      21      22      25      34      35      36      37
#       7       6       5       1       1       4       9       1
#
## tx
# [1] "Number of associations per SNP involving PGC SNPs"
#
#      0      1      2      3      4      5      6      7      8      9     10
# 742918   1346    623    353    245    222     75    129     77      5     17
#     11     12
#     41      2

hippo_nums
#      sczd_casecontrol_DE sczd_casecontrol_notDE snp_assoc_pgc snp_assoc_no_pgc
# gene 48                  24604                  2696          679518
# exon 197                 396386                 3810          891293
# jxn  41                  297140                 4357          1138472
# tx   0                   92732                  3135          742918


apply(dlpfc_nums, 1, function(x) {
    x <- unlist(x)
    print(names(x))
    print(x)
    phyper(x['sczd_casecontrol_DE'], x['snp_assoc_pgc'], x['sczd_case_control_notDE'], x['snp_assoc_no_pgc'])
})

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
