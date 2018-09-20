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
    summ <- DT[, sum(pgc), keyby = snps]
    print('Number of associations per SNP involving PGC SNPs')
    print(table(summ$V1))
    data.frame('no_pgc' = sum(summ$V1 == 0), 'pgc' = nrow(summ) - sum(summ$V1 == 0))
}

by_feat <- function(DT) {
    summ <- DT[, sum(pgc), keyby = gene]
    print('Number of associations per feature involving PGC SNPs')
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

get_nums_feat <- function(qtls, des) {
    res <- mapply(function(qtl, de) {
        snp_info <- by_feat(qtl)
        
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

## Something is wrong here
phyp_qtl <- function(x) {
    x <- unlist(x)
    phyper(x['snp_assoc_pgc'], x['sczd_casecontrol_DE'], x['sczd_casecontrol_notDE'], x['snp_assoc_no_pgc'] + x['snp_assoc_pgc'])    
}
apply(dlpfc_nums, 1, phyp_qtl)
## Should summarize at the feature level (I think)


dlpfc_nums_feat <- get_nums_feat(dlpfc,
    list(
        'gene' = outGene$DLPFC_matchQSV,
        'exon' = outFeat$DLPFC$exon,
        'jxn' = outFeat$DLPFC$jxn,
        'tx' = outFeat$DLPFC$tx
    )
)
## gene
# [1] "Number of associations per feature involving PGC SNPs"
#
#     0     1     2     3     4     5     6     7     8     9    10    11    12
# 14116    21    11     9     6     4     2     2     5     4     1     1     1
#    14    15    16    18    19    20    21    22    23    25    26    28    29
#     1     3     1     2     2     2     5     1     3     1     1     2     1
#    31    32    34    35    36    39    45    46    47    51    57    58    59
#     1     1     1     1     2     3     1     2     2     1     1     1     1
#    61    64    66    69    70    79    82    83    85    88    89    90    98
#     1     1     1     1     1     2     1     1     1     1     1     1     2
#   100   105   108   113   114   115   122   126   134   144   155   160   185
#     1     1     1     1     1     1     2     1     1     1     2     2     1
#   227   238   268   785
#     1     1     1     1
#
## exon
# [1] "Number of associations per feature involving PGC SNPs"
#
#      0      1      2      3      4      5      6      7      8      9     10
# 119243    160     90     59     42     18     26     25     27     13     15
#     11     12     13     14     15     16     17     18     19     20     21
#     13     12     15     15     14     17     14      7      9      7     12
#     22     23     24     25     26     27     28     29     30     31     32
#      7      4      9     10     11     11      8      9      6      7      1
#     33     34     35     36     37     38     39     40     41     42     43
#      3      7      3      3      6      5      4      1      1      5      7
#     44     45     46     47     49     50     51     52     53     54     56
#      1      5      1     10      1      1     11      3      3      2      2
#     57     58     59     60     61     62     63     64     65     66     67
#      2      2      3      3      4      4      3      4      2      8      3
#     68     69     70     71     73     74     75     76     77     78     79
#      5      1      2      4      2      3      2      2      2      3     17
#     80     81     83     84     85     87     88     89     90     91     92
#      8      1      7      5      4      5      5      2      2      3      7
#     93     94     95     96     97     98     99    100    102    103    105
#      1      3      2      1      1      2      2      1      1      1      1
#    107    108    109    110    111    112    113    114    115    116    117
#      1      2      4      1      3      3     10      5      2      3      5
#    118    119    120    121    122    123    125    126    127    128    132
#      3      4      2      2      6     10      2      1      2      2      1
#    137    139    140    142    146    149    154    155    159    160    162
#      1      1      1      1      1      1      1      6      1     15      1
#    163    178    182    183    201    207    212    213    214    218    223
#      1      1      2      1      2      1      1      2      2      1      1
#    226    227    239    243    260    271    295    367    368    783    784
#      1      1      1      1      1      1      1      1      1      3      5
#    785
#      3
#
## jxn
# [1] "Number of associations per feature involving PGC SNPs"
#
#     0     1     2     3     4     5     6     7     8     9    10    11    12
# 66710    72    36    25    20    16    21    14    12    10     7     7     4
#    13    14    15    16    17    18    19    20    21    22    23    24    25
#     9     5     7    12    10     3     5     4     7     5     4     6     8
#    26    27    28    29    30    31    32    34    35    39    41    42    43
#     7     4     6     2     7     4     1     3     4     5     2     2     2
#    44    45    46    47    48    49    50    51    52    54    55    56    57
#     1     3     1     4     3     2     3     2     1     2     1     2     1
#    60    61    62    63    64    67    70    71    72    73    75    76    78
#     3     4     3     1     1     3     3     2     2     2     1     1     5
#    79    80    82    83    87    89    90    92    93    95    96    97    98
#     7     1     1     3     2     1     1     2     1     1     1     2     2
#    99   104   106   107   111   112   113   114   118   120   121   122   123
#     1     2     1     3     2     3     2     1     2     2     1    10     3
#   124   125   126   130   134   137   147   150   152   153   154   155   159
#     1     1     1     2     1     2     1     1     1     1     2     1     3
#   160   165   170   216   217   226   233   257   379   555   782   783   784
#     1     1     1     1     1     1     1     1     1     1     1     1     4
#
## tx
# [1] "Number of associations per feature involving PGC SNPs"
#
#     0     1     2     3     4     5     6     7     9    10    11    12    13
# 27282    42    29    14    13     7    10     8     1     7     5     1     1
#    14    15    16    18    19    21    22    23    24    25    26    27    28
#     1     1     1     3     3     5     1     4     2     1     3     3     1
#    29    31    32    33    34    35    36    38    39    40    41    42    44
#     2     2     1     1     2     2     1     1     1     1     1     5     2
#    45    46    47    48    49    50    51    52    53    54    56    59    61
#     2     1     3     1     1     1     1     1     1     1     1     1     1
#    62    64    65    66    67    68    69    70    71    72    73    74    79
#     1     1     1     1     1     1     1     2     1     1     3     2     4
#    80    82    84    87    88    89    92    93    98   100   101   106   113
#     1     1     1     1     1     2     1     1     2     1     1     1     2
#   114   115   117   118   119   120   121   122   125   128   151   155   157
#     2     1     1     1     1     1     1     4     1     1     1     1     1
#   158   159   160   167   188   197   199   204   227   779   784
#     1     1     2     1     1     1     1     1     1     1     1

dlpfc_nums_feat
#      sczd_casecontrol_DE sczd_casecontrol_notDE snp_assoc_pgc snp_assoc_no_pgc
# gene 245                 24407                  145           14116
# exon 440                 396143                 1063          119243
# jxn  37                  297144                 527           66710
# tx   6                   92726                  270           27282

hippo_nums_feat <- get_nums_feat(hippo,
    list(
        'gene' = outGene$HIPPO_matchQSV,
        'exon' = outFeat$HIPPO$exon,
        'jxn' = outFeat$HIPPO$jxn,
        'tx' = outFeat$HIPPO$tx
    )
)

## gene
# [1] "Number of associations per feature involving PGC SNPs"
#
#     0     1     2     3     4     5     6     7     8     9    10    11    12
# 11363    14     8     3     4     5     4     1     1     3     2     1     1
#    13    15    16    17    18    19    20    21    22    24    25    29    33
#     2     1     1     1     1     2     1     2     3     1     3     1     1
#    39    43    47    56    57    60    66    68    76    79    83    89    97
#     2     1     1     1     1     1     1     1     1     2     1     1     1
#   100   107   108   113   122   131   150   155   160   200   230   280   300
#     1     1     2     1     2     2     1     1     3     1     1     1     1
#
## exon
# [1] "Number of associations per feature involving PGC SNPs"
#
#     0     1     2     3     4     5     6     7     8     9    10    11    12
# 88784   107    44    35    22    24    28    26    12    19    11    18     3
#    14    15    16    17    18    19    20    21    22    23    24    25    26
#    10    18    14     6    11     5     2     6     8     5     3     5     4
#    27    28    29    30    31    32    33    34    35    36    37    38    39
#    10     4     3     4     6     2     3     6     7     3     2     2     2
#    40    41    43    44    45    46    47    48    51    52    53    54    55
#     4     2     1     1     6     3     2     1     1     3     1     2     2
#    56    59    60    61    63    64    66    68    69    70    71    72    74
#     1     2     2     2     1     3     2     2     1     1     1     1     2
#    75    77    78    79    80    81    82    83    89    91    93    95    97
#     2     1     7     9     1     1     2     3     2     2     2     4     3
#    98    99   100   102   104   105   106   107   108   109   110   111   113
#     1     1     1     3     1     1     6     1     3     1     1     1     1
#   114   115   116   117   119   120   121   122   123   124   131   136   137
#     1     1     2     1     1     2     2     2     1     1     2     1     1
#   138   141   142   144   147   148   149   150   151   152   153   154   155
#     3     1     1     1     2     3     3     2     1     3     1     1     4
#   158   159   160   161   163   166   173   180   189   202   203   209   210
#     2     7     5     1     1     1     2     1     3     1     1     1     1
#   211   219   223   227   240   286   308   314   339   358   371   377
#     1     1     1     1     1     1     1     1     1     1     1     1
#
## jxn
# [1] "Number of associations per feature involving PGC SNPs"
#
#     0     1     2     3     4     5     6     7     8     9    10    11    12
# 51593    64    27    22    15    19    14     8     8     9    10    10     2
#    13    14    15    16    17    18    19    20    21    22    23    24    25
#     5     3    11     8     5     1     2     1     2     6     6     4     6
#    26    27    28    29    30    31    32    33    34    35    36    37    40
#     3     2     5     4     2     3     2     3     7     1     1     1     1
#    41    42    43    44    46    47    48    52    53    55    56    58    59
#     2     1     4     2     2     3     1     3     1     1     1     4     2
#    60    61    62    64    65    69    70    71    73    76    77    78    79
#     1     4     1     3     1     3     2     1     1     1     1     4     5
#    82    83    85    86    87    88    90    91   101   102   103   105   110
#     2     2     1     1     1     1     2     2     2     2     2     2     2
#   111   113   114   116   118   121   122   124   125   128   130   136   141
#     1     2     1     1     1     3     6     1     1     1     1     1     1
#   142   145   147   148   153   154   156   157   160   179   180   193   203
#     1     1     1     1     1     1     1     1     3     1     1     1     1
#   204   208   217   239   248   253   256   276   342   379
#     1     1     1     1     1     1     1     1     1     1
#
## tx
# [1] "Number of associations per feature involving PGC SNPs"
#
#     0     1     2     3     4     5     6     7     8     9    10    11    12
# 21758    28    23    14     9     2     9     3     3     3     8     5     4
#    13    16    17    18    21    22    23    24    26    29    30    31    33
#     1     1     1     1     4     1     2     3     1     1     2     1     2
#    35    36    39    40    41    45    47    48    50    51    52    53    54
#     1     1     3     3     1     1     2     1     1     1     1     2     1
#    56    57    61    62    63    66    69    72    73    74    77    79    81
#     1     1     1     1     1     2     1     2     1     2     1     5     1
#    82    89    93    95    97   101   107   110   112   119   121   122   132
#     1     2     1     2     1     1     2     1     2     1     1     4     2
#   140   144   145   152   155   156   158   159   160   165   211   226   227
#     1     3     1     1     1     1     1     1     2     1     1     1     1
#   235
#     1

hippo_nums_feat
#      sczd_casecontrol_DE sczd_casecontrol_notDE snp_assoc_pgc snp_assoc_no_pgc
# gene 48                  24604                  100           11363
# exon 197                 396386                 694           88784
# jxn  41                  297140                 426           51593
# tx   0                   92732                  209           21758

signif(apply(dlpfc_nums_feat, 1, phyp_qtl), 3)
#  gene  exon   jxn    tx
# 0.687 1.000 1.000 1.000
#signif(1 - apply(dlpfc_nums_feat, 1, phyp_qtl), 3)

signif(apply(hippo_nums_feat, 1, phyp_qtl), 3)
# gene exon  jxn   tx
#    1    1    1    1

## But hm.. in most cases (except DLPFC gene) there are more features with PGC associations that DE genes.
## Here's another version.
phyp_qtl_v2 <- function(x) {
    x <- unlist(x)
    phyper(x['sczd_casecontrol_DE'], x['snp_assoc_pgc'], x['snp_assoc_no_pgc'], x['snp_assoc_no_pgc'] + x['snp_assoc_pgc'])    
}
## Something doesn't seem right with this other one either since not all features have SNP associations.

signif(apply(dlpfc_nums_feat, 1, phyp_qtl_v2), 3)
# gene exon  jxn   tx
   # 1    0    0    0
signif(apply(hippo_nums_feat, 1, phyp_qtl_v2), 3)
# gene exon  jxn   tx
#    0    0    0    0

get_nums_featv2 <- function(qtls, des) {
    mapply(function(qtl, de) {
        
        summ <- qtl[, sum(pgc), keyby = gene]
        
        res <- data.frame(
            'feature_id' = rownames(de),
            't' = de$t,
            'sczd' = de$adj.P.Val < 0.05,
            'pgc' = NA,
            stringsAsFactors = FALSE
        )
        
        m <- match(res$feature_id, summ$gene)
        res$pgc[!is.na(m)] <- summ$V1[m[!is.na(m)]] > 0
        return(res) 
    }, qtls, des, SIMPLIFY = FALSE)
}

dlpfc_nums_featv2 <- get_nums_featv2(dlpfc,
    list(
        'gene' = outGene$DLPFC_matchQSV,
        'exon' = outFeat$DLPFC$exon,
        'jxn' = outFeat$DLPFC$jxn,
        'tx' = outFeat$DLPFC$tx
    )
)

hippo_nums_featv2 <- get_nums_featv2(hippo,
    list(
        'gene' = outGene$HIPPO_matchQSV,
        'exon' = outFeat$HIPPO$exon,
        'jxn' = outFeat$HIPPO$jxn,
        'tx' = outFeat$HIPPO$tx
    )
)

head(dlpfc_nums_featv2$gene)
#          feature_id           t  sczd pgc
# 1 ENSG00000227232.5 -0.38564882 FALSE  NA
# 2 ENSG00000278267.1  0.08617349 FALSE  NA
# 3 ENSG00000269981.1  1.11204088 FALSE  NA
# 4 ENSG00000279457.3 -0.90787551 FALSE  NA
# 5 ENSG00000228463.9 -0.56228106 FALSE  NA
# 6 ENSG00000236679.2  0.70473228 FALSE  NA
table(is.na(dlpfc_nums_featv2$gene$pgc))
# FALSE  TRUE
# 14261 10391

(dlpfc_tab <- lapply(dlpfc_nums_featv2, function(x) {
    table(
        'SCZD DE' = factor(x$sczd, levels = c('FALSE', 'TRUE')),
        'PGC assoc' = factor(x$pgc, levels = c('FALSE', 'TRUE'))
    )
}))
# $gene
#        PGC assoc
# SCZD DE FALSE  TRUE
#   FALSE 13955   144
#   TRUE    161     1
#
# $exon
#        PGC assoc
# SCZD DE  FALSE   TRUE
#   FALSE 119105   1056
#   TRUE     138      7
#
# $jxn
#        PGC assoc
# SCZD DE FALSE  TRUE
#   FALSE 66704   526
#   TRUE      6     1
#
# $tx
#        PGC assoc
# SCZD DE FALSE  TRUE
#   FALSE 27280   269
#   TRUE      2     1

(hippo_tab <- lapply(hippo_nums_featv2, function(x) {
    table(
        'SCZD DE' = factor(x$sczd, levels = c('FALSE', 'TRUE')),
        'PGC assoc' = factor(x$pgc, levels = c('FALSE', 'TRUE'))
    )
}))
# $gene
#        PGC assoc
# SCZD DE FALSE  TRUE
#   FALSE 11340   100
#   TRUE     23     0
#
# $exon
#        PGC assoc
# SCZD DE FALSE  TRUE
#   FALSE 88742   694
#   TRUE     42     0
#
# $jxn
#        PGC assoc
# SCZD DE FALSE  TRUE
#   FALSE 51590   426
#   TRUE      3     0
#
# $tx
#        PGC assoc
# SCZD DE FALSE  TRUE
#   FALSE 21758   209
#   TRUE      0     0

## I think this is the right hypergeometric test:
# * number of white balls draw: DE gene among PGC associated genes: 1
# * number of white balls in urn: DE genes: 161 + 1
# * number of black balls in urn: non DE genes: 13955 + 144
# * number of balls drawn: PGC associated genes: 1 + 144
phyper(1, 161 + 1, 13955 + 144, 144 + 1)
# [1] 0.507626
sapply(dlpfc_tab, function(x) {
    phyper(x[2, 2], sum(x[2, ]), sum(x[1, ]), sum(x[, 2]))
})
#      gene      exon       jxn        tx
# 0.5076260 0.9999503 0.9987455 0.9997148
sapply(hippo_tab, function(x) {
    phyper(x[2, 2], sum(x[2, ]), sum(x[1, ]), sum(x[, 2]))
})
#      gene      exon       jxn        tx
# 0.8173235 0.7210108 0.9756322 1.0000000

lapply(dlpfc_tab, chisq.test)
# $gene
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 0.013433, df = 1, p-value = 0.9077
#
#
# $exon
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 21.474, df = 1, p-value = 3.587e-06
#
#
# $jxn
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 3.6404, df = 1, p-value = 0.05639
#
#
# $tx
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 7.6085, df = 1, p-value = 0.005809

lapply(hippo_tab, chisq.test)
# $gene
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 5.7023e-29, df = 1, p-value = 1
#
#
# $exon
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 3.8966e-26, df = 1, p-value = 1
#
#
# $jxn
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 7.1069e-23, df = 1, p-value = 1
#
#
# $tx
#
#     Pearson's Chi-squared test
#
# data:  X[[i]]
# X-squared = NaN, df = 1, p-value = NA

x <- dlpfc_nums_featv2$gene
library('limma')


geneSetTest(index = x$pgc[!is.na(x$pgc)], x$t[!is.na(x$pgc)])

lapply(dlpfc_nums_featv2, function(x) 


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
