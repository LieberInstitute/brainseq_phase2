# Based on https://github.com/LieberInstitute/brainseq_phase2/blob/master/check_sex/casecontrol/check_sex_casecontrol.R

library('SummarizedExperiment')
library('purrr')
library('jaffelab')
library('sessioninfo')

## For the gene annotation info
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_gene_unfiltered.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_exon_unfiltered.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_jxn_unfiltered.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_tx_unfiltered.Rdata", verbose = TRUE)

gr <- list(
    'gene' = rowRanges(rse_gene),
    'exon' = rowRanges(rse_exon),
    'jxn' = rowRanges(rse_jxn),
    'tx' = rowRanges(rse_tx)
)

gr$jxn$gene_type <- gr$gene$gene_type[match(gr$jxn$newGeneID, gr$gene$gencodeID)]


compute_overlap <- function(feat) {
    table(
        'Expressed' = factor(feat$passExprsCut, levels = c('FALSE', 'TRUE')),
        'Protein Coding' = factor(feat$gene_type == 'protein_coding', levels = c('FALSE', 'TRUE')),
        useNA = 'ifany'
    )
}

map(gr, ~ addmargins(compute_overlap(.x)))
# $gene
#          Protein Coding
# Expressed FALSE  TRUE   Sum
#     FALSE 28088  5297 33385
#     TRUE   9999 14653 24652
#     Sum   38087 19950 58037
#
# $exon
#          Protein Coding
# Expressed  FALSE   TRUE    Sum
#     FALSE  81822  93218 175040
#     TRUE   27355 369228 396583
#     Sum   109177 462446 571623
#
# $jxn
#          Protein Coding
# Expressed  FALSE   TRUE   <NA>    Sum
#     FALSE  36366 258420 245110 539896
#     TRUE   16446 236243  44492 297181
#     Sum    52812 494663 289602 837077
#
# $tx
#          Protein Coding
# Expressed  FALSE   TRUE    Sum
#     FALSE  32034  73327 105361
#     TRUE   20604  72128  92732
#     Sum    52638 145455 198093

make_table <- function(stats) {
    ov <- map(stats, ~ compute_overlap(.x))
    
    res <- map_dfr(ov,
        ~ as.data.frame(matrix(as.vector(.x[1:2, 1:2]), nrow = 1, dimnames = list(1, c('Null_both', 'Expressed_only', 'PC_only', 'Expressed_PC'))))
    )
    res$feature <- names(stats)
    res$OR <- map_dbl(ov, ~ getOR(.x[1:2, 1:2]))
    res$pval <- map_dbl(ov, ~ chisq.test(.x[1:2, 1:2])$p.value)
    res$pval_bonf <- p.adjust(res$pval, 'bonf')
    
    ## Re-order by feature
    # res <- res[c(1,5,2,6,3,7,4,8), ]
    return(res)
}

options(width = 120)
table_expr_pc <- make_table(gr)
table_expr_pc
#   Null_both Expressed_only PC_only Expressed_PC feature        OR pval pval_bonf
# 1     28088           9999    5297        14653    gene  7.770712    0         0
# 2     81822          27355   93218       369228    exon 11.847541    0         0
# 3     36366          16446  258420       236243     jxn  2.021474    0         0
# 4     32034          20604   73327        72128      tx  1.529324    0         0


gene_types <- data.frame(
    sort(table(gr$gene$gene_type), decreasing = TRUE)
)
colnames(gene_types) <- c('gene_type', 'frequency_all')
gene_types$gene_type <- as.character(gene_types$gene_type)
gene_types$percent_all <- gene_types$frequency_all / sum(gene_types$frequency_all) * 100
gene_types$frequency_expressed <- table(gr$gene$gene_type[gr$gene$passExprsCut])[gene_types$gene_type]
gene_types$frequency_expressed[is.na(gene_types$frequency_expressed)] <- 0
gene_types$percent_expressed <- gene_types$frequency_expressed / sum(gene_types$frequency_expressed) * 100
gene_types$ratio_expressed <- gene_types$frequency_expressed / gene_types$frequency_all
options(width = 200)
gene_types
#                             gene_type frequency_all  percent_all frequency_expressed percent_expressed ratio_expressed
# 1                      protein_coding         19950 34.374623085               14653      59.439396398      0.73448622
# 2                processed_pseudogene         10275 17.704223168                1799       7.297582346      0.17508516
# 3                             lincRNA          7539 12.989989145                1965       7.970955703      0.26064465
# 4                           antisense          5530  9.528404294                1618       6.563361999      0.29258590
# 5              unprocessed_pseudogene          2663  4.588452194                 339       1.375141976      0.12730004
# 6                            misc_RNA          2213  3.813084756                 511       2.072854129      0.23090827
# 7                               snRNA          1900  3.273773627                 537       2.178322246      0.28263158
# 8                               miRNA          1569  2.703447801                 728       2.953107253      0.46398980
# 9                                 TEC          1048  1.805744611                 393       1.594191141      0.37500000
# 10                             snoRNA           944  1.626548581                 370       1.500892423      0.39194915
# 11                     sense_intronic           909  1.566242225                 720       2.920655525      0.79207921
# 12 transcribed_unprocessed_pseudogene           737  1.269879560                 302       1.225052734      0.40976934
# 13                               rRNA           544  0.937333081                 112       0.454324193      0.20588235
# 14               processed_transcript           516  0.889087996                 204       0.827519065      0.39534884
# 15   transcribed_processed_pseudogene           450  0.775367438                 162       0.657147493      0.36000000
# 16                    IG_V_pseudogene           193  0.332546479                   2       0.008112932      0.01036269
# 17                  sense_overlapping           187  0.322208246                 108       0.438098329      0.57754011
# 18                 unitary_pseudogene           154  0.265347968                  21       0.085185786      0.13636364
# 19                          IG_V_gene           145  0.249840619                   3       0.012169398      0.02068966
# 20                          TR_V_gene           108  0.186088185                   0       0.000000000      0.00000000
# 21                          TR_J_gene            79  0.136120061                   5       0.020282330      0.06329114
# 22     transcribed_unitary_pseudogene            60  0.103382325                  22       0.089242252      0.36666667
# 23             polymorphic_pseudogene            51  0.087874976                   5       0.020282330      0.09803922
# 24                             scaRNA            49  0.084428899                  19       0.077072854      0.38775510
# 25                          IG_D_gene            37  0.063752434                   0       0.000000000      0.00000000
# 26           3prime_overlapping_ncRNA            30  0.051691163                  11       0.044621126      0.36666667
# 27                    TR_V_pseudogene            30  0.051691163                   2       0.008112932      0.06666667
# 28                            Mt_tRNA            22  0.037906853                  22       0.089242252      1.00000000
# 29                         pseudogene            21  0.036183814                   2       0.008112932      0.09523810
# 30                          IG_J_gene            18  0.031014698                   1       0.004056466      0.05555556
# 31                          IG_C_gene            14  0.024122543                   2       0.008112932      0.14285714
# 32                    IG_C_pseudogene             9  0.015507349                   1       0.004056466      0.11111111
# 33                           ribozyme             8  0.013784310                   2       0.008112932      0.25000000
# 34                          TR_C_gene             6  0.010338233                   4       0.016225864      0.66666667
# 35                               sRNA             5  0.008615194                   1       0.004056466      0.20000000
# 36      bidirectional_promoter_lncRNA             4  0.006892155                   1       0.004056466      0.25000000
# 37                          TR_D_gene             4  0.006892155                   0       0.000000000      0.00000000
# 38                    TR_J_pseudogene             4  0.006892155                   0       0.000000000      0.00000000
# 39                    IG_J_pseudogene             3  0.005169116                   0       0.000000000      0.00000000
# 40                         non_coding             3  0.005169116                   2       0.008112932      0.66666667
# 41                            Mt_rRNA             2  0.003446078                   2       0.008112932      1.00000000
# 42                      IG_pseudogene             1  0.001723039                   0       0.000000000      0.00000000
# 43                       macro_lncRNA             1  0.001723039                   0       0.000000000      0.00000000
# 44                              scRNA             1  0.001723039                   1       0.004056466      1.00000000
# 45                           vaultRNA             1  0.001723039                   0       0.000000000      0.00000000

colSums(gene_types[, -c(1, 6)])
# frequency_all         percent_all frequency_expressed   percent_expressed
#         58037                 100               24652                 100

## List the files for each brain region
f_stub <- c(
    'DLPFC' = 'rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda',
    'HIPPO' = 'rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda'
)
f_ori <- paste0('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/', f_stub)
names(f_ori) <- names(f_stub)
stopifnot(all(file.exists(f_ori)))

load_stats <- function(stats_files) {
    final <- do.call(c, map2(stats_files, names(stats_files), function(f, region) {
        message(paste(Sys.time(), 'loading', f))
        load(f, verbose = TRUE)
        
        stopifnot(identical(rownames(outGene), names(gr$gene[gr$gene$passExprsCut])))
        stopifnot(identical(rownames(outExon), names(gr$exon[gr$exon$passExprsCut])))
        stopifnot(identical(rownames(outJxn), names(gr$jxn[gr$jxn$passExprsCut])))
        stopifnot(identical(rownames(outTx), names(gr$tx[gr$tx$passExprsCut])))
        
        res <- list(
            'Gene' = cbind(outGene, gene_type = gr$gene$gene_type[gr$gene$passExprsCut]),
            'Exon' = cbind(outExon, gene_type = gr$exon$gene_type[gr$exon$passExprsCut]),
            'Jxn' = cbind(outJxn, gene_type = gr$jxn$gene_type[gr$jxn$passExprsCut]),
            'Tx' = cbind(outTx, gene_type = gr$tx$gene_type[gr$tx$passExprsCut])
        )
        names(res) <- paste0(region, '_', names(res))
        
        ## Add Bonferroni
        res <- map(res, function(x) {
            x$p_bonf = p.adjust(x$P.Value, method = "bonf")
        	return(x)
        })
        return(res)
    }))
    names(final) <- gsub('.*\\.', '', names(final))
    return(final)
}

sz_stats <- load_stats(f_ori)


compute_overlap2 <- function(stats) {
    table(
        'SCZD DE' = factor(stats$adj.P.Val < 0.05, levels = c('FALSE', 'TRUE')),
        'Protein Coding' = factor(stats$gene_type == 'protein_coding', levels = c('FALSE', 'TRUE')),
        useNA = 'ifany'
    )
}

map(sz_stats, ~ addmargins(compute_overlap2(.x)))
# $DLPFC_Gene
#        Protein Coding
# SCZD DE FALSE  TRUE   Sum
#   FALSE  9966 14441 24407
#   TRUE     33   212   245
#   Sum    9999 14653 24652
#
# $DLPFC_Exon
#        Protein Coding
# SCZD DE  FALSE   TRUE    Sum
#   FALSE  27340 368803 396143
#   TRUE      15    425    440
#   Sum    27355 369228 396583
#
# $DLPFC_Jxn
#        Protein Coding
# SCZD DE  FALSE   TRUE   <NA>    Sum
#   FALSE  16446 236221  44477 297144
#   TRUE       0     22     15     37
#   Sum    16446 236243  44492 297181
#
# $DLPFC_Tx
#        Protein Coding
# SCZD DE FALSE  TRUE   Sum
#   FALSE 20603 72123 92726
#   TRUE      1     5     6
#   Sum   20604 72128 92732
#
# $HIPPO_Gene
#        Protein Coding
# SCZD DE FALSE  TRUE   Sum
#   FALSE  9992 14612 24604
#   TRUE      7    41    48
#   Sum    9999 14653 24652
#
# $HIPPO_Exon
#        Protein Coding
# SCZD DE  FALSE   TRUE    Sum
#   FALSE  27350 369036 396386
#   TRUE       5    192    197
#   Sum    27355 369228 396583
#
# $HIPPO_Jxn
#        Protein Coding
# SCZD DE  FALSE   TRUE   <NA>    Sum
#   FALSE  16446 236215  44479 297140
#   TRUE       0     28     13     41
#   Sum    16446 236243  44492 297181
#
# $HIPPO_Tx
#        Protein Coding
# SCZD DE FALSE  TRUE   Sum
#   FALSE 20604 72128 92732
#   TRUE      0     0     0
#   Sum   20604 72128 92732

make_table2 <- function(stats) {
    ov <- map(stats, ~ compute_overlap2(.x))
    
    res <- map_dfr(ov,
        ~ as.data.frame(matrix(as.vector(.x[1:2, 1:2]), nrow = 1, dimnames = list(1, c('Null_both', 'SCZD_DE_only', 'PC_only', 'SCZD_DE_PC'))))
    )
    res$region <- ss(names(ov), '_', 1)
    res$feature <- tolower(ss(names(ov), '_', 2))
    res$OR <- map_dbl(ov, ~ getOR(.x[1:2, 1:2]))
    set.seed(20180308)
    res$pval <- map_dbl(ov, ~ chisq.test(.x[1:2, 1:2], simulate.p.value = TRUE, B = 1e5)$p.value)
    res$pval_bonf <- p.adjust(res$pval, 'bonf')
    
    ## Re-order by feature
    res <- res[c(1,5,2,6,3,7,4,8), ]
    return(res)
}

table_sz <- make_table2(sz_stats)
table_sz
#   Null_both SCZD_DE_only PC_only SCZD_DE_PC region feature       OR         pval    pval_bonf
# 1      9966           33   14441        212  DLPFC    gene 4.433488 0.0000099999 0.0000699993
# 5      9992            7   14612         41  HIPPO    gene 4.005240 0.0004699953 0.0032899671
# 2     27340           15  368803        425  DLPFC    exon 2.100399 0.0046699533 0.0326896731
# 6     27350            5  369036        192  HIPPO    exon 2.845901 0.0157898421 0.1105288947
# 3     16446            0  236221         22  DLPFC     jxn      Inf 0.3968560314 1.0000000000
# 7     16446            0  236215         28  HIPPO     jxn      Inf 0.2568474315 1.0000000000
# 4     20603            1   72123          5  DLPFC      tx 1.428324 1.0000000000 1.0000000000
# 8     20604            0   72128          0  HIPPO      tx      NaN          NaN          Na

## Which are the SCZD DE genes that are not protein coding?
map(sz_stats[c(1, 5)], ~ sort(table(.x$gene_type[.x$gene_type != 'protein_coding' & .x$adj.P.Val < 0.05])))
# $DLPFC_Gene
#
#                              miRNA                           misc_RNA
#                                  1                                  1
#                             scaRNA             unprocessed_pseudogene
#                                  1                                  1
#                                TEC               processed_pseudogene
#                                  2                                  3
# transcribed_unprocessed_pseudogene                          antisense
#                                  4                                  9
#                            lincRNA
#                                 11
#
# $HIPPO_Gene
#
#                  miRNA   processed_transcript         sense_intronic      sense_overlapping
#                      1                      1                      1                      1
# unprocessed_pseudogene              antisense
#                      1                      2

options(width = 300)
map(sz_stats[c(1, 5)], ~ .x[.x$gene_type != 'protein_coding' & .x$adj.P.Val < 0.05, ])
# $DLPFC_Gene
#                    Length          gencodeID       ensemblID                          gene_type      Symbol  EntrezID Class   meanExprs NumTx    gencodeTx passExprsCut       logFC     AveExpr         t      P.Value   adj.P.Val            B                          gene_type     p_bonf
# ENSG00000217801.9    2171  ENSG00000217801.9 ENSG00000217801 transcribed_unprocessed_pseudogene                    NA InGen   0.6467002     6 ENST0000....         TRUE  0.18364495  1.05634480  4.269374 2.517046e-05 0.011932733  2.316607094 transcribed_unprocessed_pseudogene 0.62050213
# ENSG00000249087.6    2655  ENSG00000249087.6 ENSG00000249087                          antisense  ZNF436-AS1        NA InGen   1.8058076     6 ENST0000....         TRUE  0.09097166  2.56181311  3.587634 3.803550e-04 0.042815123 -0.105286557                          antisense 1.00000000
# ENSG00000275131.2    9026  ENSG00000275131.2 ENSG00000275131 transcribed_unprocessed_pseudogene                    NA InGen   8.6455383     2 ENST0000....         TRUE  0.08229392  6.43989901  3.994339 7.884084e-05 0.022099458  0.888014379 transcribed_unprocessed_pseudogene 1.00000000
# ENSG00000226067.6    2075  ENSG00000226067.6 ENSG00000226067                            lincRNA   LINC00623        NA InGen   1.2757396     5 ENST0000....         TRUE  0.15359368  1.84952938  3.917310 1.073445e-04 0.026334409  1.069540504                            lincRNA 1.00000000
# ENSG00000231752.5    4333  ENSG00000231752.5 ENSG00000231752 transcribed_unprocessed_pseudogene       EMBP1    647121 InGen   0.5348729     3 ENST0000....         TRUE -0.10872351  1.44381892 -3.644484 3.077971e-04 0.040116920  0.149106001 transcribed_unprocessed_pseudogene 1.00000000
# ENSG00000232973.11   7447 ENSG00000232973.11 ENSG00000232973                          antisense  CYP1B1-AS1    285154 InGen   0.5298032    35 ENST0000....         TRUE  0.16494630  2.53270968  3.729530 2.231258e-04 0.035265431  0.374246268                          antisense 1.00000000
# ENSG00000252010.1     276  ENSG00000252010.1 ENSG00000252010                             scaRNA     SCARNA5    677775 InGen 577.3816840     1 ENST0000....         TRUE -0.15138377  7.46049234 -4.161513 3.967920e-05 0.017160905  1.466610579                             scaRNA 0.97817157
# ENSG00000241684.5    7232  ENSG00000241684.5 ENSG00000241684                          antisense ADAMTS9-AS2        NA InGen   0.9304402     4 ENST0000....         TRUE  0.18055897  3.19925243  5.044884 7.246952e-07 0.001988977  5.601290624                          antisense 0.01786519
# ENSG00000277118.1     286  ENSG00000277118.1 ENSG00000277118                           misc_RNA                    NA InGen   1.8601384     1 ENST0000....         TRUE -0.36607651 -0.58423708 -4.832425 2.008307e-06 0.002912282  3.852797890                           misc_RNA 0.04950879
# ENSG00000250476.1     627  ENSG00000250476.1 ENSG00000250476               processed_pseudogene     ENPP7P9        NA InGen   0.3874621     1 ENST0000....         TRUE  0.31884104 -1.75868566  3.577266 3.952095e-04 0.044057714 -0.609684337               processed_pseudogene 1.00000000
# ENSG00000152931.7    9752  ENSG00000152931.7 ENSG00000152931                            lincRNA       PART1     25859 InGen   3.6312970     3 ENST0000....         TRUE -0.09983194  5.98261356 -3.701476 2.482722e-04 0.036430984 -0.138481891                            lincRNA 1.00000000
# ENSG00000248909.1     581  ENSG00000248909.1 ENSG00000248909               processed_pseudogene    HMGB1P21        NA InGen   1.0288115     1 ENST0000....         TRUE -0.21492348 -0.83667300 -3.876015 1.264022e-04 0.028072666  0.549837192               processed_pseudogene 1.00000000
# ENSG00000220563.1     715  ENSG00000220563.1 ENSG00000220563               processed_pseudogene       PKMP3        NA InGen   5.5245251     1 ENST0000....         TRUE  0.11218601  2.20974878  3.617944 3.398804e-04 0.040367706  0.025308125               processed_pseudogene 1.00000000
# ENSG00000223414.2    4376  ENSG00000223414.2 ENSG00000223414                            lincRNA   LINC00473     90632 InGen   2.0767367     5 ENST0000....         TRUE -0.20041607  3.73342100 -3.624044 3.322388e-04 0.040367706 -0.129728544                            lincRNA 1.00000000
# ENSG00000227544.8    2849  ENSG00000227544.8 ENSG00000227544                            lincRNA                    NA InGen   3.2219357     5 ENST0000....         TRUE -0.22548045  3.09792826 -3.518841 4.896089e-04 0.049466550 -0.397681666                            lincRNA 1.00000000
# ENSG00000251521.2     977  ENSG00000251521.2 ENSG00000251521 transcribed_unprocessed_pseudogene     IMPA1P1        NA InGen   0.2569406     2 ENST0000....         TRUE  0.36890038 -2.15215596  3.873201 1.278108e-04 0.028132061 -0.002190503 transcribed_unprocessed_pseudogene 1.00000000
# ENSG00000279670.1    1670  ENSG00000279670.1 ENSG00000279670                                TEC                    NA InGen   0.5485052     1 ENST0000....         TRUE -0.17277719  0.06890644 -3.609957 3.501340e-04 0.040714640 -0.037577285                                TEC 1.00000000
# ENSG00000229065.1     529  ENSG00000229065.1 ENSG00000229065                          antisense                    NA InGen   0.9368437     1 ENST0000....         TRUE -0.23391000 -0.84880121 -3.599599 3.638652e-04 0.041720956 -0.259983816                          antisense 1.00000000
# ENSG00000228914.2     937  ENSG00000228914.2 ENSG00000228914             unprocessed_pseudogene      OR1H1P        NA InGen   0.2935881     1 ENST0000....         TRUE -0.42820820 -1.50632266 -3.942395 9.713303e-05 0.024685808  0.455140511             unprocessed_pseudogene 1.00000000
# ENSG00000234156.1     636  ENSG00000234156.1 ENSG00000234156                          antisense                    NA InGen   0.2841293     2 ENST0000....         TRUE -0.56689784 -2.09413456 -4.845357 1.889398e-06 0.002911090  2.777480431                          antisense 0.04657744
# ENSG00000278266.1    2933  ENSG00000278266.1 ENSG00000278266                            lincRNA                    NA InGen   0.7068402     1 ENST0000....         TRUE -0.17667453  1.34420235 -3.996720 7.808630e-05 0.022099458  1.344036866                            lincRNA 1.00000000
# ENSG00000256064.1     415  ENSG00000256064.1 ENSG00000256064                            lincRNA                    NA InGen   2.2728502     1 ENST0000....         TRUE  0.17632293 -0.19431469  3.663540 2.865416e-04 0.039026647  0.092485504                            lincRNA 1.00000000
# ENSG00000279146.1    1308  ENSG00000279146.1 ENSG00000279146                                TEC                    NA InGen   1.6114912     1 ENST0000....         TRUE  0.14456658  0.94053622  4.108875 4.937751e-05 0.019289635  1.723417176                                TEC 1.00000000
# ENSG00000204603.6    2562  ENSG00000204603.6 ENSG00000204603                            lincRNA   LINC01257    116437 InGen   0.3197739     2 ENST0000....         TRUE -0.29131971  0.35557030 -4.000668 7.685014e-05 0.022099458  1.250872117                            lincRNA 1.00000000
# ENSG00000226519.1     340  ENSG00000226519.1 ENSG00000226519                            lincRNA   LINC00390        NA InGen   2.5933369     1 ENST0000....         TRUE  0.29485012  0.26748864  4.899043 1.464450e-06 0.002729704  4.600821125                            lincRNA 0.03610161
# ENSG00000199088.5      67  ENSG00000199088.5 ENSG00000199088                              miRNA      MIR379    494328 InGen   1.2355836     1 ENST0000....         TRUE  0.46195520 -3.64864556  3.751286 2.052988e-04 0.034278814 -1.017444495                              miRNA 1.00000000
# ENSG00000256802.2    2736  ENSG00000256802.2 ENSG00000256802                          antisense             100130111 InGen   0.6040155     2 ENST0000....         TRUE  0.16494409  1.00955664  4.275368 2.453492e-05 0.011859508  2.332137421                          antisense 0.60483489
# ENSG00000259673.5    6732  ENSG00000259673.5 ENSG00000259673                            lincRNA    IQCH-AS1 100506686 InGen   0.6180747     5 ENST0000....         TRUE -0.08459616  2.40691054 -3.629616 3.254007e-04 0.040116920  0.052872923                            lincRNA 1.00000000
# ENSG00000261997.1    3828  ENSG00000261997.1 ENSG00000261997                            lincRNA                    NA InGen   0.3365153     1 ENST0000....         TRUE -0.26288434  0.18513818 -3.869047 1.299181e-04 0.028342836  0.793514675                            lincRNA 1.00000000
# ENSG00000277825.1     721  ENSG00000277825.1 ENSG00000277825                            lincRNA                    NA InGen   0.2889462     1 ENST0000....         TRUE  0.33137994 -2.37634158  3.672320 2.772195e-04 0.038288811 -0.612338656                            lincRNA 1.00000000
# ENSG00000268001.1    1969  ENSG00000268001.1 ENSG00000268001                          antisense   CARD8-AS1 100505812 InGen   0.2970611     1 ENST0000....         TRUE -0.24802535 -0.85613493 -3.783806 1.811393e-04 0.032358309  0.261179161                          antisense 1.00000000
# ENSG00000225106.1    1354  ENSG00000225106.1 ENSG00000225106                          antisense                    NA InGen   0.2976448     2 ENST0000....         TRUE  0.30300036 -0.76068425  3.548641 4.390945e-04 0.046457324 -0.354605762                          antisense 1.00000000
# ENSG00000234161.1     426  ENSG00000234161.1 ENSG00000234161                          antisense  PABPC5-AS1 102724167 InGen   2.4388142     1 ENST0000....         TRUE  0.16516086 -0.07136440  3.986108 8.150338e-05 0.022099458  1.129525831                          antisense 1.00000000
#
# $HIPPO_Gene
#                   Length         gencodeID       ensemblID              gene_type    Symbol  EntrezID Class meanExprs NumTx    gencodeTx passExprsCut      logFC     AveExpr        t      P.Value  adj.P.Val            B              gene_type    p_bonf
# ENSG00000260526.1   2036 ENSG00000260526.1 ENSG00000260526              antisense           105377371 InGen 0.7053111     1 ENST0000....         TRUE 0.18476967 -0.08079761 4.071596 5.940253e-05 0.03960635  1.444816763              antisense 1.0000000
# ENSG00000221430.1    142 ENSG00000221430.1 ENSG00000221430                  miRNA   MIR1294 100302181 InGen 0.2991786     1 ENST0000....         TRUE 0.58598408 -4.58597282 4.360758 1.768923e-05 0.02295131 -0.005430596                  miRNA 0.4360748
# ENSG00000253671.1   1053 ENSG00000253671.1 ENSG00000253671   processed_transcript                  NA InGen 1.7397323     2 ENST0000....         TRUE 0.20616894  0.70477982 4.286179 2.433105e-05 0.02856233  2.329098793   processed_transcript 0.5998090
# ENSG00000197882.3    910 ENSG00000197882.3 ENSG00000197882 unprocessed_pseudogene   OR7E13P        NA InGen 0.5256227     1 ENST0000....         TRUE 0.35803343 -2.05526096 3.974474 8.789916e-05 0.04616611  0.496793379 unprocessed_pseudogene 1.0000000
# ENSG00000261799.1   3170 ENSG00000261799.1 ENSG00000261799      sense_overlapping                  NA InGen 3.7991942     1 ENST0000....         TRUE 0.09373977  3.07627516 4.372106 1.684511e-05 0.02295131  2.685540719      sense_overlapping 0.4152658
# ENSG00000270048.1    263 ENSG00000270048.1 ENSG00000270048         sense_intronic                  NA InGen 6.3770150     1 ENST0000....         TRUE 0.16038547  0.44721958 4.380628 1.623676e-05 0.02295131  2.639757779         sense_intronic 0.4002686
# ENSG00000260992.1    807 ENSG00000260992.1 ENSG00000260992              antisense DOCK9-AS2        NA InGen 0.8688187     1 ENST0000....         TRUE 0.27398946 -0.91054465 4.202414 3.462431e-05 0.03282917  1.659505846              antisense 0.8535584



# export_table <- function(tab, filename) {
    write.csv(tab, file = paste0(filename, '.csv'), row.names = FALSE, quote = FALSE)
}
export_table(table_expr_pc, 'expressed_vs_protein_coding')
export_table(table_sz, 'sz_vs_protein_coding')

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
#  date     2019-03-11
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  bindr                  0.1.1     2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 [1] CRAN (R 3.5.0)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.5    2019-01-04 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                  0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
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
#  limma                  3.38.3    2018-12-02 [1] Bioconductor
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.0)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                * 0.2.5     2018-05-29 [2] CRAN (R 3.5.0)
#  R6                     2.3.0     2018-10-04 [2] CRAN (R 3.5.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr                  0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
