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

## Match the rest of the names used for the other browser files
names(deres) <- c("development", "regionspecific", "sczd_casecontrol")

## Explore the columns for each model
lapply(deres, function(x) head(x$gene))
# $development
#                        Age.RegionHIPPO RegionHIPPO.fetal RegionHIPPO.birth
# gene.ENSG00000227232.5       0.2697299        0.15581216         0.1322596
# gene.ENSG00000278267.1       1.9601751        1.26193327        -1.4895673
# gene.ENSG00000269981.1       4.2690953        2.80205028        -2.8921325
# gene.ENSG00000279457.3       2.8975233        1.43030974        -3.1954019
# gene.ENSG00000228463.9       0.3618521        0.05269397        -1.2220409
# gene.ENSG00000236679.2       4.4252513        1.33663350        -5.1996594
#                        RegionHIPPO.infant RegionHIPPO.child RegionHIPPO.teen
# gene.ENSG00000227232.5         -0.4566381        0.05835068      0.003025649
# gene.ENSG00000278267.1         -0.4741538       -0.04049151      0.047792118
# gene.ENSG00000269981.1         -1.3697974       -0.04027343      0.027785197
# gene.ENSG00000279457.3          0.3590451       -0.12967817      0.086913714
# gene.ENSG00000228463.9          0.9067541       -0.06455060      0.019558403
# gene.ENSG00000236679.2          0.7982818       -0.01164576     -0.019399668
#                        RegionHIPPO.adult    AveExpr          F       P.Value
# gene.ENSG00000227232.5      -0.016199796  1.0287983  376.45464 1.292818e-213
# gene.ENSG00000278267.1      -0.003501927 -1.9897147  131.47785 2.722490e-116
# gene.ENSG00000269981.1       0.001060586 -3.3995867 1433.76142  0.000000e+00
# gene.ENSG00000279457.3      -0.021362841  1.5311529   23.99013  8.411078e-29
# gene.ENSG00000228463.9      -0.002402971  1.3012205   80.80571  4.092264e-82
# gene.ENSG00000236679.2       0.014307537 -0.6546682   26.96296  3.513726e-32
#                            adj.P.Val type        P.Bonf span_Age.RegionHIPPO
# gene.ENSG00000227232.5 5.855326e-213 gene 3.187054e-209            8.0021738
# gene.ENSG00000278267.1 8.087096e-116 gene 6.711481e-112           -0.3196381
# gene.ENSG00000269981.1  0.000000e+00 gene  0.000000e+00            7.0747345
# gene.ENSG00000279457.3  1.583670e-28 gene  2.073499e-24           -0.6987235
# gene.ENSG00000228463.9  1.040884e-81 gene  1.008825e-77            7.1261466
# gene.ENSG00000236679.2  6.780460e-32 gene  8.662038e-28           -6.7109491
#                        span_RegionHIPPO.fetal span_RegionHIPPO.birth
# gene.ENSG00000227232.5              4.8358656             -7.6366303
# gene.ENSG00000278267.1             -0.9955068             -0.9892802
# gene.ENSG00000269981.1              2.9981554             -8.1490961
# gene.ENSG00000279457.3             -0.5870097              0.6102499
# gene.ENSG00000228463.9              2.8389013             -8.1404169
# gene.ENSG00000236679.2             -2.9819266              5.0293539
#                        span_RegionHIPPO.infant span_RegionHIPPO.child
# gene.ENSG00000227232.5             -0.34271082              0.3053811
# gene.ENSG00000278267.1              1.55054669             -0.3984407
# gene.ENSG00000269981.1              1.16419566             -0.0592514
# gene.ENSG00000279457.3              0.01141038              0.2138010
# gene.ENSG00000228463.9              0.94316799              0.1536882
# gene.ENSG00000236679.2              2.00529792             -0.4493747
#                        span_RegionHIPPO.teen span_AveExpr   span_F span_P.Value
# gene.ENSG00000227232.5           -0.48973287   -3.0590940 4.880797 3.885935e-04
# gene.ENSG00000278267.1            0.21143315   -4.5325174 7.143530 8.320897e-06
# gene.ENSG00000269981.1           -0.04672582   -4.5919306 3.312900 6.848164e-03
# gene.ENSG00000279457.3           -0.15968076    0.2943755 1.505577 1.913726e-01
# gene.ENSG00000228463.9           -0.08851521   -1.4835077 3.196008 8.521479e-03
# gene.ENSG00000236679.2            0.11331415   -2.3098149 6.320305 3.224673e-05
#                        span_adj.P.Val span_type span_P.Bonf seqnames  start
# gene.ENSG00000227232.5   1.397463e-03      gene   1.0000000     chr1  14404
# gene.ENSG00000278267.1   4.274365e-05      gene   0.2051268     chr1  17369
# gene.ENSG00000269981.1   1.665558e-02      gene   1.0000000     chr1 137682
# gene.ENSG00000279457.3   2.713203e-01      gene   1.0000000     chr1 184923
# gene.ENSG00000228463.9   2.005073e-02      gene   1.0000000     chr1 257864
# gene.ENSG00000236679.2   1.468587e-04      gene   0.7949464     chr1 347982
#                           end width strand Length         gencodeID
# gene.ENSG00000227232.5  29570 15167      -   1351 ENSG00000227232.5
# gene.ENSG00000278267.1  17436    68      -     68 ENSG00000278267.1
# gene.ENSG00000269981.1 137965   284      -    284 ENSG00000269981.1
# gene.ENSG00000279457.3 200322 15400      -   1982 ENSG00000279457.3
# gene.ENSG00000228463.9 297502 39639      -   4039 ENSG00000228463.9
# gene.ENSG00000236679.2 348366   385      -    385 ENSG00000236679.2
#                              ensemblID              gene_type    Symbol
# gene.ENSG00000227232.5 ENSG00000227232 unprocessed_pseudogene    WASH7P
# gene.ENSG00000278267.1 ENSG00000278267                  miRNA MIR6859-1
# gene.ENSG00000269981.1 ENSG00000269981   processed_pseudogene
# gene.ENSG00000279457.3 ENSG00000279457         protein_coding
# gene.ENSG00000228463.9 ENSG00000228463                lincRNA
# gene.ENSG00000236679.2 ENSG00000236679   processed_pseudogene RPL23AP24
#                         EntrezID Class meanExprs NumTx
# gene.ENSG00000227232.5        NA InGen 1.6969802     1
# gene.ENSG00000278267.1 102466751 InGen 4.3552350     1
# gene.ENSG00000269981.1        NA InGen 0.4905604     1
# gene.ENSG00000279457.3 102723897 InGen 1.5866593     3
# gene.ENSG00000228463.9        NA InGen 0.7330683     5
# gene.ENSG00000236679.2        NA InGen 1.9370674     1
#                                                                                                                         gencodeTx
# gene.ENSG00000227232.5                                                                                          ENST00000488147.1
# gene.ENSG00000278267.1                                                                                          ENST00000619216.1
# gene.ENSG00000269981.1                                                                                          ENST00000595919.1
# gene.ENSG00000279457.3                                           c("ENST00000623834.3"; "ENST00000623083.3"; "ENST00000624735.1")
# gene.ENSG00000228463.9 c("ENST00000442116.1"; "ENST00000448958.1"; "ENST00000634344.1"; "ENST00000424587.6"; "ENST00000335577.4")
# gene.ENSG00000236679.2                                                                                          ENST00000458203.2
#                        passExprsCut replicates_in_BrainSpan
# gene.ENSG00000227232.5         TRUE                    TRUE
# gene.ENSG00000278267.1         TRUE                    TRUE
# gene.ENSG00000269981.1         TRUE                    TRUE
# gene.ENSG00000279457.3         TRUE                   FALSE
# gene.ENSG00000228463.9         TRUE                    TRUE
# gene.ENSG00000236679.2         TRUE                    TRUE
#
# $regionspecific
#                                    logFC    AveExpr          t      P.Value
# adult_gene.ENSG00000227232.5 -0.21909569  0.9292268 -2.7371981 6.440476e-03
# adult_gene.ENSG00000278267.1 -0.02520089 -2.0587727 -0.1611352 8.720588e-01
# adult_gene.ENSG00000269981.1  0.07622330 -3.2875290  0.4109298 6.813184e-01
# adult_gene.ENSG00000279457.3 -0.16417396  1.4547027 -2.2039861 2.802793e-02
# adult_gene.ENSG00000228463.9 -0.38693129  1.4140680 -4.1036979 4.822949e-05
# adult_gene.ENSG00000236679.2 -0.03582149 -0.6486502 -0.3492594 7.270569e-01
#                                 adj.P.Val          B   age type P.Bonf
# adult_gene.ENSG00000227232.5 0.0112205381 -3.8909690 adult gene      1
# adult_gene.ENSG00000278267.1 0.8979572516 -7.0932027 adult gene      1
# adult_gene.ENSG00000269981.1 0.7327479344 -6.8540451 adult gene      1
# adult_gene.ENSG00000279457.3 0.0433547361 -5.2812551 adult gene      1
# adult_gene.ENSG00000228463.9 0.0001165221  0.5859877 adult gene      1
# adult_gene.ENSG00000236679.2 0.7741871316 -7.3002362 adult gene      1
#                               span_logFC span_AveExpr     span_t span_P.Value
# adult_gene.ENSG00000227232.5  0.48701908  -2.87280639  0.3855199    0.7102749
# adult_gene.ENSG00000278267.1 -0.16781978  -4.52479062 -0.3009702    0.7713967
# adult_gene.ENSG00000269981.1 -0.55551977  -3.91233552 -1.2149204    0.2603047
# adult_gene.ENSG00000279457.3 -0.10717824  -0.02337325 -0.3096328    0.7650463
# adult_gene.ENSG00000228463.9  0.12622082  -0.51840652  0.2286770    0.8250709
# adult_gene.ENSG00000236679.2  0.05113398  -2.25119522  0.1042499    0.9196329
#                              span_adj.P.Val    span_B span_age span_type
# adult_gene.ENSG00000227232.5      0.9741138 -5.688121    adult      gene
# adult_gene.ENSG00000278267.1      0.9955806 -5.704816    adult      gene
# adult_gene.ENSG00000269981.1      0.6960457 -5.066631    adult      gene
# adult_gene.ENSG00000279457.3      0.9929675 -6.035248    adult      gene
# adult_gene.ENSG00000228463.9      0.9976167 -5.927596    adult      gene
# adult_gene.ENSG00000236679.2      0.9976167 -5.750355    adult      gene
#                              span_P.Bonf seqnames  start    end width strand
# adult_gene.ENSG00000227232.5           1     chr1  14404  29570 15167      -
# adult_gene.ENSG00000278267.1           1     chr1  17369  17436    68      -
# adult_gene.ENSG00000269981.1           1     chr1 137682 137965   284      -
# adult_gene.ENSG00000279457.3           1     chr1 184923 200322 15400      -
# adult_gene.ENSG00000228463.9           1     chr1 257864 297502 39639      -
# adult_gene.ENSG00000236679.2           1     chr1 347982 348366   385      -
#                              Length         gencodeID       ensemblID
# adult_gene.ENSG00000227232.5   1351 ENSG00000227232.5 ENSG00000227232
# adult_gene.ENSG00000278267.1     68 ENSG00000278267.1 ENSG00000278267
# adult_gene.ENSG00000269981.1    284 ENSG00000269981.1 ENSG00000269981
# adult_gene.ENSG00000279457.3   1982 ENSG00000279457.3 ENSG00000279457
# adult_gene.ENSG00000228463.9   4039 ENSG00000228463.9 ENSG00000228463
# adult_gene.ENSG00000236679.2    385 ENSG00000236679.2 ENSG00000236679
#                                           gene_type    Symbol  EntrezID Class
# adult_gene.ENSG00000227232.5 unprocessed_pseudogene    WASH7P        NA InGen
# adult_gene.ENSG00000278267.1                  miRNA MIR6859-1 102466751 InGen
# adult_gene.ENSG00000269981.1   processed_pseudogene                  NA InGen
# adult_gene.ENSG00000279457.3         protein_coding           102723897 InGen
# adult_gene.ENSG00000228463.9                lincRNA                  NA InGen
# adult_gene.ENSG00000236679.2   processed_pseudogene RPL23AP24        NA InGen
#                              meanExprs NumTx
# adult_gene.ENSG00000227232.5 1.6969802     1
# adult_gene.ENSG00000278267.1 4.3552350     1
# adult_gene.ENSG00000269981.1 0.4905604     1
# adult_gene.ENSG00000279457.3 1.5866593     3
# adult_gene.ENSG00000228463.9 0.7330683     5
# adult_gene.ENSG00000236679.2 1.9370674     1
#                                                                                                                               gencodeTx
# adult_gene.ENSG00000227232.5                                                                                          ENST00000488147.1
# adult_gene.ENSG00000278267.1                                                                                          ENST00000619216.1
# adult_gene.ENSG00000269981.1                                                                                          ENST00000595919.1
# adult_gene.ENSG00000279457.3                                           c("ENST00000623834.3"; "ENST00000623083.3"; "ENST00000624735.1")
# adult_gene.ENSG00000228463.9 c("ENST00000442116.1"; "ENST00000448958.1"; "ENST00000634344.1"; "ENST00000424587.6"; "ENST00000335577.4")
# adult_gene.ENSG00000236679.2                                                                                          ENST00000458203.2
#                              passExprsCut replicates_in_BrainSpan
# adult_gene.ENSG00000227232.5         TRUE                   FALSE
# adult_gene.ENSG00000278267.1         TRUE                   FALSE
# adult_gene.ENSG00000269981.1         TRUE                   FALSE
# adult_gene.ENSG00000279457.3         TRUE                   FALSE
# adult_gene.ENSG00000228463.9         TRUE                   FALSE
# adult_gene.ENSG00000236679.2         TRUE                   FALSE
#
# $sczd_casecontrol
#                   Length         gencodeID       ensemblID
# ENSG00000227232.5   1351 ENSG00000227232.5 ENSG00000227232
# ENSG00000278267.1     68 ENSG00000278267.1 ENSG00000278267
# ENSG00000269981.1    284 ENSG00000269981.1 ENSG00000269981
# ENSG00000279457.3   1982 ENSG00000279457.3 ENSG00000279457
# ENSG00000228463.9   4039 ENSG00000228463.9 ENSG00000228463
# ENSG00000236679.2    385 ENSG00000236679.2 ENSG00000236679
#                                gene_type    Symbol  EntrezID Class meanExprs
# ENSG00000227232.5 unprocessed_pseudogene    WASH7P        NA InGen 1.6969802
# ENSG00000278267.1                  miRNA MIR6859-1 102466751 InGen 4.3552350
# ENSG00000269981.1   processed_pseudogene                  NA InGen 0.4905604
# ENSG00000279457.3         protein_coding           102723897 InGen 1.5866593
# ENSG00000228463.9                lincRNA                  NA InGen 0.7330683
# ENSG00000236679.2   processed_pseudogene RPL23AP24        NA InGen 1.9370674
#                   NumTx    gencodeTx passExprsCut        logFC    AveExpr
# ENSG00000227232.5     1 ENST0000....         TRUE -0.021942394  1.3455791
# ENSG00000278267.1     1 ENST0000....         TRUE  0.007710115 -1.6794227
# ENSG00000269981.1     1 ENST0000....         TRUE  0.142483641 -2.8225824
# ENSG00000279457.3     3 ENST0000....         TRUE -0.045899219  1.8750810
# ENSG00000228463.9     5 ENST0000....         TRUE -0.040515525  1.8814260
# ENSG00000236679.2     1 ENST0000....         TRUE  0.046356081 -0.2071824
#                             t   P.Value adj.P.Val         B region type
# ENSG00000227232.5 -0.38564882 0.6999870 0.9100562 -5.820263  DLPFC gene
# ENSG00000278267.1  0.08617349 0.9313769 0.9819371 -5.287980  DLPFC gene
# ENSG00000269981.1  1.11204088 0.2668707 0.6613949 -4.726874  DLPFC gene
# ENSG00000279457.3 -0.90787551 0.3645579 0.7387243 -5.618657  DLPFC gene
# ENSG00000228463.9 -0.56228106 0.5742783 0.8548830 -5.858495  DLPFC gene
# ENSG00000236679.2  0.70473228 0.4814374 0.8086745 -5.358103  DLPFC gene


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
# 2020-05-04 22:02:06 processing development at the gene level
# 2020-05-04 22:02:08 processing development at the exon level
# 2020-05-04 22:02:34 processing development at the jxn level
# 2020-05-04 22:02:56 processing development at the tx level
# 2020-05-04 22:03:03 processing regionspecific at the gene level
# 2020-05-04 22:03:05 processing regionspecific at the exon level
# 2020-05-04 22:03:42 processing regionspecific at the jxn level
# 2020-05-04 22:04:14 processing regionspecific at the tx level
# 2020-05-04 22:04:24 processing sczd_casecontrol at the gene level
# 2020-05-04 22:04:28 processing sczd_casecontrol at the exon level
# 2020-05-04 22:05:01 processing sczd_casecontrol at the jxn level
# 2020-05-04 22:05:34 processing sczd_casecontrol at the tx level

system("chmod 770 BrainSeqPhaseII_stats_*")
system("ls -lh BrainSeqPhaseII_stats_*")
# -rwxrwx--- 1 lcollado lieber_jaffe 239M May  4 22:02 BrainSeqPhaseII_stats_development_exon.csv
# -rwxrwx--- 1 lcollado lieber_jaffe  18M May  4 22:02 BrainSeqPhaseII_stats_development_gene.csv
# -rwxrwx--- 1 lcollado lieber_jaffe 198M May  4 22:02 BrainSeqPhaseII_stats_development_jxn.csv
# -rwxrwx--- 1 lcollado lieber_jaffe  64M May  4 22:03 BrainSeqPhaseII_stats_development_tx.csv
# -rwxrwx--- 1 lcollado lieber_jaffe 363M May  4 22:03 BrainSeqPhaseII_stats_regionspecific_exon.csv
# -rwxrwx--- 1 lcollado lieber_jaffe  27M May  4 22:03 BrainSeqPhaseII_stats_regionspecific_gene.csv
# -rwxrwx--- 1 lcollado lieber_jaffe 304M May  4 22:04 BrainSeqPhaseII_stats_regionspecific_jxn.csv
# -rwxrwx--- 1 lcollado lieber_jaffe 101M May  4 22:04 BrainSeqPhaseII_stats_regionspecific_tx.csv
# -rwxrwx--- 1 lcollado lieber_jaffe 219M May  4 22:05 BrainSeqPhaseII_stats_sczd_casecontrol_exon.csv
# -rwxrwx--- 1 lcollado lieber_jaffe  18M May  4 22:04 BrainSeqPhaseII_stats_sczd_casecontrol_gene.csv
# -rwxrwx--- 1 lcollado lieber_jaffe 197M May  4 22:05 BrainSeqPhaseII_stats_sczd_casecontrol_jxn.csv
# -rwxrwx--- 1 lcollado lieber_jaffe  71M May  4 22:05 BrainSeqPhaseII_stats_sczd_casecontrol_tx.csv
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
#  version  R version 3.6.1 Patched (2019-10-31 r77350)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2020-05-04
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version  date       lib source
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.6    2020-04-05 [1] CRAN (R 3.6.1)
#  Biobase              * 2.46.0   2019-10-29 [2] Bioconductor
#  BiocGenerics         * 0.32.0   2019-10-29 [1] Bioconductor
#  BiocParallel         * 1.20.1   2019-12-21 [1] Bioconductor
#  bitops                 1.0-6    2013-08-17 [2] CRAN (R 3.6.1)
#  cli                    2.0.2    2020-02-28 [1] CRAN (R 3.6.1)
#  codetools              0.2-16   2018-12-24 [3] CRAN (R 3.6.1)
#  colorout             * 1.2-2    2019-10-31 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1    2019-03-18 [2] CRAN (R 3.6.1)
#  crayon                 1.3.4    2017-09-16 [1] CRAN (R 3.6.1)
#  DelayedArray         * 0.12.2   2020-01-06 [2] Bioconductor
#  digest                 0.6.25   2020-02-23 [1] CRAN (R 3.6.1)
#  dplyr                  0.8.5    2020-03-07 [1] CRAN (R 3.6.1)
#  ellipsis               0.3.0    2019-09-20 [1] CRAN (R 3.6.1)
#  fansi                  0.4.1    2020-01-08 [1] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.22.1   2020-03-27 [1] Bioconductor
#  GenomeInfoDbData       1.2.2    2019-10-28 [2] Bioconductor
#  GenomicRanges        * 1.38.0   2019-10-29 [1] Bioconductor
#  ggplot2                3.3.0    2020-03-05 [1] CRAN (R 3.6.1)
#  glue                   1.4.0    2020-04-03 [1] CRAN (R 3.6.1)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 3.6.1)
#  here                   0.1      2017-05-28 [1] CRAN (R 3.6.1)
#  htmltools              0.4.0    2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5.1    2019-10-08 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2    2019-09-11 [1] CRAN (R 3.6.1)
#  IRanges              * 2.20.2   2020-01-13 [1] Bioconductor
#  jsonlite               1.6.1    2020-02-02 [2] CRAN (R 3.6.1)
#  later                  1.0.0    2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38  2018-11-04 [3] CRAN (R 3.6.1)
#  lifecycle              0.2.0    2020-03-06 [1] CRAN (R 3.6.1)
#  magrittr               1.5      2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17   2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.56.0   2020-03-13 [1] CRAN (R 3.6.1)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 3.6.1)
#  pillar                 1.4.3    2019-12-20 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 3.6.1)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 3.6.1)
#  promises               1.1.0    2019-10-04 [1] CRAN (R 3.6.1)
#  pryr                   0.1.4    2018-02-18 [2] CRAN (R 3.6.1)
#  purrr                  0.3.4    2020-04-17 [1] CRAN (R 3.6.1)
#  R6                     2.4.1    2019-11-12 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.4    2020-03-17 [1] CRAN (R 3.6.1)
#  RCurl                  1.98-1.1 2020-01-19 [2] CRAN (R 3.6.1)
#  rlang                  0.4.5    2020-03-01 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4    2019-10-31 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              1.3-2    2018-01-03 [2] CRAN (R 3.6.1)
#  S4Vectors            * 0.24.3   2020-01-18 [1] Bioconductor
#  scales                 1.1.0    2019-11-18 [2] CRAN (R 3.6.1)
#  servr                  0.16     2020-03-02 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1    2018-11-05 [1] CRAN (R 3.6.1)
#  stringi                1.4.6    2020-02-17 [2] CRAN (R 3.6.1)
#  stringr                1.4.0    2019-02-10 [1] CRAN (R 3.6.1)
#  SummarizedExperiment * 1.16.1   2019-12-19 [1] Bioconductor
#  tibble                 3.0.1    2020-04-20 [1] CRAN (R 3.6.1)
#  tidyselect             1.0.0    2020-01-27 [2] CRAN (R 3.6.1)
#  vctrs                  0.2.4    2020-03-10 [1] CRAN (R 3.6.1)
#  withr                  2.2.0    2020-04-20 [1] CRAN (R 3.6.1)
#  xfun                   0.13     2020-04-13 [1] CRAN (R 3.6.1)
#  XVector                0.26.0   2019-10-29 [1] Bioconductor
#  zlibbioc               1.32.0   2019-10-29 [2] Bioconductor
#
# [1] /users/lcollado/R/3.6.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library
