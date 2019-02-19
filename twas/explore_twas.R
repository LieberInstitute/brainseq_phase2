library('tibble')
library('sessioninfo')
library('purrr')
library('SummarizedExperiment')
library('dplyr')

## Load twas data (from read_twas.R)
load('rda/twas.Rdata', verbose = TRUE)

## Read in the RSE info
## adapted from https://github.com/LieberInstitute/brainseq_phase2/blob/master/development/load_funs.R
load_rse <- function(type) {
    load_file <- file.path(
        '/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff',
        paste0('rse_', type, '.Rdata'))
    stopifnot(file.exists(load_file))
    load(load_file)

    ## Get the appropriate object
    if(type == 'gene') {
        rse <- rse_gene
    } else if (type == 'exon') {
        ## Drop those 4 exons not present in BrainSpan
        # rse <- rse_exon[-c(175584, 175585, 175586, 175604), ]
        rse <- rse_exon
    } else if (type == 'jxn') {
        rse <- rse_jxn
    } else if (type == 'tx') {
        rse <- rse_tx
    }
    ## Keep controls only
    # rse <- rse[, colData(rse)$Dx == 'Control']

    ## Set as factor
    colData(rse)$Region <- relevel(factor(colData(rse)$Region), 'DLPFC')
    colData(rse)$Race <- relevel(factor(colData(rse)$Race), ref = 'CAUC')
    colData(rse)$Sex <- relevel(factor(colData(rse)$Sex), ref = 'F')

    ## Add age linear splines
    fetal <- ifelse(colData(rse)$Age < 0, 1,0)
    birth <- colData(rse)$Age
    birth[birth < 0] <- 0 # linear spline
    infant <- colData(rse)$Age - 1
    infant[infant < 0] <- 0 # linear spline
    child <- colData(rse)$Age - 10
    child[child < 0] <- 0 # linear spline
    teen <- colData(rse)$Age - 20
    teen[teen < 0] <- 0 # linear spline
    adult <- colData(rse)$Age - 50
    adult[adult < 0] <- 0 # linear spline

    colData(rse)$fetal <- fetal
    colData(rse)$birth <- birth
    colData(rse)$infant <- infant
    colData(rse)$child <- child
    colData(rse)$teen <- teen
    colData(rse)$adult <- adult

    ## Add means
    colData(rse)$mean_mitoRate <- mean(colData(rse)$mitoRate)
    colData(rse)$mean_totalAssignedGene <- mean(colData(rse)$totalAssignedGene)
    ## Makes the design matrix not full rank in one of the models
    #    colData(rse)$mean_rRNA_rate <- mean(colData(rse)$rRNA_rate)
    colData(rse)$mean_RIN <- mean(colData(rse)$RIN)

    return(rse)
}

rse <- map(unique(twas$all$feature), load_rse)
names(rse) <- unique(twas$all$feature)

## Add the gene id
twas_exp <- map(twas, function(tw) {
    ## For testing the inner part of the function
    # tw <- twas$included
    
    by_feat <- split(tw, tw$feature)
    ## Make sure it's in the right order
    if(!identical(names(by_feat), names(rse))) {
        message(paste(Sys.time(), 'fixing order'))
        by_feat <- by_feat[names(rse)]
    }
    
    ## Now add the gene gencode ID and symbol
    result <- pmap_dfr(
        list(by_feat, rse, names(rse)),
        function(info, rs, feature) {
            
            ## Find the appropriate variables
            gene_var <- case_when(
                feature == 'gene' ~ 'gencodeID',
                feature == 'exon' ~ 'gencodeID',
                feature == 'jxn' ~ 'newGeneID',
                feature == 'tx' ~ 'gene_id'
            )
            symbol_var <- case_when(
                feature == 'gene' ~ 'Symbol',
                feature == 'exon' ~ 'Symbol',
                feature == 'jxn' ~ 'newGeneSymbol',
                feature == 'tx' ~ 'gene_name'
            )
            
            ## Match by id
            m <- match(info$ID, names(rowRanges(rs)))
            stopifnot(!is.na(m))
            
            ## Add the gene id/symbol
            info$geneid <- mcols(rowRanges(rs))[[gene_var]][m]
            info$genesymbol <- mcols(rowRanges(rs))[[symbol_var]][m]
            
            ## Done
            return(info)   
        }
    )
    
    return(result)
})
names(twas_exp) <- names(twas)

## Save the data for later use
dir.create('rda', showWarnings = FALSE)
save(twas_exp, file = 'rda/twas_exp.Rdata')




map(twas_exp, colnames)
# $all
#  [1] "PANEL"        "FILE"         "ID"           "CHR"          "P0"
#  [6] "P1"           "HSQ"          "BEST.GWAS.ID" "BEST.GWAS.Z"  "EQTL.ID"
# [11] "EQTL.R2"      "EQTL.Z"       "EQTL.GWAS.Z"  "NSNP"         "NWGT"
# [16] "MODEL"        "MODELCV.R2"   "MODELCV.PV"   "TWAS.Z"       "TWAS.P"
# [21] "chr"          "region"       "feature"      "type"         "geneid"
# [26] "genesymbol"
#
# $included
#  [1] "FILE"          "ID"            "TWAS.Z"        "TWAS.P"
#  [5] "JOINT.BETA"    "JOINT.BETA.SE" "JOINT.Z"       "JOINT.P"
#  [9] "chr"           "region"        "feature"       "type"
# [13] "geneid"        "genesymbol"
#
# $dropped
#  [1] "FILE"         "ID"           "TWAS.Z"       "TWAS.P"       "COND.BETA"
#  [6] "COND.BETA.SE" "COND.Z"       "COND.P"       "chr"          "region"
# [11] "feature"      "type"         "geneid"       "genesymbol"

walk(twas_exp, print, width = 120)
# # A tibble: 408,043 x 26
#    PANEL FILE  ID      CHR     P0     P1    HSQ BEST.GWAS.ID BEST.GWAS.Z EQTL.ID
#    <lgl> <chr> <chr> <dbl>  <dbl>  <dbl>  <dbl> <chr>              <dbl> <chr>
#  1 NA    /dcl… ENSG…     1 8.96e7 8.98e7 0.0498 rs10801784         -3.09 rs2065…
#  2 NA    /dcl… ENSG…     1 9.08e7 9.09e7 0.114  rs10922910         -3.58 rs1616…
#  3 NA    /dcl… ENSG…     1 9.13e7 9.14e7 0.453  rs13447450         -3.71 rs2819…
#  4 NA    /dcl… ENSG…     1 9.15e7 9.15e7 0.505  rs13447450         -3.71 rs1344…
#  5 NA    /dcl… ENSG…     1 9.22e7 9.22e7 0.168  rs4847377          -4.08 rs1737…
#  6 NA    /dcl… ENSG…     1 6.22e6 6.24e6 0.295  NA                 NA    NA
#  7 NA    /dcl… ENSG…     1 9.23e7 9.24e7 0.421  rs4847377          -4.08 rs3131…
#  8 NA    /dcl… ENSG…     1 9.25e7 9.28e7 0.0459 rs4847377          -4.08 rs6603…
#  9 NA    /dcl… ENSG…     1 9.28e7 9.30e7 0.157  rs4847377          -4.08 rs2391…
# 10 NA    /dcl… ENSG…     1 9.28e7 9.28e7 0.0318 rs4847377          -4.08 rs3424…
#     EQTL.R2 EQTL.Z EQTL.GWAS.Z  NSNP  NWGT MODEL MODELCV.R2 MODELCV.PV  TWAS.Z
#       <dbl>  <dbl>       <dbl> <dbl> <dbl> <chr>      <dbl>      <dbl>   <dbl>
#  1  0.0106    4.41      -1.96    380    23 enet      0.0179   4.10e- 3  -1.47
#  2  0.0171   -5.1        2.99    421    30 enet      0.0283   4.03e- 4  -2.56
#  3  0.313    11.8        0.363   545    58 enet      0.385    1.65e-44  -0.616
#  4  0.212    -9.72      -2.35    509    44 enet      0.407    9.19e-48   2.14
#  5  0.0358    6.13      -3.24    355    27 enet      0.123    2.31e-13  -1.93
#  6 NA        NA         NA         0     0 blump     0.167    6.16e-18  NA
#  7  0.0683   -5.77       2.28    313    63 enet      0.137    8.11e-15  -3.53
#  8  0.0294    4.35      -4.00    319    14 enet      0.0323   1.67e- 4  -3.17
#  9  0.0505    5.82      -2.23    267     5 lasso     0.0985   6.55e-11  -1.18
# 10  0.00988   4.17      -1.03    244     4 lasso     0.016    6.25e- 3  -1.59
# # … with 408,033 more rows, and 7 more variables: TWAS.P <dbl>, chr <chr>,
# #   region <chr>, feature <chr>, type <chr>, geneid <chr>, genesymbol <chr>
# # A tibble: 1,572 x 14
#    FILE      ID    TWAS.Z  TWAS.P JOINT.BETA JOINT.BETA.SE JOINT.Z JOINT.P chr
#    <chr>     <chr>  <dbl>   <dbl>      <dbl>         <dbl>   <dbl>   <dbl> <chr>
#  1 /dcl01/l… ENSG…   -4   5.50e-5       -4               1    -4   5.60e-5 1
#  2 /dcl01/l… ENSG…   -4.1 3.90e-5       -4.1             1    -4.1 3.80e-5 1
#  3 /dcl01/l… ENSG…   -5.4 8.30e-8       -5.4             1    -5.4 8.30e-8 1
#  4 /dcl01/l… ENSG…   -5.3 8.70e-8       -5.3             1    -5.3 8.80e-8 1
#  5 /dcl01/l… ENSG…   -5.5 3.20e-8       -5.5             1    -5.5 3.20e-8 1
#  6 /dcl01/l… ENSG…    4.5 5.40e-6        4.5             1     4.5 5.40e-6 1
#  7 /dcl01/l… ENSG…    4   6.20e-5        4               1     4   6.10e-5 1
#  8 /dcl01/l… ENSG…    4.2 2.70e-5        4.2             1     4.2 2.70e-5 1
#  9 /dcl01/l… ENSG…   -4.2 2.30e-5       -4.2             1    -4.2 2.20e-5 1
# 10 /dcl01/l… ENSG…    4.1 4.40e-5        4.1             1     4.1 4.30e-5 1
#    region feature type  geneid             genesymbol
#    <chr>  <chr>   <chr> <chr>              <chr>
#  1 DLPFC  gene    psycm ENSG00000007923.15 DNAJC11
#  2 DLPFC  gene    psycm ENSG00000155363.18 MOV10
#  3 DLPFC  gene    psycm ENSG00000143374.15 TARS2
#  4 DLPFC  gene    psycm ENSG00000143603.18 KCNN3
#  5 DLPFC  gene    psycm ENSG00000143537.13 ADAM15
#  6 DLPFC  gene    psycm ENSG00000117593.9  DARS2
#  7 DLPFC  gene    psycm ENSG00000023572.8  GLRX2
#  8 DLPFC  gene    psycm ENSG00000117335.19 CD46
#  9 DLPFC  gene    psycm ENSG00000185495.10 ""
# 10 DLPFC  gene    psycm ENSG00000054282.15 SDCCAG8
# # … with 1,562 more rows
# # A tibble: 254,181 x 14
#    FILE  ID    TWAS.Z  TWAS.P COND.BETA COND.BETA.SE COND.Z  COND.P chr   region
#    <chr> <chr>  <dbl>   <dbl>     <dbl>        <dbl>  <dbl>   <dbl> <chr> <chr>
#  1 /dcl… ENSG… -1.47  0.14       -1.47             1 -1.47  0.142   1     DLPFC
#  2 /dcl… ENSG… -2.56  0.01       -2.56             1 -2.56  0.0105  1     DLPFC
#  3 /dcl… ENSG… -0.616 0.54       -0.616            1 -0.616 0.538   1     DLPFC
#  4 /dcl… ENSG…  2.14  0.032       2.14             1  2.14  0.0323  1     DLPFC
#  5 /dcl… ENSG… -1.93  0.054      -1.93             1 -1.93  0.0536  1     DLPFC
#  6 /dcl… ENSG… -3.53  0.00042    -3.53             1 -3.53  0.00042 1     DLPFC
#  7 /dcl… ENSG… -3.17  0.0015     -3.17             1 -3.17  0.00152 1     DLPFC
#  8 /dcl… ENSG… -1.18  0.24       -1.18             1 -1.18  0.238   1     DLPFC
#  9 /dcl… ENSG… -1.59  0.11       -1.59             1 -1.59  0.112   1     DLPFC
# 10 /dcl… ENSG… -0.656 0.51       -0.656            1 -0.656 0.512   1     DLPFC
#    feature type  geneid             genesymbol
#    <chr>   <chr> <chr>              <chr>
#  1 gene    psycm ENSG00000171488.14 LRRC8C
#  2 gene    psycm ENSG00000233593.8  ""
#  3 gene    psycm ENSG00000162669.15 HFM1
#  4 gene    psycm ENSG00000097046.12 CDC7
#  5 gene    psycm ENSG00000273487.1  ""
#  6 gene    psycm ENSG00000122484.8  RPAP2
#  7 gene    psycm ENSG00000067208.14 EVI5
#  8 gene    psycm ENSG00000154511.11 FAM69A
#  9 gene    psycm ENSG00000207523.1  SNORA66
# 10 gene    psycm ENSG00000158292.6  GPR153
# # … with 254,171 more rows


# > summary(twas_exp$included$TWAS.P)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.000e+00 2.775e-07 5.300e-06 2.687e-05 2.600e-05 1.100e-03
# > summary(twas_exp$all$TWAS.P)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#    0.00    0.09    0.33    0.39    0.66    1.00  147856


## load eQTL raggr results
raggr_files <- c(
    'hippo_raggr' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/eqtl_tables/mergedEqtl_output_hippo_raggr_4features.rda',
    'dlpfc_raggr' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/eqtl_tables/mergedEqtl_output_dlpfc_raggr_4features.rda'
)
raggr <- map(raggr_files, function(x) {
    load(x, verbose = TRUE)
    return(allEqtl)
})
names(raggr) <- c('HIPPO', 'DLPFC')


map_dbl(raggr, ~ sum(.x$FDR < 0.01 ))
# HIPPO  DLPFC
# 66923 106438

# riskLoci <- read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/pgc_riskLoci.csv", stringsAsFactors=FALSE)


## Read in the files that Emily cleaned up at
## https://github.com/LieberInstitute/brainseq_phase2/blob/master/eQTL_GWAS_riskSNPs/create_eqtl_table_indexInfo.R
raggr_clean_files <- c(
    'HIPPO' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/raggr_179_snps_hippo_eqtls_fdr01.csv',
    'DLPFC' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/raggr_179_snps_dlp_eqtls_fdr01.csv'
)
raggr_clean <- map(raggr_clean_files, read.csv, stringsAsFactors = FALSE)
names(raggr_clean) <- names(raggr_clean_files)


## Start looking into how to match tables by SNP id..
hmm <- twas_exp$all$BEST.GWAS.ID
length(hmm)
# [1] 408043
length(unique(hmm))
# [1] 7976
hmm <- hmm[!is.na(hmm)]
length(hmm)
# [1] 260187
length(unique(hmm))
# [1] 7975

h <- unique(hmm)
table(grepl(':', h))
# FALSE  TRUE
 # 6037  1938
head(h[grepl(':', h)])

h_eqtl <- unique(twas_exp$all$EQTL.ID)
h_eqtl <- h_eqtl[!is.na(h_eqtl)]



table(raggr_clean$HIPPO$SNP %in% h)
# FALSE  TRUE
# 64971  1952
table(raggr_clean$DLPFC$SNP %in% h)
#  FALSE   TRUE
# 103381   3057

table(raggr_clean$HIPPO$SNP %in% h_eqtl)
# FALSE  TRUE
# 62369  4554
table(raggr_clean$DLPFC$SNP %in% h_eqtl)
# FALSE  TRUE
# 99802  6636


table(unique(raggr_clean$HIPPO$SNP) %in% h)
# FALSE  TRUE
#  5353   157
table(unique(raggr_clean$DLPFC$SNP) %in% h)
# FALSE  TRUE
#  6600   180

table(unique(raggr_clean$HIPPO$SNP) %in% h_eqtl)
# FALSE  TRUE
#  5302   208
table(unique(raggr_clean$DLPFC$SNP) %in% h_eqtl)
# FALSE  TRUE
#  6555   225


table(unique(raggr_clean$HIPPO$IndexSNP) %in% h)
# FALSE  TRUE
#   100     3
table(unique(raggr_clean$DLPFC$IndexSNP) %in% h)
# FALSE  TRUE
#   112     4

table(unique(raggr_clean$HIPPO$IndexSNP) %in% h_eqtl)
# FALSE
  # 103
table(unique(raggr_clean$DLPFC$IndexSNP) %in% h_eqtl)
# FALSE
#   116

## Read in the BIM files used for computing the weights (one per region)
## https://github.com/LieberInstitute/brainseq_phase2/blob/master/twas/filter_data/filter_snps.R#L128-L174
library('data.table')
bims <- map(c('DLPFC', 'HIPPO'), function(region) {
    bimfile <- paste0(
       '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered_',
       region,
       '.bim'
   )
   fread(bimfile,
       col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2'),
       colClasses = c('character', 'character', 'numeric', 'integer', 'character', 'character')
   )
})
names(bims) <- c('DLPFC', 'HIPPO')


## Check matches by SNP name
table(unique(raggr_clean$HIPPO$SNP) %in% bims$HIPPO$snp)
# FALSE  TRUE
#  4483  1027
table(unique(raggr_clean$DLPFC$SNP) %in% bims$DLPFC$snp)
# FALSE  TRUE
#  5586  1194


table(unique(raggr_clean$HIPPO$IndexSNP) %in% bims$HIPPO$snp)
# FALSE  TRUE
#   100     3
table(unique(raggr_clean$DLPFC$IndexSNP) %in% bims$DLPFC$snp)
# FALSE  TRUE
#   112     4

## Check with pairs of chr:position
## instead of SNP names
make_pairs <- function(rag, bim) {
    pairs_rag <- paste(
        gsub('chr', '', rag$chr_hg38),
        ':',
        rag$pos_hg38,
        sep = ''
    )
    pairs_bim <- paste(
        bim$chr,
        ':',
        bim$basepair,
        sep = ''
    )
    return(list(rag = pairs_rag, bim = pairs_bim))
}
pairs <- map2(raggr_clean, bims, make_pairs)

## Values of "trues" match the searches by SNP name...
with(pairs$HIPPO, table(unique(rag) %in% unique(bim)))
# FALSE  TRUE
#  4470  1027
with(pairs$DLPFC, table(unique(rag) %in% unique(bim)))
# FALSE  TRUE
#  5564  1194


## Read in the original BSP2 bim file with hg19 coordinates
bfile <- '/dcl01/lieber/ajaffe/Brain/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2.bim'
bsp2_bim <- fread(
    bfile,
    col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2')
)

## Check by SNP name...
table(unique(raggr_clean$HIPPO$SNP) %in% bsp2_bim$snp)
# TRUE
# 5510
table(unique(raggr_clean$DLPFC$SNP) %in% bsp2_bim$snp)
# TRUE
# 6780

## Ok, not all index snps were in our data to begin with
## which is why we did the raggr analysis in the
## first place
table(unique(raggr_clean$HIPPO$IndexSNP) %in% bsp2_bim$snp)
# FALSE  TRUE
#    30    73
table(unique(raggr_clean$DLPFC$IndexSNP) %in% bsp2_bim$snp)
# FALSE  TRUE
#    31    85


table(bims$HIPPO$snp %in% bsp2_bim$snp)
#    TRUE
# 1022527
table(bims$DLPFC$snp %in% bsp2_bim$snp)
#    TRUE
# 1022527



check_by_locus <- function(rag, ref) {
    by_loc <- split(rag$SNP, rag$IndexSNP)
    map_dbl(by_loc, ~ sum(.x %in% ref))
}

by_locus <- map(raggr_clean, check_by_locus, ref = h)
map(by_locus, ~ table(.x > 0))
# $HIPPO
#
# FALSE  TRUE
#    40    63
#
# $DLPFC
#
# FALSE  TRUE
#    48    68

by_locus_considered <- map2(
    raggr_clean,
    list(HIPPO = bims$HIPPO$snp, DLPFC = bims$DLPFC$snp),
    check_by_locus
)
map(by_locus_considered, ~ table(.x > 0))
# $HIPPO
#
# FALSE  TRUE
#    16    87
#
# $DLPFC
#
# FALSE  TRUE
#    17    99

map2(
    by_locus,
    by_locus_considered,
    ~ addmargins(table(
        'among TWAS-all best' = .x > 0,
        'among SNPs for TWAS weights' = .y > 0
    ))
)
# $HIPPO
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    16   24  40
#               TRUE      0   63  63
#               Sum      16   87 103
#
# $DLPFC
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    17   31  48
#               TRUE      0   68  68
#               Sum      17   99 116


by_locus_eqtl <- map(raggr_clean, check_by_locus, ref = h_eqtl)
map(by_locus_eqtl, ~ table(.x > 0))
# $HIPPO
#
# FALSE  TRUE
#    57    46
#
# $DLPFC
#
# FALSE  TRUE
#    67    49

map2(
    by_locus_eqtl,
    by_locus_considered,
    ~ addmargins(table(
        'among TWAS-all best' = .x > 0,
        'among SNPs for TWAS weights' = .y > 0
    ))
)
# $HIPPO
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    16   41  57
#               TRUE      0   46  46
#               Sum      16   87 103
#
# $DLPFC
#                    among SNPs for TWAS weights
# among TWAS-all best FALSE TRUE Sum
#               FALSE    17   50  67
#               TRUE      0   49  49
#               Sum      17   99 116



gene_by_locus <- function(rag, ref) {
    by_loc <- split(rag$gene, rag$IndexSNP)
    map_dbl(by_loc, ~ sum(.x %in% ref))
}


g_by_locus <- map(names(rse), function(feature) {
    map(
        map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
        gene_by_locus,
        ref = twas_exp$all$ID[twas_exp$all$feature == feature]
    )
})
names(g_by_locus) <- names(rse)
map(g_by_locus, function(x) {
    r <- map_dfr(x, ~ table(.x > 0))
    r$state <- c(FALSE, TRUE)
    return(r)
})
# $gene
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1     5    11 FALSE
# 2    45    56 TRUE
#
# $exon
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1     9    12 FALSE
# 2    68    75 TRUE
#
# $jxn
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    12    19 FALSE
# 2    76    79 TRUE
#
# $tx
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1     9     9 FALSE
# 2    56    67 TRUE


### Ehem.... matching gene ids vs snp ids... ehem... T_T
# g_by_locus_considered <- map(names(rse), function(feature) {
#     map2(
#         map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
#         list(HIPPO = bims$HIPPO$snp, DLPFC = bims$DLPFC$snp),
#         gene_by_locus
#     )
# })
# names(g_by_locus_considered) <- names(rse)
# map(g_by_locus_considered, function(x) {
#     r <- map_dfr(x, ~ table(factor(.x > 0, levels = c('FALSE', 'TRUE'))))
#     r$state <- c(FALSE, TRUE)
#     return(r)
# })
#

clean_tabs <- function(l) {
    map2_dfr(l, names(l), function(x, y) {
        x$feature <- y
        return(x)
    })
}

g_by_locus_reg <- map(names(rse), function(feature) {
    map2(
        map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
        map(names(raggr_clean), ~ twas_exp$all$ID[twas_exp$all$feature == feature & twas_exp$all$region == .x]),
        gene_by_locus
    )
})
names(g_by_locus_reg) <- names(rse)
clean_tabs(map(g_by_locus_reg, function(x) {
    r <- map_dfr(x, ~ table(.x > 0))
    r$state <- c(FALSE, TRUE)
    return(r)
}))
# # A tibble: 8 x 4
#   HIPPO DLPFC state feature
#   <int> <int> <lgl> <chr>
# 1     6    11 FALSE gene
# 2    44    56 TRUE  gene
# 3    11    13 FALSE exon
# 4    66    74 TRUE  exon
# 5    13    19 FALSE jxn
# 6    75    79 TRUE  jxn
# 7    12     9 FALSE tx
# 8    53    67 TRUE  tx

clean_tabs_type <- function(l) {
    map2_dfr(l, names(l), function(x, y) {
        x$type <- y
        return(x)
    })
}


g_by_type <- map(c('psycm', 'pgc2'), function(type) {
    g_reg <- map(names(rse), function(feature) {
        map2(
            map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
            map(names(raggr_clean), ~ twas_exp$all$ID[twas_exp$all$feature == feature & twas_exp$all$region == .x & twas_exp$all$type == type]),
            gene_by_locus
        )
    })
    names(g_reg) <- names(rse)
    map(g_reg, function(x) {
        r <- map_dfr(x, ~ table(.x > 0))
        r$state <- c(FALSE, TRUE)
        return(r)
    })
})
names(g_by_type) <- c('psycm', 'pgc2')
clean_tabs_type(map(g_by_type, clean_tabs))
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1     6    11 FALSE gene    psycm
#  2    44    56 TRUE  gene    psycm
#  3    17    13 FALSE exon    psycm
#  4    60    74 TRUE  exon    psycm
#  5    13    19 FALSE jxn     psycm
#  6    75    79 TRUE  jxn     psycm
#  7    12     9 FALSE tx      psycm
#  8    53    67 TRUE  tx      psycm
#  9     6    11 FALSE gene    pgc2
# 10    44    56 TRUE  gene    pgc2
# 11    11    13 FALSE exon    pgc2
# 12    66    74 TRUE  exon    pgc2
# 13    13    19 FALSE jxn     pgc2
# 14    75    79 TRUE  jxn     pgc2
# 15    12     9 FALSE tx      pgc2
# 16    53    67 TRUE  tx      pgc2


g_by_type_inc <- map(c('psycm', 'pgc2'), function(type) {
    g_reg <- map(names(rse), function(feature) {
        map2(
            map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
            map(names(raggr_clean), ~ twas_exp$included$ID[twas_exp$included$feature == feature & twas_exp$included$region == .x & twas_exp$included$type == type]),
            gene_by_locus
        )
    })
    names(g_reg) <- names(rse)
    map(g_reg, function(x) {
        r <- map_dfr(x, ~ table(.x > 0))
        r$state <- c(FALSE, TRUE)
        return(r)
    })
})
names(g_by_type_inc) <- c('psycm', 'pgc2')
clean_tabs_type(map(g_by_type_inc, clean_tabs))
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1    30    38 FALSE gene    psycm
#  2    20    29 TRUE  gene    psycm
#  3    46    45 FALSE exon    psycm
#  4    31    42 TRUE  exon    psycm
#  5    48    48 FALSE jxn     psycm
#  6    40    50 TRUE  jxn     psycm
#  7    39    41 FALSE tx      psycm
#  8    26    35 TRUE  tx      psycm
#  9    32    42 FALSE gene    pgc2
# 10    18    25 TRUE  gene    pgc2
# 11    51    58 FALSE exon    pgc2
# 12    26    29 TRUE  exon    pgc2
# 13    61    57 FALSE jxn     pgc2
# 14    27    41 TRUE  jxn     pgc2
# 15    43    44 FALSE tx      pgc2
# 16    22    32 TRUE  tx      pgc2




check_by_feature <- function(tw, cut) {
    res <- map(c('psycm', 'pgc2'), function(type) {
        g_reg <- map(names(rse), function(feature) {
            map2(
                map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
                map(names(raggr_clean), ~ tw$ID[tw$feature == feature & tw$region == .x & tw$type == type & tw$TWAS.P < cut]),
                ~ match(unique(.x$gene), .y)
            )
        })
        names(g_reg) <- names(rse)
        map(g_reg, function(x) {
            # r <- map_dfr(x, ~ table(factor(!is.na(.x), levels = c('FALSE', 'TRUE'))))
            r <- map_dfr(x, ~ table(!is.na(.x)))
            r$state <- c(FALSE, TRUE)
            return(r)
        })
    })
    names(res) <- c('psycm', 'pgc2')
    
    clean_tabs_type(map(res, clean_tabs))
}

check_by_feature(twas_exp$all, cut = 1.01)
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1    53    66 FALSE gene    psycm
#  2    70   105 TRUE  gene    psycm
#  3   408   669 FALSE exon    psycm
#  4   449   694 TRUE  exon    psycm
#  5   232   295 FALSE jxn     psycm
#  6   275   364 TRUE  jxn     psycm
#  7   104   161 FALSE tx      psycm
#  8   140   171 TRUE  tx      psycm
#  9    53    66 FALSE gene    pgc2
# 10    70   105 TRUE  gene    pgc2
# 11   395   669 FALSE exon    pgc2
# 12   462   694 TRUE  exon    pgc2
# 13   232   295 FALSE jxn     pgc2
# 14   275   364 TRUE  jxn     pgc2
# 15   104   161 FALSE tx      pgc2
# 16   140   171 TRUE  tx      pgc2

check_by_feature(twas_exp$all, cut = 0.01)
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1    75    98 FALSE gene    psycm
#  2    48    73 TRUE  gene    psycm
#  3   511   841 FALSE exon    psycm
#  4   346   522 TRUE  exon    psycm
#  5   307   386 FALSE jxn     psycm
#  6   200   273 TRUE  jxn     psycm
#  7   140   210 FALSE tx      psycm
#  8   104   122 TRUE  tx      psycm
#  9    83   106 FALSE gene    pgc2
# 10    40    65 TRUE  gene    pgc2
# 11   520   907 FALSE exon    pgc2
# 12   337   456 TRUE  exon    pgc2
# 13   331   419 FALSE jxn     pgc2
# 14   176   240 TRUE  jxn     pgc2
# 15   148   225 FALSE tx      pgc2
# 16    96   107 TRUE  tx      pgc2

check_by_feature(twas_exp$included, cut = 0.01)
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1   105   141 FALSE gene    psycm
#  2    18    30 TRUE  gene    psycm
#  3   827  1321 FALSE exon    psycm
#  4    30    42 TRUE  exon    psycm
#  5   467   610 FALSE jxn     psycm
#  6    40    49 TRUE  jxn     psycm
#  7   220   300 FALSE tx      psycm
#  8    24    32 TRUE  tx      psycm
#  9   107   145 FALSE gene    pgc2
# 10    16    26 TRUE  gene    pgc2
# 11   833  1335 FALSE exon    pgc2
# 12    24    28 TRUE  exon    pgc2
# 13   480   619 FALSE jxn     pgc2
# 14    27    40 TRUE  jxn     pgc2
# 15   224   302 FALSE tx      pgc2
# 16    20    30 TRUE  tx      pgc2



check_by_feature_rev <- function(tw, cut) {
    res <- map(c('psycm', 'pgc2'), function(type) {
        g_reg <- map(names(rse), function(feature) {
            map2(
                map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
                map(names(raggr_clean), ~ tw$ID[tw$feature == feature & tw$region == .x & tw$type == type & tw$TWAS.P < cut]),
                ~ match(unique(.y[!is.na(.y)]), .x$gene)
            )
        })
        names(g_reg) <- names(rse)
        map(g_reg, function(x) {
            # r <- map_dfr(x, ~ table(factor(!is.na(.x), levels = c('FALSE', 'TRUE'))))
            r <- map_dfr(x, ~ table(!is.na(.x)))
            r$state <- c(FALSE, TRUE)
            return(r)
        })
    })
    names(res) <- c('psycm', 'pgc2')
    
    clean_tabs_type(map(res, clean_tabs))
}

check_by_feature_rev(twas_exp$all, cut = 1.01)
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1  3907  5377 FALSE gene    psycm
#  2    70   105 TRUE  gene    psycm
#  3 25237 38437 FALSE exon    psycm
#  4   449   694 TRUE  exon    psycm
#  5 15832 20289 FALSE jxn     psycm
#  6   275   364 TRUE  jxn     psycm
#  7  7014  8890 FALSE tx      psycm
#  8   140   171 TRUE  tx      psycm
#  9  3960  5450 FALSE gene    pgc2
# 10    70   105 TRUE  gene    pgc2
# 11 29404 39031 FALSE exon    pgc2
# 12   462   694 TRUE  exon    pgc2
# 13 16087 20556 FALSE jxn     pgc2
# 14   275   364 TRUE  jxn     pgc2
# 15  7132  9035 FALSE tx      pgc2
# 16   140   171 TRUE  tx      pgc2

check_by_feature_rev(twas_exp$all, cut = 0.01)
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1   354   514 FALSE gene    psycm
#  2    48    73 TRUE  gene    psycm
#  3  2389  3689 FALSE exon    psycm
#  4   346   522 TRUE  exon    psycm
#  5  1533  1913 FALSE jxn     psycm
#  6   200   273 TRUE  jxn     psycm
#  7   720   897 FALSE tx      psycm
#  8   104   122 TRUE  tx      psycm
#  9   264   391 FALSE gene    pgc2
# 10    40    65 TRUE  gene    pgc2
# 11  2142  2904 FALSE exon    pgc2
# 12   337   456 TRUE  exon    pgc2
# 13  1235  1519 FALSE jxn     pgc2
# 14   176   240 TRUE  jxn     pgc2
# 15   608   691 FALSE tx      pgc2
# 16    96   107 TRUE  tx      pgc2

check_by_feature_rev(twas_exp$included, cut = 0.01)
# # A tibble: 16 x 5
#    HIPPO DLPFC state feature type
#    <int> <int> <lgl> <chr>   <chr>
#  1    60    96 FALSE gene    psycm
#  2    18    30 TRUE  gene    psycm
#  3    62    74 FALSE exon    psycm
#  4    30    42 TRUE  exon    psycm
#  5    95    93 FALSE jxn     psycm
#  6    40    49 TRUE  jxn     psycm
#  7    89    98 FALSE tx      psycm
#  8    24    32 TRUE  tx      psycm
#  9    40    55 FALSE gene    pgc2
# 10    16    26 TRUE  gene    pgc2
# 11    45    57 FALSE exon    pgc2
# 12    24    28 TRUE  exon    pgc2
# 13    67    53 FALSE jxn     pgc2
# 14    27    40 TRUE  jxn     pgc2
# 15    58    54 FALSE tx      pgc2
# 16    20    30 TRUE  tx      pgc2








## put this on pause...


library('UpSetR')

tw <- twas_exp$included

library('VennDiagram')
library('RColorBrewer')
venn_cols <- brewer.pal('Set1', n = 4)

genes <- with(subset(tw, type == 'psycm'), split(geneid, feature))
genes <- genes[names(rse)]
genes <- map(genes, ~ .x[!is.na(.x)])

names(venn_cols) <- names(genes)

make_venn <- function(genes, title = 'DE features grouped by gene id') {
    v <- venn.diagram(genes, filename = NULL,
        main = title,
        col = 'transparent', fill = venn_cols[names(genes)],
        alpha = 0.5, margin = 0,
        main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
        cat.col = venn_cols[names(genes)])
    grid.newpage()
    grid.draw(v)
}

dir.create('pdf', showWarnings = FALSE)

pdf('pdf/test.pdf', useDingbats = FALSE)
make_venn(genes)
dev.off()
system('rm VennDiagram*')

library('gplots')
genes2 <- map(genes, unique)
names(genes2) <- names(genes)

x <- venn(genes2, show.plot = FALSE)

## Get the matrix out
y <- matrix(x, ncol = ncol(x), dimnames = attr(x, 'dimnames'))
z <- map_dbl(names(genes), ~ sum(y[ y[, .x] > 0, 'num']))
names(z) <- names(genes)
z



tw_up <- map_dfr(split(tw, tw$geneid), function(i) {
    
    data.frame(
        geneid = unique(i$geneid),
        genesymbol = unique(i$genesymbol),
        gene = ifelse('gene' %in% i$feature, 1, 0),
        exon = ifelse('exon' %in% i$feature, 1, 0),
        jxn = ifelse('jxn' %in% i$feature, 1, 0),
        tx = ifelse('tx' %in% i$feature, 1, 0),
        DLPFC = ifelse('DLPFC' %in% i$region, 1, 0),
        HIPPO = ifelse('HIPPO' %in% i$region, 1, 0),
        psycm = ifelse('psycm' %in% i$type, 1, 0),
        pgc2 = ifelse('pgc2' %in% i$type, 1, 0),
        stringsAsFactors = FALSE
    )
    
})

set_meta <- data.frame(
    sets = c(names(rse), 'DLPFC', 'HIPPO', 'psycm', 'pgc2'),
    colors = c(brewer.pal('Set1', n = 4), 'dark orange', 'skyblue3', 'sienna1', 'springgreen4'),
    stringsAsFactors = FALSE
)
cols <- set_meta$colors
names(cols) <- set_meta$colors


pdf('pdf/test.pdf', useDingbats = FALSE, width = 14)
upset(tw_up,
    nset = 8,
    sets = c(names(rse), 'DLPFC', 'HIPPO', 'psycm', 'pgc2'),
    # empty.intersections = 'on',
    group.by = 'sets',
    # nintersects = NA,
    nintersects = 70,
    keep.order = TRUE,
    set.metadata = list(
        data = set_meta,
        plots = list(
            list(
                type = 'matrix_rows',
                column = 'colors',
                colors = cols
            )
        )
    )
)
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()