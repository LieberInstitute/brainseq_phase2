## A cleaner and more focused script that explore_twas.R

library('tibble')
library('sessioninfo')
library('purrr')
# library('readr')
library('ggplot2')
library('gplots')
library('VennDiagram')
library('RColorBrewer')
library('readxl')
source('twas_functions.R')

load('rda/twas_exp.Rdata', verbose = TRUE)

## Andrew's exploration code that focuses on the 'all' part
tt <- twas_exp$all
## Drop TWAS NA p-values
tt <- tt[!is.na(tt$TWAS.P), ]
## Focus on PGC2 GWAS
tt <- tt[which(tt$type == "pgc2"),]

## Load big snpMap table
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda", verbose = TRUE)
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)
snpMap$pos_hg38_info <- paste0(gsub('^chr', '', snpMap$chr_hg38), ":", snpMap$pos_hg38)
## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap = snpMap[-snpInd,]

## Match to the big snpMap table
m_to_fullMap <- match(tt$BEST.GWAS.ID, snpMap$SNP)
stopifnot(!any(is.na(m_to_fullMap)))
m_to_fullMap_qtl <- match(tt$EQTL.ID, snpMap$SNP)
stopifnot(!any(is.na(m_to_fullMap_qtl)))

## Add pos hg19 and pos hg38
tt$BEST.GWAS.pos_hg19 <- snpMap$pos_hg19[m_to_fullMap]
tt$BEST.GWAS.pos_hg38 <- snpMap$pos_hg38_info[m_to_fullMap]

tt$EQTL.pos_hg19 <- snpMap$pos_hg19[m_to_fullMap_qtl]
tt$EQTL.pos_hg38 <- snpMap$pos_hg38_info[m_to_fullMap_qtl]

## Add computed P-values for BEST GWAS and EQTL
tt$BEST.GWAS.P.computed <- 2*(pnorm( abs(tt$BEST.GWAS.Z ) , lower.tail=F ))
tt$EQTL.P.computed <- 2*(pnorm( abs(tt$EQTL.GWAS.Z ) , lower.tail=F ))

## Compute FDR/Bonf by region for each feature 4 features
tt <- map_dfr(split(tt, tt$region), function(reg) {
    res <- map_dfr(split(reg, reg$feature), function(reg_feat) {
        reg_feat$TWAS.FDR <- p.adjust(reg_feat$TWAS.P, 'fdr')
        reg_feat$TWAS.Bonf <- p.adjust(reg_feat$TWAS.P, 'bonf')
        reg_feat$BEST.GWAS.FDR.computed <- p.adjust(reg_feat$BEST.GWAS.P.computed, 'fdr')
        reg_feat$EQTL.FDR.computed <- p.adjust(reg_feat$EQTL.P.computed, 'fdr')
        return(reg_feat)
    })
    return(res[order(res$TWAS.P), ])
})

## Add raggr data
## Code based on https://github.com/LieberInstitute/brainseq_phase2/blob/master/eQTL_GWAS_riskSNPs/create_eqtl_table_indexInfo.R


## risk loci from PGC paper
indexLoci <- read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/pgc_riskLoci.csv", stringsAsFactors=FALSE)
indexLoci$hg19POS = paste0(indexLoci$Chromosome, ":", indexLoci$snp_pos_hg19)

## risk loci from PGC paper + rAggr proxy markers
riskLoci <- read.csv("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/rAggr_results_179.csv", stringsAsFactors=FALSE)
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
length(unique(riskLoci$SNP2_Name))
# [1] 10981
riskLoci$hg19POS1 = paste0(riskLoci$SNP1_Chr, ":", riskLoci$SNP1_Pos) 
riskLoci$hg19POS2 = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos)
length(unique(riskLoci$hg19POS2))
# [1] 10975

snpMap$Status <- 'Other'
snpMap$Status[snpMap$pos_hg19 %in% riskLoci$hg19POS2] <- 'Proxy'
snpMap$Status[snpMap$pos_hg19 %in% indexLoci$hg19POS] <- 'Index'

get_proxy_info <- function(pos_hg19, prefix) {
    status <- rep('Other', length(pos_hg19))
    status[pos_hg19 %in% riskLoci$hg19POS2] <- 'Proxy'
    status[pos_hg19 %in% indexLoci$hg19POS] <- 'Index'
    print(table(status))
    
    indexSNP <- rep(NA, length(pos_hg19))
    m <- match(pos_hg19, riskLoci$hg19POS2)
    stopifnot(sum(is.na(m)) == sum(status == 'Other'))
    indexSNP[!is.na(m)] <- riskLoci$SNP1_Name[m[!is.na(m)]]
    
    indexSNP_pos_hg19 <- rep(NA, length(pos_hg19))
    indexSNP_pos_hg19[!is.na(m)] <- riskLoci$hg19POS1[m[!is.na(m)]]
    
    distance <- rep(NA, length(pos_hg19))
    distance[!is.na(m)] <- riskLoci$Distance[m[!is.na(m)]]
    
    
    
    res <- tibble(
        status = status,
        indexSNP = indexSNP,
        indexSNP_pos_hg19 = indexSNP_pos_hg19,
        indexSNP_distance = distance
    )
    
    print(length(unique(res$indexSNP[!is.na(res$indexSNP)])))
    print(summary(res$indexSNP_distance[!is.na(res$indexSNP_distance)]))
    
    colnames(res) <- paste0(prefix, colnames(res))
    return(res)
}
tt <- as_tibble(cbind(tt,
    get_proxy_info(tt$BEST.GWAS.pos_hg19, 'BEST.GWAS.'),
    get_proxy_info(tt$EQTL.pos_hg19, 'EQTL.')
))
## BEST.GWAS info
# status
#  Index  Other  Proxy
#   1317 127039   4580
# [1] 93
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -331128   -1922       0    4371   21445  322056
## EQTL info
# status
#  Index  Other  Proxy
#      3 132415    518
# [1] 51
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -339607  -72415  -14170  -16679   34378  437328

print(tt, width = 600)
head(as.data.frame(tt))


table(table(tt$ID))
#     1     2
# 69672 31632

ids <- unique(tt$ID)
is_DLPFC <- tt$region == 'DLPFC'

i_DLPFC <- which(is_DLPFC)[match(ids, tt$ID[is_DLPFC])]
i_HIPPO <- which(!is_DLPFC)[match(ids, tt$ID[!is_DLPFC])]
stopifnot(
    all(table(tt$region) - c(
        sum(tt$region[i_DLPFC] == 'DLPFC', na.rm = TRUE),
        sum(tt$region[i_HIPPO] == 'HIPPO', na.rm = TRUE)
    ) == 0)
)

ttReg_map <- data.frame(
    ID = ids,
    i_DLPFC = i_DLPFC,
    i_HIPPO = i_HIPPO,
    stringsAsFactors = FALSE
)


get_variable_by_region <- function(var, NAs_0 = FALSE) {
    
    m <- match(ttReg_map$ID, tt$ID)
    
    res <- data.frame(
        ID = ttReg_map$ID,
        feature = tt$feature[m],
        geneid = tt$geneid[m],
        genesymbol = tt$genesymbol[m],
        DLPFC = tt[ttReg_map$i_DLPFC, var, drop = TRUE],
        HIPPO = tt[ttReg_map$i_HIPPO, var, drop = TRUE],
        in_both = !is.na(ttReg_map$i_DLPFC) & !is.na(ttReg_map$i_HIPPO),
        TWAS.FDR_DLPFC = tt$TWAS.FDR[ttReg_map$i_DLPFC],
        TWAS.FDR_HIPPO = tt$TWAS.FDR[ttReg_map$i_HIPPO],
        TWAS.Bonf_DLPFC = tt$TWAS.Bonf[ttReg_map$i_DLPFC],
        TWAS.Bonf_HIPPO = tt$TWAS.Bonf[ttReg_map$i_HIPPO],
        BEST.GWAS.status_DLPFC = tt$BEST.GWAS.status[ttReg_map$i_DLPFC],
        BEST.GWAS.status_HIPPO = tt$BEST.GWAS.status[ttReg_map$i_HIPPO],
        stringsAsFactors = FALSE
    )

    res$FDR.5perc <- 'None'
    res$FDR.5perc[res$TWAS.FDR_DLPFC < 0.05] <- 'DLPFC'
    res$FDR.5perc[res$TWAS.FDR_HIPPO < 0.05] <- 'HIPPO'
    res$FDR.5perc[res$TWAS.FDR_DLPFC < 0.05 & res$TWAS.FDR_HIPPO < 0.05] <- 'Both'
    
    res$Bonf.5perc <- 'None'
    res$Bonf.5perc[res$TWAS.Bonf_DLPFC < 0.05] <- 'DLPFC'
    res$Bonf.5perc[res$TWAS.Bonf_HIPPO < 0.05] <- 'HIPPO'
    res$Bonf.5perc[res$TWAS.Bonf_DLPFC < 0.05 & res$TWAS.Bonf_HIPPO < 0.05] <- 'Both'

    res$BEST.GWAS.status <- 'Other'
    res$BEST.GWAS.status[c(which(res$BEST.GWAS.status_DLPFC != 'Other'),  which(res$BEST.GWAS.status_HIPPO != 'Other'))] <- 'Risk Locus'
    
    if(NAs_0 == TRUE) {
        res$DLPFC[is.na(res$DLPFC)] <- 0
        res$HIPPO[is.na(res$HIPPO)] <- 0
    }
    
    ## Make the features as factor, so its looks ok when plotting
    res$feature <- factor(res$feature, levels = c('gene', 'exon', 'jxn', 'tx'))
    res$FDR.5perc <- factor(res$FDR.5perc, levels = c('None', 'DLPFC', 'HIPPO', 'Both'))
    res$Bonf.5perc <- factor(res$Bonf.5perc, levels = c('None', 'DLPFC', 'HIPPO', 'Both'))
    res$BEST.GWAS.status <- factor(res$BEST.GWAS.status, levels = c('Other', 'Risk Locus'))
    
    return(res)
}


region_twas_z <- get_variable_by_region('TWAS.Z', NAs_0 = TRUE)

table(region_twas_z$in_both)
# FALSE  TRUE
# 69672 31632
table(region_twas_z$in_both) / nrow(region_twas_z) * 100
#    FALSE     TRUE
# 68.77517 31.22483
addmargins(table(
    'in both' = region_twas_z$in_both,
    'FDR < 0.05' = region_twas_z$FDR.5perc,
    useNA = 'ifany'
))
#        FDR < 0.05
# in both   None  DLPFC  HIPPO   Both    Sum
#   FALSE  66336   2112   1224      0  69672
#   TRUE   29417    487    481   1247  31632
#   Sum    95753   2599   1705   1247 101304

addmargins(table(
    'in both' = region_twas_z$in_both,
    'Bonf < 0.05' = region_twas_z$Bonf.5perc,
    useNA = 'ifany'
))
#        Bonf < 0.05
# in both   None  DLPFC  HIPPO   Both    Sum
#   FALSE  69279    240    153      0  69672
#   TRUE   31353     87     71    121  31632
#   Sum   100632    327    224    121 101304

table(region_twas_z$BEST.GWAS.status)
# Other Risk Locus
# 96753       4551



pdf('pdf/pgc2_twas_z.pdf', useDingbats = FALSE, width = 24, height = 14)
ggplot(region_twas_z,
    aes(x = DLPFC, y = HIPPO, color = FDR.5perc, shape = in_both)) +
    geom_point() +
    facet_grid(BEST.GWAS.status ~ feature) +
    coord_fixed() +
    theme_bw(base_size = 30) +
    ggtitle('TWAS Z by brain region') +
    scale_color_manual(values = c('grey80', 'dark orange', 'skyblue3', 'purple'))
    
ggplot(region_twas_z,
    aes(x = DLPFC, y = HIPPO, color = Bonf.5perc, shape = in_both)) +
    geom_point() +
    facet_grid(BEST.GWAS.status ~ feature) +
    coord_fixed() +
    theme_bw(base_size = 30) +
    ggtitle('TWAS Z by brain region') +
    scale_color_manual(values = c('grey80', 'dark orange', 'skyblue3', 'purple'))
dev.off()


## Get the numbers of points in the different parts of the plot
map(split(region_twas_z, region_twas_z$feature),
    ~ map(split(.x, .x$BEST.GWAS.status), 
        ~ addmargins(table('FDR <5%' = .x$FDR.5perc, 'In both' = .x$in_both))
    )
)
# $gene
# $gene$Other
#        In both
# FDR <5% FALSE TRUE  Sum
#   None   3633 2557 6190
#   DLPFC    95   54  149
#   HIPPO    31   33   64
#   Both      0   70   70
#   Sum    3759 2714 6473
#
# $gene$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE Sum
#   None    144   78 222
#   DLPFC    28    7  35
#   HIPPO     8    4  12
#   Both      0   20  20
#   Sum     180  109 289
#
#
# $exon
# $exon$Other
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None  34312 14515 48827
#   DLPFC   846   221  1067
#   HIPPO   488   216   704
#   Both      0   431   431
#   Sum   35646 15383 51029
#
# $exon$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE  Sum
#   None   1412  465 1877
#   DLPFC   291   35  326
#   HIPPO   122   30  152
#   Both      0  147  147
#   Sum    1825  677 2502
#
#
# $jxn
# $jxn$Other
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None  18572  7680 26252
#   DLPFC   463    98   561
#   HIPPO   314    94   408
#   Both      0   308   308
#   Sum   19349  8180 27529
#
# $jxn$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE  Sum
#   None    607  266  873
#   DLPFC   133   13  146
#   HIPPO    73   22   95
#   Both      0   79   79
#   Sum     813  380 1193
#
#
# $tx
# $tx$Other
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None   7379  3727 11106
#   DLPFC   188    53   241
#   HIPPO   146    76   222
#   Both      0   153   153
#   Sum    7713  4009 11722
#
# $tx$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE Sum
#   None    277  129 406
#   DLPFC    68    6  74
#   HIPPO    42    6  48
#   Both      0   39  39
#   Sum     387  180 567


## Now for Bonf
map(split(region_twas_z, region_twas_z$feature),
    ~ map(split(.x, .x$BEST.GWAS.status), 
        ~ addmargins(table('Bonf <5%' = .x$Bonf.5perc, 'In both' = .x$in_both))
    )
)
# $gene
# $gene$Other
#         In both
# Bonf <5% FALSE TRUE  Sum
#    None   3740 2699 6439
#    DLPFC    14    7   21
#    HIPPO     5    3    8
#    Both      0    5    5
#    Sum    3759 2714 6473
#
# $gene$`Risk Locus`
#         In both
# Bonf <5% FALSE TRUE Sum
#    None    166   92 258
#    DLPFC    10    3  13
#    HIPPO     4    5   9
#    Both      0    9   9
#    Sum     180  109 289
#
#
# $exon
# $exon$Other
#         In both
# Bonf <5% FALSE  TRUE   Sum
#    None  35552 15348 50900
#    DLPFC    55     9    64
#    HIPPO    39    16    55
#    Both      0    10    10
#    Sum   35646 15383 51029
#
# $exon$`Risk Locus`
#         In both
# Bonf <5% FALSE TRUE  Sum
#    None   1721  602 2323
#    DLPFC    68   24   92
#    HIPPO    36   13   49
#    Both      0   38   38
#    Sum    1825  677 2502
#
#
# $jxn
# $jxn$Other
#         In both
# Bonf <5% FALSE  TRUE   Sum
#    None  19297  8139 27436
#    DLPFC    32    14    46
#    HIPPO    20    14    34
#    Both      0    13    13
#    Sum   19349  8180 27529
#
# $jxn$`Risk Locus`
#         In both
# Bonf <5% FALSE TRUE  Sum
#    None    765  342 1107
#    DLPFC    29   11   40
#    HIPPO    19    8   27
#    Both      0   19   19
#    Sum     813  380 1193
#
#
# $tx
# $tx$Other
#         In both
# Bonf <5% FALSE  TRUE   Sum
#    None   7689  3978 11667
#    DLPFC    11    13    24
#    HIPPO    13     5    18
#    Both      0    13    13
#    Sum    7713  4009 11722
#
# $tx$`Risk Locus`
#         In both
# Bonf <5% FALSE TRUE Sum
#    None    349  153 502
#    DLPFC    21    6  27
#    HIPPO    17    7  24
#    Both      0   14  14
#    Sum     387  180 567


## Load SCZD case-control data
## From https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/explore_case_control.R#L65-L71
outFeat <- lapply(c('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda', '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda'), function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    outTx$ensemblID <- gsub('\\..*', '', outTx$gene_id)
    return(list('gene' = outGene, 'exon' = outExon, 'jxn' = outJxn, 'tx' = outTx))
})
names(outFeat) <- c('DLPFC', 'HIPPO')


## Add the SCZD info back to the tt object
tt$SCZD_FDR <- tt$SCZD_pvalue <- tt$SCZD_t <- NA
for(feat in features) {
    for(region in names(outFeat)) {
        message(paste(Sys.time(), 'processing', region, 'at the', feat, 'level'))
        i <- which(tt$feature == feat & tt$region == region)
        current <- outFeat[[region]][[feat]]
        print(length(i))
        # table(tt$region[i], tt$feature[i])
        m <- match(tt$ID[i], rownames(current))
        stopifnot(!any(is.na(m)))
        tt$SCZD_t[i] <- current$t[m]
        tt$SCZD_pvalue[i] <- current$P.Value[m]
        tt$SCZD_FDR[i] <- current$adj.P.Val[m]
    }
}
## For checking the code
# sapply(tt[, 45:47], function(x) { sum(is.na(x)) })


## Subset by significant (TWAS FDR < 5%)
ttSig <- map(split(tt, tt$region), ~ .x[.x$TWAS.FDR < 0.05, ])
map_dfr(ttSig, dim)
# # A tibble: 2 x 2
#   DLPFC HIPPO
#   <int> <int>
# 1  3846  2952
# 2    47    47

ttSig_bonf <- map(split(tt, tt$region), ~ .x[.x$TWAS.Bonf < 0.05, ])
map_dfr(ttSig_bonf, dim)
# # A tibble: 2 x 2
#   DLPFC HIPPO
#   <int> <int>
# 1   448   345
# 2    47    47


map_int(ttSig, ~ length(unique(.x$geneid)))
# DLPFC HIPPO
#  1032   910

map_int(ttSig_bonf, ~ length(unique(.x$geneid)))
# DLPFC HIPPO
#   146   130

map_int(ttSig, ~ length(unique(.x$geneid[.x$TWAS.P < 5e-08])))
# DLPFC HIPPO
#    51    47

dim(tt)
# [1] 132936     47
dim(region_twas_z)
# [1] 101304     16

## save for later
save(tt, ttSig, ttSig_bonf, get_variable_by_region, ttReg_map, region_twas_z, file = 'rda/pgc2_tt_objects.Rdata')


## Continue
load('rda/pgc2_tt_objects.Rdata', verbose = TRUE)




## Check correlations among FDR corrected p-values between TWAS
## and either BEST GWAS FDR or EQTL GWAS FDR
check_cor <- function(x, y) {
    cor(-log10(x), -log10(y))
}


with(tt, check_cor(TWAS.P, BEST.GWAS.P.computed))
# [1] 0.3941553
with(tt, check_cor(TWAS.P, EQTL.P.computed))
# [1] 0.8042339

map_dbl(ttSig, ~ with(.x, check_cor(TWAS.P, BEST.GWAS.P.computed)))
#    DLPFC    HIPPO
# 0.580336 0.544264
map_dbl(ttSig, ~ with(.x, check_cor(TWAS.P, EQTL.P.computed)))
#     DLPFC     HIPPO
# 0.6949813 0.6957536

map_dbl(ttSig_bonf, ~ with(.x, check_cor(TWAS.P, BEST.GWAS.P.computed)))
#     DLPFC     HIPPO
# 0.5443409 0.5150044
map_dbl(ttSig_bonf, ~ with(.x, check_cor(TWAS.P, EQTL.P.computed)))
#     DLPFC     HIPPO
# 0.6294699 0.6180254

map_dbl(split(tt, tt$BEST.GWAS.status), ~ with(.x, check_cor(TWAS.P, BEST.GWAS.P.computed)))
#      Index      Other      Proxy
# -0.1043150  0.3533135  0.2944056
map_dbl(split(tt, tt$BEST.GWAS.status), ~ with(.x, check_cor(TWAS.P, EQTL.P.computed)))
#     Index     Other     Proxy
# 0.7201989 0.7873373 0.8291230


tt_sigonly <- tt[tt$TWAS.FDR < 0.05, ]
with(tt_sigonly, addmargins(table(BEST.GWAS.status, EQTL.status, useNA = 'ifany')))
#                 EQTL.status
# BEST.GWAS.status Index Other Proxy  Sum
#            Index     1   174    50  225
#            Other     0  5306    34 5340
#            Proxy     1   866   366 1233
#            Sum       2  6346   450 6798

tt_sigonly_bonf <- tt[tt$TWAS.Bonf < 0.05, ]
with(tt_sigonly_bonf, addmargins(table(BEST.GWAS.status, EQTL.status, useNA = 'ifany')))
#                 EQTL.status
# BEST.GWAS.status Index Other Proxy Sum
#            Index     1    36     7  44
#            Other     0   339    13 352
#            Proxy     0   217   180 397
#            Sum       1   592   200 793

create_gwas_or_eqtl <- function(tt_sigonly, filename = 'pdf/twas_fdr5perc_vs_gwas_or_eqtl.pdf', titleslug = 'FDR') {
    pdf(filename, useDingbats = FALSE, width = 21, height = 14)
    print(ggplot(tt_sigonly, aes(
        x = -log10(TWAS.P),
        y = -log10(EQTL.P.computed),
        color = BEST.GWAS.P.computed < 5e-08
    )) + geom_point() +
        facet_grid(region * 
            ifelse(EQTL.status == 'Other', 'Other', 'Risk Locus') ~
            factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
        ) +
        theme_bw(base_size = 30) +
        ggtitle(paste0('TWAS (', titleslug, '<5%) vs EQTL p-values')) +
        labs(caption = 'Risk Loci by EQTL')
    )

    print(ggplot(tt_sigonly, aes(
        x = -log10(TWAS.P),
        y = -log10(EQTL.P.computed),
        color = BEST.GWAS.P.computed < 5e-08
    )) + geom_point() +
        facet_grid(region * 
            ifelse(BEST.GWAS.status == 'Other', 'Other', 'Risk Locus') ~
            factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
        ) +
        theme_bw(base_size = 30) +
        ggtitle(paste0('TWAS (', titleslug, '<5%) vs EQTL p-values')) +
        labs(caption = 'Risk Loci by BEST GWAS')
    )
    
    print(ggplot(tt_sigonly, aes(
        x = TWAS.Z,
        y = EQTL.Z,
        color = BEST.GWAS.P.computed < 5e-08
    )) + geom_point() +
        facet_grid(region * 
            ifelse(EQTL.status == 'Other', 'Other', 'Risk Locus') ~
            factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
        ) +
        theme_bw(base_size = 30) +
        ggtitle(paste0('TWAS (', titleslug, '<5%) vs EQTL z-scores')) +
        labs(caption = 'Risk Loci by EQTL')
    )

    print(ggplot(tt_sigonly, aes(
        x = TWAS.Z,
        y = EQTL.Z,
        color = BEST.GWAS.P.computed < 5e-08
    )) + geom_point() +
        facet_grid(region * 
            ifelse(BEST.GWAS.status == 'Other', 'Other', 'Risk Locus') ~
            factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
        ) +
        theme_bw(base_size = 30) +
        ggtitle(paste0('TWAS (', titleslug, '<5%) vs EQTL z-scores')) +
        labs(caption = 'Risk Loci by BEST GWAS')
    )
    
    dev.off()
}
create_gwas_or_eqtl(tt_sigonly, 'pdf/pgc2_twas_fdr5perc_vs_gwas_or_eqtl.pdf', 'FDR')
create_gwas_or_eqtl(tt_sigonly_bonf, 'pdf/pgc2_twas_bonf5perc_vs_gwas_or_eqtl.pdf', 'Bonf')


create_by_status <- function(tt_sigonly, filename = 'pdf/twas_fdr5perc_by_status.pdf', titleslug = 'FDR') {
    pdf(filename, useDingbats = FALSE, width = 28, height = 14)
    print(ggplot(tt_sigonly, aes(
        y = -log10(TWAS.P),
        x = ifelse(EQTL.status == 'Other', 'Other', 'Risk Locus'),
        fill = BEST.GWAS.P.computed < 5e-08
    )) + geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_point(aes(fill = BEST.GWAS.P.computed < 5e-08), shape = 21, position = position_jitterdodge(jitter.width = 0.2)) +
        facet_grid(region ~ factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
        ) +
        theme_bw(base_size = 30) +
        ggtitle(paste0('TWAS (', titleslug, '<5%) by locus')) +
        xlab('Risk Loci assignment by EQTL SNP') +
        ylim(c(0, max(-log10(tt_sigonly$TWAS.P))))
    )

    print(ggplot(tt_sigonly, aes(
        y = -log10(TWAS.P),
        x = ifelse(BEST.GWAS.status == 'Other', 'Other', 'Risk Locus'),
        fill = BEST.GWAS.P.computed < 5e-08
    )) + geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_point(aes(fill = BEST.GWAS.P.computed < 5e-08), shape = 21, position = position_jitterdodge(jitter.width = 0.2)) +
        facet_grid(region ~ factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
        ) +
        theme_bw(base_size = 30) +
        ggtitle(paste0('TWAS (', titleslug, '<5%) by locus')) +
        xlab('Risk Loci assignment by BEST GWAS SNP') +
        ylim(c(0, max(-log10(tt_sigonly$TWAS.P))))
    )
    dev.off()
    
}
create_by_status(tt_sigonly, 'pdf/pgc2_twas_fdr5perc_by_status.pdf', 'FDR')
create_by_status(tt_sigonly_bonf, 'pdf/pgc2_twas_bonf5perc_by_status.pdf', 'Bonf')


map(ttSig, ~ map(split(.x, .x$feature), ~
    addmargins(table(
        'BEST GWAS P < 5e-08' = .x$BEST.GWAS.P.computed < 5e-08,
        'Risk Locus (by BEST GWAS)' =  .x$BEST.GWAS.status != 'Other',
        useNA = 'ifany'
    ))
))
# $DLPFC
# $DLPFC$exon
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE  1292  203 1495
#               TRUE    206  270  476
#               Sum    1498  473 1971
#
# $DLPFC$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   195   20 215
#               TRUE     24   35  59
#               Sum     219   55 274
#
# $DLPFC$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE   768   87  855
#               TRUE    101  138  239
#               Sum     869  225 1094
#
# $DLPFC$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   359   39 398
#               TRUE     35   74 109
#               Sum     394  113 507
#
#
# $HIPPO
# $HIPPO$exon
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE   981  118 1099
#               TRUE    154  181  335
#               Sum    1135  299 1434
#
# $HIPPO$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   119   12 131
#               TRUE     15   20  35
#               Sum     134   32 166
#
# $HIPPO$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   639   69 708
#               TRUE     77  105 182
#               Sum     716  174 890
#
# $HIPPO$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   335   24 359
#               TRUE     40   63 103
#               Sum     375   87 462

map(ttSig_bonf, ~ map(split(.x, .x$feature), ~
    addmargins(table(
        'BEST GWAS P < 5e-08' = .x$BEST.GWAS.P.computed < 5e-08,
        'Risk Locus (by BEST GWAS)' =  .x$BEST.GWAS.status != 'Other',
        useNA = 'ifany'
    ))
))
# $DLPFC
# $DLPFC$exon
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    21   28  49
#               TRUE     53  102 155
#               Sum      74  130 204
#
# $DLPFC$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    17    2  19
#               TRUE      9   20  29
#               Sum      26   22  48
#
# $DLPFC$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    16    9  25
#               TRUE     43   50  93
#               Sum      59   59 118
#
# $DLPFC$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    15    7  22
#               TRUE     22   34  56
#               Sum      37   41  78
#
#
# $HIPPO
# $HIPPO$exon
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    10    9  19
#               TRUE     55   78 133
#               Sum      65   87 152
#
# $HIPPO$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE     8    4  12
#               TRUE      5   14  19
#               Sum      13   18  31
#
# $HIPPO$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    13   14  27
#               TRUE     34   32  66
#               Sum      47   46  93
#
# $HIPPO$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    13    6  19
#               TRUE     18   32  50
#               Sum      31   38  69


## Add locus considered section

## Read in the files that Emily cleaned up at
## https://github.com/LieberInstitute/brainseq_phase2/blob/master/eQTL_GWAS_riskSNPs/create_eqtl_table_indexInfo.R
raggr_clean_files <- c(
    'HIPPO' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/raggr_179_snps_hippo_eqtls_fdr01.csv',
    'DLPFC' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/raggr_179_snps_dlp_eqtls_fdr01.csv'
)
raggr_clean <- map(raggr_clean_files, read.csv, stringsAsFactors = FALSE)
names(raggr_clean) <- names(raggr_clean_files)

check_by_locus <- function(rag, ref) {
    by_loc <- split(rag$SNP, rag$IndexSNP)
    map_dbl(by_loc, ~ sum(.x %in% ref))
}


clean_by_state <- function(x) {
    r <- map_dfr(x, ~ .x)
    r$state <- c(FALSE, TRUE)
    return(r)
}

by_locus <- function(cut, list = FALSE, var = 'TWAS.FDR') {
    by_locus <- map2(
        raggr_clean,
        map(names(raggr_clean), ~ tt$BEST.GWAS.ID[tt[, var, drop = TRUE] < cut & tt$region == .x]),
        check_by_locus
    )
    if(list) return(by_locus)
    clean_by_state(map(by_locus, ~ table(.x > 0)))
}
perc_locus <- function(cut, var = 'TWAS.FDR') {
    x <- by_locus(cut, var = var)
    x[2, 1:2] / colSums(x[, 1:2]) * 100
}

stopifnot(identical(by_locus(1.1), by_locus(1.1, var = 'TWAS.Bonf')))

## raggr eQTL locus considered in the TWAS analysis
by_locus(1.1)
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    53    58 FALSE
# 2    50    58 TRUE
perc_locus(1.1)
#      HIPPO DLPFC
# 1 48.54369    50

## raggr eQTL locus that have a TWAS FDR <5% result
by_locus(0.05)
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    63    67 FALSE
# 2    40    49 TRUE
perc_locus(0.05)
#      HIPPO    DLPFC
# 1 38.83495 42.24138
perc_locus(0.05) / perc_locus(1.1) * 100
#   HIPPO    DLPFC
# 1    80 84.48276

## raggr eQTL locus that have a TWAS Bonf <5% result
by_locus(0.05, var = 'TWAS.Bonf')
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    79    85 FALSE
# 2    24    31 TRUE
perc_locus(0.05, var = 'TWAS.Bonf')
#      HIPPO    DLPFC
# 1 23.30097 26.72414
perc_locus(0.05, var = 'TWAS.Bonf') / perc_locus(1.1) * 100
#   HIPPO    DLPFC
# 1    48 53.44828


get_matrix <- function(x) {
    matrix(x, ncol = ncol(x), dimnames = attr(x, 'dimnames'))
}

get_venn_info <- function(cut, var = 'TWAS.FDR') {
    map(
        by_locus(cut, list = TRUE, var = var),
        ~ names(which(.x > 0))
    )
}
venn_by_locus <- function(cut, var = 'TWAS.FDR') {
    venn(
        get_venn_info(cut, var = var),
        show.plot = FALSE
    )
}

shared_by_locus <- function(cut, var = 'TWAS.FDR') {
    get_matrix(
        venn_by_locus(cut, var = var)
    )
}

## Overlap by region of all loci considered
shared_by_locus(1.1)
#    num HIPPO DLPFC
# 00   0     0     0
# 01  12     0     1
# 10   4     1     0
# 11  46     1     1

## Now with the ones that are TWAS FDR < 5%
shared_by_locus(0.05)
#    num HIPPO DLPFC
# 00   0     0     0
# 01  11     0     1
# 10   2     1     0
# 11  38     1     1

## Now with the ones that are TWAS Bonf < 5%
shared_by_locus(0.05, var = 'TWAS.Bonf')
#    num HIPPO DLPFC
# 00   0     0     0
# 01   9     0     1
# 10   2     1     0
# 11  22     1     1

make_pretty_venn <- function(cut, title = '', var = 'TWAS.FDR') {
    info <- get_venn_info(cut, var = var)
    cols <- c('DLPFC' = 'dark orange', 'HIPPO' = 'skyblue3')
    v <- venn.diagram(info, filename = NULL,
        main = title,
        col = 'transparent', fill = rev(cols),
        alpha = 0.5, margin = 0,
        main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
        cat.col = rev(cols))
    grid.newpage()
    grid.draw(v)
}

pdf('pdf/pgc2_venn_by_locus.pdf', useDingbats = FALSE)
make_pretty_venn(1.1, 'rAggr loci considered in TWAS')
make_pretty_venn(0.05, 'rAggr loci with TWAS FDR<5%')
make_pretty_venn(0.05, 'rAggr loci with TWAS Bonf<5%', var = 'TWAS.Bonf')
dev.off()

system('rm VennDiagram*')


gene_by_locus <- function(rag, ref) {
    by_loc <- split(rag$gene, rag$IndexSNP)
    map_dbl(by_loc, ~ sum(.x %in% ref))
}

features <- c('gene', 'exon', 'jxn', 'tx')

by_feature <- function(cut, list = FALSE, var = 'TWAS.FDR') {
    g_by_locus <- map(features, function(feature) {
        map2(
            map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
            map(names(raggr_clean), ~ tt$ID[tt$feature == feature & tt$region == .x & tt[, var, drop = TRUE] < cut]),
            gene_by_locus
        )
    })
    names(g_by_locus) <- features
    if(list) return(g_by_locus)
    map(g_by_locus, function(x) {
        r <- map_dfr(x, ~ table(.x > 0))
        r$state <- c(FALSE, TRUE)
        return(r)
    })
}

clean_tabs <- function(l) {
    map2_dfr(l, names(l), function(x, y) {
        x$feature <- y
        return(x)
    })
}

## Features by locus that were considered in the TWAS analysis
clean_tabs(by_feature(cut = 1.1))
## Numbers differ from https://github.com/LieberInstitute/brainseq_phase2/blob/master/twas/explore_twas.R#L774-L784
## because in this script we already dropped the results with NA TWAS.P values
# # A tibble: 8 x 4
#   HIPPO DLPFC state feature
#   <int> <int> <lgl> <chr>
# 1    16    16 FALSE gene
# 2    34    51 TRUE  gene
# 3    16    18 FALSE exon
# 4    61    69 TRUE  exon
# 5    22    24 FALSE jxn
# 6    66    74 TRUE  jxn
# 7    19    17 FALSE tx
# 8    46    59 TRUE  tx

perc_feature <- function(cut, var = 'TWAS.FDR') {
    y <- clean_tabs(by_feature(cut, var = var))
    clean_tabs(
        map(
            split(y, factor(y$feature, levels = features)), 
            ~ .x[2, 1:2] / colSums(.x[, 1:2]) * 100
    ))
}

perc_feature(1.1)
#      HIPPO    DLPFC feature
# 1 68.00000 76.11940    gene
# 2 79.22078 79.31034    exon
# 3 75.00000 75.51020     jxn
# 4 70.76923 77.63158      tx


## Features by locus that have a TWAS FDR<5% result
clean_tabs(by_feature(cut = 0.05))
# # A tibble: 8 x 4
#   HIPPO DLPFC state feature
#   <int> <int> <lgl> <chr>
# 1    29    34 FALSE gene
# 2    21    33 TRUE  gene
# 3    31    31 FALSE exon
# 4    46    56 TRUE  exon
# 5    39    35 FALSE jxn
# 6    49    63 TRUE  jxn
# 7    32    30 FALSE tx
# 8    33    46 TRUE  tx

perc_feature(0.05)
#      HIPPO    DLPFC feature
# 1 42.00000 49.25373    gene
# 2 59.74026 64.36782    exon
# 3 55.68182 64.28571     jxn
# 4 50.76923 60.52632      tx

## Compute the percent using as denominator the number of features considered
cbind(perc_feature(0.05)[, 1:2] / perc_feature(1.1)[, 1:2] * 100, feature = features)
#      HIPPO    DLPFC feature
# 1 61.76471 64.70588    gene
# 2 75.40984 81.15942    exon
# 3 74.24242 85.13514     jxn
# 4 71.73913 77.96610      tx

## Now with Bonf < 5%
clean_tabs(by_feature(cut = 0.05, var = 'TWAS.Bonf'))
# # A tibble: 8 x 4
#   HIPPO DLPFC state feature
#   <int> <int> <lgl> <chr>
# 1    34    49 FALSE gene
# 2    16    18 TRUE  gene
# 3    56    59 FALSE exon
# 4    21    28 TRUE  exon
# 5    67    68 FALSE jxn
# 6    21    30 TRUE  jxn
# 7    45    53 FALSE tx
# 8    20    23 TRUE  tx

perc_feature(0.05, var = 'TWAS.Bonf')
#      HIPPO    DLPFC feature
# 1 32.00000 26.86567    gene
# 2 27.27273 32.18391    exon
# 3 23.86364 30.61224     jxn
# 4 30.76923 30.26316      tx

cbind(perc_feature(0.05, var = 'TWAS.Bonf')[, 1:2] / perc_feature(1.1)[, 1:2] * 100, feature = features)
#      HIPPO    DLPFC feature
# 1 47.05882 35.29412    gene
# 2 34.42623 40.57971    exon
# 3 31.81818 40.54054     jxn
# 4 43.47826 38.98305      tx

get_venn_info_by_feature <- function(cut, var = 'TWAS.FDR') {
    map(
        by_feature(cut, list = TRUE, var = var),
        ~ map(.x, 
            ~ names(which(.x > 0))
        )
    )
}
venn_by_locus_by_feature <- function(cut, var = 'TWAS.FDR') {
    map(get_venn_info_by_feature(cut, var = var), venn, show.plot = FALSE)
}

shared_by_locus_by_feature <- function(cut, var = 'TWAS.FDR') {
    map(venn_by_locus_by_feature(cut, var = var), get_matrix)
}

## Find the locus that have shared features
## between rAggr and TWAS (given a FDR cut)
## then find the locus overlap across regions.
shared_by_locus_by_feature(1.1)
# $gene
#    num HIPPO DLPFC
# 00   0     0     0
# 01  22     0     1
# 10   5     1     0
# 11  29     1     1
#
# $exon
#    num HIPPO DLPFC
# 00   0     0     0
# 01  15     0     1
# 10   7     1     0
# 11  54     1     1
#
# $jxn
#    num HIPPO DLPFC
# 00   0     0     0
# 01  16     0     1
# 10   8     1     0
# 11  58     1     1
#
# $tx
#    num HIPPO DLPFC
# 00   0     0     0
# 01  22     0     1
# 10   9     1     0
# 11  37     1     1

## TWAS FDR <5%
shared_by_locus_by_feature(0.05)
# $gene
#    num HIPPO DLPFC
# 00   0     0     0
# 01  16     0     1
# 10   4     1     0
# 11  17     1     1
#
# $exon
#    num HIPPO DLPFC
# 00   0     0     0
# 01  15     0     1
# 10   5     1     0
# 11  41     1     1
#
# $jxn
#    num HIPPO DLPFC
# 00   0     0     0
# 01  20     0     1
# 10   6     1     0
# 11  43     1     1
#
# $tx
#    num HIPPO DLPFC
# 00   0     0     0
# 01  19     0     1
# 10   6     1     0
# 11  27     1     1

## TWAS Bonf <5%
shared_by_locus_by_feature(0.05, var = 'TWAS.Bonf')
# $gene
#    num HIPPO DLPFC
# 00   0     0     0
# 01   7     0     1
# 10   5     1     0
# 11  11     1     1
#
# $exon
#    num HIPPO DLPFC
# 00   0     0     0
# 01  10     0     1
# 10   3     1     0
# 11  18     1     1
#
# $jxn
#    num HIPPO DLPFC
# 00   0     0     0
# 01  12     0     1
# 10   3     1     0
# 11  18     1     1
#
# $tx
#    num HIPPO DLPFC
# 00   0     0     0
# 01   7     0     1
# 10   4     1     0
# 11  16     1     1


make_pretty_venn_by_feature <- function(cut, title = '', var = 'TWAS.FDR') {
    info_all <- get_venn_info_by_feature(cut, var = var)
    cols <- c('DLPFC' = 'dark orange', 'HIPPO' = 'skyblue3')
    map2(
        info_all,
        names(info_all),
        function(info, feature) {
        v <- venn.diagram(info, filename = NULL,
            main = paste0(title, ' - ', feature),
            col = 'transparent', fill = rev(cols),
            alpha = 0.5, margin = 0,
            main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
            cat.col = rev(cols))
        grid.newpage()
        grid.draw(v)
        }
    )
}

pdf('pdf/pgc2_venn_by_locus_by_feature.pdf', useDingbats = FALSE)
make_pretty_venn_by_feature(1.1, 'rAggr loci considered in TWAS')
make_pretty_venn_by_feature(0.05, 'rAggr loci with TWAS FDR<5%')
make_pretty_venn_by_feature(0.05, 'rAggr loci with TWAS Bonf<5%', var = 'TWAS.Bonf')
dev.off()

system('rm VennDiagram*')


## Venn diagrams of features by region, then joint (grouped by gene id)

## Number of features that have TWAS weights
n_feat_considered <- cbind(map_dfr(split(tt, tt$region), ~ 
    map_dbl(
        split(.x, factor(.x$feature, levels = features)),
        ~ nrow(.x)
    )
), features)
n_feat_considered
#   DLPFC HIPPO features
# 1  5555  4030     gene
# 2 39725 29866     exon
# 3 20920 16362      jxn
# 4  9206  7272       tx

## Percent from total number of features considered
# $ grep -A 1 "Final RSE" logs/compute_weights_DLPFC_*
# compute_weights_DLPFC_exon.txt:[1] "Final RSE feature dimensions:"
# compute_weights_DLPFC_exon.txt-[1] 381196    397
# --
# compute_weights_DLPFC_gene.txt:[1] "Final RSE feature dimensions:"
# compute_weights_DLPFC_gene.txt-[1] 23402   397
# --
# compute_weights_DLPFC_jxn.txt:[1] "Final RSE feature dimensions:"
# compute_weights_DLPFC_jxn.txt-[1] 269261    397
# --
# compute_weights_DLPFC_tx.txt:[1] "Final RSE feature dimensions:"
# compute_weights_DLPFC_tx.txt-[1] 88969   397

## In percent of features passed to the weight computation step
n_feat_total <- c('gene' = 23402, 'exon' = 381196, 'jxn' = 269261, 'tx' = 88969)
cbind(map_dfr(n_feat_considered[, 1:2], ~ .x / n_feat_total * 100), features)
#       DLPFC     HIPPO features
# 1 23.737287 17.220750     gene
# 2 10.421148  7.834815     exon
# 3  7.769413  6.076632      jxn
# 4 10.347424  8.173634       tx

## Number of features with TWAS FDR<5%
n_feat_twas5perc <- cbind(map_dfr(ttSig, ~ 
    map_dbl(
        split(.x, factor(.x$feature, levels = features)),
        ~ nrow(.x)
    )
), features)
n_feat_twas5perc
#   DLPFC HIPPO features
# 1   274   166     gene
# 2  1971  1434     exon
# 3  1094   890      jxn
# 4   507   462       tx

cbind(n_feat_twas5perc[, 1:2] / n_feat_considered[, 1:2] * 100, features)
#      DLPFC    HIPPO features
# 1 4.932493 4.119107     gene
# 2 4.961611 4.801446     exon
# 3 5.229446 5.439433      jxn
# 4 5.507278 6.353135       tx

## Now with TWAS Bonf <5%
n_feat_twas5perc_bonf <- cbind(map_dfr(ttSig_bonf, ~ 
    map_dbl(
        split(.x, factor(.x$feature, levels = features)),
        ~ nrow(.x)
    )
), features)
n_feat_twas5perc_bonf
#   DLPFC HIPPO features
# 1    48    31     gene
# 2   204   152     exon
# 3   118    93      jxn
# 4    78    69       tx

cbind(n_feat_twas5perc_bonf[, 1:2] / n_feat_considered[, 1:2] * 100, features)
#       DLPFC     HIPPO features
# 1 0.8640864 0.7692308     gene
# 2 0.5135305 0.5089399     exon
# 3 0.5640535 0.5683902      jxn
# 4 0.8472735 0.9488449       tx

get_feat <- function(cut, var = 'TWAS.FDR') {
    x <- tt[ tt[, var, drop = TRUE]< cut, ]
    map(
        split(x, factor(x$feature, levels = features)),
        ~ map(split(.x, .x$region), ~.x$ID)
    )
}

shared_by_feature <- function(cut, var = 'TWAS.FDR') {
    map(
        get_feat(cut, var = var),
        ~ get_matrix(venn(.x, show.plot = FALSE))
    )
}

shared_by_feature(1.1)
# $gene
#     num DLPFC HIPPO
# 00    0     0     0
# 01 1207     0     1
# 10 2732     1     0
# 11 2823     1     1
#
# $exon
#      num DLPFC HIPPO
# 00     0     0     0
# 01 13806     0     1
# 10 23665     1     0
# 11 16060     1     1
#
# $jxn
#      num DLPFC HIPPO
# 00     0     0     0
# 01  7802     0     1
# 10 12360     1     0
# 11  8560     1     1
#
# $tx
#     num DLPFC HIPPO
# 00    0     0     0
# 01 3083     0     1
# 10 5017     1     0
# 11 4189     1     1

## TWAS FDR <5%
shared_by_feature(0.05)
# $gene
#    num DLPFC HIPPO
# 00   0     0     0
# 01  76     0     1
# 10 184     1     0
# 11  90     1     1
#
# $exon
#     num DLPFC HIPPO
# 00    0     0     0
# 01  856     0     1
# 10 1393     1     0
# 11  578     1     1
#
# $jxn
#    num DLPFC HIPPO
# 00   0     0     0
# 01 503     0     1
# 10 707     1     0
# 11 387     1     1
#
# $tx
#    num DLPFC HIPPO
# 00   0     0     0
# 01 270     0     1
# 10 315     1     0
# 11 192     1     1

## TWAS Bonf <5%
shared_by_feature(0.05, var = 'TWAS.Bonf')
# $gene
#    num DLPFC HIPPO
# 00   0     0     0
# 01  17     0     1
# 10  34     1     0
# 11  14     1     1
#
# $exon
#    num DLPFC HIPPO
# 00   0     0     0
# 01 104     0     1
# 10 156     1     0
# 11  48     1     1
#
# $jxn
#    num DLPFC HIPPO
# 00   0     0     0
# 01  61     0     1
# 10  86     1     0
# 11  32     1     1
#
# $tx
#    num DLPFC HIPPO
# 00   0     0     0
# 01  42     0     1
# 10  51     1     0
# 11  27     1     1

make_pretty_venn_shared_by_feature <- function(cut, title = '', var = 'TWAS.FDR') {
    info_all <- get_feat(cut, var = var)
    cols <- c('DLPFC' = 'dark orange', 'HIPPO' = 'skyblue3')
    map2(
        info_all,
        names(info_all),
        function(info, feature) {
        v <- venn.diagram(info, filename = NULL,
            main = paste0(title, ' - ', feature),
            col = 'transparent', fill = cols,
            alpha = 0.5, margin = 0,
            main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
            cat.col = cols)
        grid.newpage()
        grid.draw(v)
        }
    )
}

pdf('pdf/pgc2_venn_by_feature.pdf', useDingbats = FALSE)
make_pretty_venn_shared_by_feature(1.1, 'Features with TWAS weights')
make_pretty_venn_shared_by_feature(0.05, 'Features with TWAS FDR<5%')
make_pretty_venn_shared_by_feature(0.05, 'Features with TWAS Bonf<5%', var = 'TWAS.Bonf')
dev.off()

system('rm VennDiagram*')


## Now by gene id
get_feat_geneid <- function(cut, var = 'TWAS.FDR') {
    x <- tt[ tt[, var, drop = TRUE]< cut, ]
    map(
        split(x, factor(x$feature, levels = features)),
        ~ map(split(.x, .x$region), ~ .x$geneid[!is.na(.x$geneid)])
    )
}

shared_by_geneid <- function(cut, var = 'TWAS.FDR') {
    map(
        get_feat_geneid(cut, var = var),
        ~ get_matrix(venn(.x, show.plot = FALSE))
    )
}

shared_by_geneid(1.1)
# $gene
#     num DLPFC HIPPO
# 00    0     0     0
# 01 1207     0     1
# 10 2732     1     0
# 11 2823     1     1
#
# $exon
#      num DLPFC HIPPO
# 00     0     0     0
# 01  3068     0     1
# 10  6821     1     0
# 11 32904     1     1
#
# $jxn
#      num DLPFC HIPPO
# 00     0     0     0
# 01  1865     0     1
# 10  3409     1     0
# 11 13342     1     1
#
# $tx
#     num DLPFC HIPPO
# 00    0     0     0
# 01 1688     0     1
# 10 3115     1     0
# 11 6091     1     1

## TWAS FDR <5%
shared_by_geneid(0.05)
# $gene
#    num DLPFC HIPPO
# 00   0     0     0
# 01  76     0     1
# 10 184     1     0
# 11  90     1     1
#
# $exon
#     num DLPFC HIPPO
# 00    0     0     0
# 01  392     0     1
# 10  723     1     0
# 11 1248     1     1
#
# $jxn
#    num DLPFC HIPPO
# 00   0     0     0
# 01 201     0     1
# 10 319     1     0
# 11 576     1     1
#
# $tx
#    num DLPFC HIPPO
# 00   0     0     0
# 01 205     0     1
# 10 244     1     0
# 11 263     1     1

## TWAS Bonf <5%
shared_by_geneid(0.05, var = 'TWAS.Bonf')
# $gene
#    num DLPFC HIPPO
# 00   0     0     0
# 01  17     0     1
# 10  34     1     0
# 11  14     1     1
#
# $exon
#    num DLPFC HIPPO
# 00   0     0     0
# 01  30     0     1
# 10  85     1     0
# 11 119     1     1
#
# $jxn
#    num DLPFC HIPPO
# 00   0     0     0
# 01  32     0     1
# 10  44     1     0
# 11  54     1     1
#
# $tx
#    num DLPFC HIPPO
# 00   0     0     0
# 01  30     0     1
# 10  42     1     0
# 11  36     1     1

make_pretty_venn_shared_by_geneid <- function(cut, title = '', var = 'TWAS.FDR') {
    info_all <- get_feat_geneid(cut, var = var)
    cols <- c('DLPFC' = 'dark orange', 'HIPPO' = 'skyblue3')
    map2(
        info_all,
        names(info_all),
        function(info, feature) {
        v <- venn.diagram(info, filename = NULL,
            main = paste0(title, ' - ', feature),
            col = 'transparent', fill = cols,
            alpha = 0.5, margin = 0,
            main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
            cat.col = cols)
        grid.newpage()
        grid.draw(v)
        }
    )
}

pdf('pdf/pgc2_venn_by_feature_using_geneid.pdf', useDingbats = FALSE)
make_pretty_venn_shared_by_geneid(1.1, 'Features with TWAS weights (by gene ID)')
make_pretty_venn_shared_by_geneid(0.05, 'Features with TWAS FDR<5% (by gene ID)')
make_pretty_venn_shared_by_geneid(0.05, 'Features with TWAS Bonf<5% (by gene ID)', var = 'TWAS.Bonf')
dev.off()

system('rm VennDiagram*')


get_feat_geneid2 <- function(cut, var = 'TWAS.FDR') {
    x <- tt[ tt[, var, drop = TRUE]< cut, ]
    map(
        split(x, x$region),
        ~ map(split(.x, factor(.x$feature, levels = features)), ~ .x$geneid[!is.na(.x$geneid)])
    )
}

make_pretty_venn_shared_by_geneid2 <- function(cut, title = '', var = 'TWAS.FDR') {
    info_all <- get_feat_geneid2(cut, var = var)
    cols <- brewer.pal('Set1', n = 4)
    map2(
        info_all,
        names(info_all),
        function(info, feature) {
        v <- venn.diagram(info, filename = NULL,
            main = paste0(title, ' - ', feature),
            col = 'transparent', fill = cols,
            alpha = 0.5, margin = 0,
            main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
            cat.col = cols)
        grid.newpage()
        grid.draw(v)
        }
    )
}

pdf('pdf/pgc2_venn_by_feature_using_geneid_across_features.pdf', useDingbats = FALSE)
make_pretty_venn_shared_by_geneid2(1.1, 'Features with TWAS weights (by gene ID)')
make_pretty_venn_shared_by_geneid2(0.05, 'Features with TWAS FDR<5% (by gene ID)')
make_pretty_venn_shared_by_geneid2(0.05, 'Features with TWAS Bonf<5% (by gene ID)', var = 'TWAS.Bonf')
dev.off()

system('rm VennDiagram*')

cbind(map_dfr(split(tt, tt$region), ~ map_dfr(split(.x, .x$feature), ~ sum(.x$TWAS.FDR < 0.05))), region = c('DLPFC', 'HIPPO'))
#   exon gene  jxn  tx region
# 1 1971  274 1094 507  DLPFC
# 2 1434  166  890 462  HIPPO

cbind(map_dfr(split(tt, tt$region), ~ map_dfr(split(.x, .x$feature), ~ sum(.x$TWAS.P < 5e-08))), region = c('DLPFC', 'HIPPO'))
#   exon gene jxn tx region
# 1   81   10  51 23  DLPFC
# 2   65   11  27 17  HIPPO

cbind(map_dfr(split(tt, tt$region), ~ map_dfr(split(.x, .x$feature), ~ sum(p.adjust(.x$TWAS.P, 'bonf') < 0.05))), region = c('DLPFC', 'HIPPO'))
#   exon gene jxn tx region
# 1  204   48 118 78  DLPFC
# 2  152   31  93 69  HIPPO



## Compare against https://www.nature.com/articles/s41588-018-0092-1#Sec27
gusev_twas <- read_xlsx('41588_2018_92_MOESM3_ESM.xlsx', sheet = 1)
gusev_twas <- subset(gusev_twas, Expression %in% c('CMC', 'CMC-splicing'))
dim(gusev_twas)
# [1] 124  13

length(unique(gusev_twas$Gene))
# [1] 83

gusev_twas_byset <- split(gusev_twas, gusev_twas$Expression)
map_dbl(gusev_twas_byset, ~ length(unique(.x$Gene)))
# CMC CMC-splicing
#  44           46

gusev_gene <- data.frame(
    gene = unique(gusev_twas$Gene),
    stringsAsFactors = FALSE
)
## Categorize the genes by whether they show up in both CMC and CMC-splicing or only one
gusev_gene$set_evidence <- map_chr(gusev_gene$gene, function(g) {
    res <- map_lgl(gusev_twas_byset, ~ g %in% .x$Gene)
    if(all(res)) {
        return('Both')
    } else {
        names(gusev_twas_byset)[res]
    }
})
table(gusev_gene$set_evidence)
# Both          CMC CMC-splicing
#    7           37           39

twas_ov <- function(sig, prefix) {
    genes <- tolower(gusev_gene$gene)
    res <- map2_dfc(sig, names(ttSig), function(tsub, region) {
        res2 <- map_dfc(split(tsub, factor(tsub$feature, levels = features)), ~ genes %in% tolower(.x$genesymbol))
        res2$any_feature <- pmap_lgl(res2, any)
        colnames(res2) <- paste0(region, '_', colnames(res2))
        return(res2)
    })
    colnames(res) <- paste0(prefix, colnames(res))
    return(res)
}

## Load the CLOZUK+PGC2 TWAS results
get_psycm <- function() {
    e <- new.env()
    load('rda/tt_objects.Rdata', verbose = TRUE, envir = e)
    e
}

psycm <- get_psycm()
stopifnot(!identical(ttSig, psycm$ttSig))
stopifnot(!identical(ttSig_bonf, psycm$ttSig_bonf))


## Find whether the Gusev et al TWAS genes show up in our data
gusev_gene <- cbind(
    gusev_gene,
    twas_ov(ttSig, 'pgc2_FDR_'),
    twas_ov(ttSig_bonf, 'pgc2_Bonf_'),
    twas_ov(psycm$ttSig, 'psycm_FDR_'),
    twas_ov(psycm$ttSig_bonf, 'psycm_Bonf_')
)

## Is the gene in any of our results?
gusev_gene$in_any <- pmap_lgl(gusev_gene[, -(1:2)], any)
gusev_gene$in_any_FDR <- pmap_lgl(gusev_gene[, grep('FDR', colnames(gusev_gene))], any)
gusev_gene$in_any_Bonf <- pmap_lgl(gusev_gene[, grep('Bonf', colnames(gusev_gene))], any)

dim(gusev_gene)
# [1] 83 45

## Run a quick check
stopifnot(all(with(gusev_gene,
    pgc2_FDR_DLPFC_gene | pgc2_FDR_DLPFC_exon | pgc2_FDR_DLPFC_jxn | pgc2_FDR_DLPFC_tx == pgc2_FDR_DLPFC_any_feature
)))

save(gusev_gene, file = 'rda/gusev_gene.Rdata')

## Birds eye view across all features and all TWAS we did
with(gusev_gene, addmargins(table(set_evidence, in_any)))
#               in_any
# set_evidence   FALSE TRUE Sum
#   Both             0    7   7
#   CMC             10   27  37
#   CMC-splicing     3   36  39
#   Sum             13   70  83

## Now using TWAS FDR<5% only (exact same numbers as above)
with(gusev_gene, addmargins(table(set_evidence, in_any_FDR)))
#               in_any_FDR
# set_evidence   FALSE TRUE Sum
#   Both             0    7   7
#   CMC             10   27  37
#   CMC-splicing     3   36  39
#   Sum             13   70  83

## Or TWAS Bonf <5% only (more restrictive)
with(gusev_gene, addmargins(table(set_evidence, in_any_Bonf)))
#               in_any_Bonf
# set_evidence   FALSE TRUE Sum
#   Both             1    6   7
#   CMC             16   21  37
#   CMC-splicing     9   30  39
#   Sum             26   57  83

## Main motivation behind this comparison:
## DLPFC with PGC2 GWAS vs Gusev et al
## Either with TWAS FDR <5%
with(gusev_gene, addmargins(table(set_evidence, pgc2_FDR_DLPFC_any_feature)))
# set_evidence   FALSE TRUE Sum
#   Both             1    6   7
#   CMC             12   25  37
#   CMC-splicing     6   33  39
#   Sum             19   64  83

## by feature
map_int(gusev_gene[, grep('pgc2_FDR_DLPFC', colnames(gusev_gene))], ~ sum(.x))
#        pgc2_FDR_DLPFC_gene        pgc2_FDR_DLPFC_exon
#                         23                         45
#         pgc2_FDR_DLPFC_jxn          pgc2_FDR_DLPFC_tx
#                         41                         43
# pgc2_FDR_DLPFC_any_feature
#                         64

## or detailed comparison by feature
map(gusev_gene[, grep('pgc2_FDR_DLPFC', colnames(gusev_gene))], ~ addmargins(table(set_evidence = gusev_gene$set_evidence, 'present?' = .x)))
# $pgc2_FDR_DLPFC_gene
#               present?
# set_evidence   FALSE TRUE Sum
#   Both             6    1   7
#   CMC             21   16  37
#   CMC-splicing    33    6  39
#   Sum             60   23  83
#
# $pgc2_FDR_DLPFC_exon
#               present?
# set_evidence   FALSE TRUE Sum
#   Both             2    5   7
#   CMC             16   21  37
#   CMC-splicing    20   19  39
#   Sum             38   45  83
#
# $pgc2_FDR_DLPFC_jxn
#               present?
# set_evidence   FALSE TRUE Sum
#   Both             2    5   7
#   CMC             26   11  37
#   CMC-splicing    14   25  39
#   Sum             42   41  83
#
# $pgc2_FDR_DLPFC_tx
#               present?
# set_evidence   FALSE TRUE Sum
#   Both             1    6   7
#   CMC             19   18  37
#   CMC-splicing    20   19  39
#   Sum             40   43  83
#
# $pgc2_FDR_DLPFC_any_feature
#               present?
# set_evidence   FALSE TRUE Sum
#   Both             1    6   7
#   CMC             12   25  37
#   CMC-splicing     6   33  39
#   Sum             19   64  83

                         
## or with TWAS Bonf <5%
with(gusev_gene, addmargins(table(set_evidence, pgc2_Bonf_DLPFC_any_feature)))
# set_evidence   FALSE TRUE Sum
#   Both             3    4   7
#   CMC             28    9  37
#   CMC-splicing    21   18  39
#   Sum             52   31  83

## by feature
map_int(gusev_gene[, grep('pgc2_Bonf_DLPFC', colnames(gusev_gene))], ~ sum(.x))
#        pgc2_Bonf_DLPFC_gene        pgc2_Bonf_DLPFC_exon
#                           9                          19
#         pgc2_Bonf_DLPFC_jxn          pgc2_Bonf_DLPFC_tx
#                          19                          18
# pgc2_Bonf_DLPFC_any_feature
#                          31

## or detailed comparison by feature
map(gusev_gene[, grep('pgc2_Bonf_DLPFC', colnames(gusev_gene))], ~ addmargins(table(set_evidence = gusev_gene$set_evidence, 'present?' = .x)))
# $pgc2_Bonf_DLPFC_gene
#               present?
# set_evidence   FALSE TRUE Sum
#   Both             6    1   7
#   CMC             32    5  37
#   CMC-splicing    36    3  39
#   Sum             74    9  83
#
# $pgc2_Bonf_DLPFC_exon
#               present?
# set_evidence   FALSE TRUE Sum
#   Both             4    3   7
#   CMC             29    8  37
#   CMC-splicing    31    8  39
#   Sum             64   19  83
#
# $pgc2_Bonf_DLPFC_jxn
#               present?
# set_evidence   FALSE TRUE Sum
#   Both             3    4   7
#   CMC             33    4  37
#   CMC-splicing    28   11  39
#   Sum             64   19  83
#
# $pgc2_Bonf_DLPFC_tx
#               present?
# set_evidence   FALSE TRUE Sum
#   Both             5    2   7
#   CMC             30    7  37
#   CMC-splicing    30    9  39
#   Sum             65   18  83
#
# $pgc2_Bonf_DLPFC_any_feature
#               present?
# set_evidence   FALSE TRUE Sum
#   Both             3    4   7
#   CMC             28    9  37
#   CMC-splicing    21   18  39
#   Sum             52   31  83


## Per each TWAS we made, what's the overlap across regions?
## Start with PGC2
## TWAS FDR <5%
with(gusev_gene, addmargins(table(pgc2_FDR_DLPFC_any_feature, pgc2_FDR_HIPPO_any_feature)))
#                           pgc2_FDR_HIPPO_any_feature
# pgc2_FDR_DLPFC_any_feature FALSE TRUE Sum
#                      FALSE    16    3  19
#                      TRUE     16   48  64
#                      Sum      32   51  83

## TWAS Bonf <5%
with(gusev_gene, addmargins(table(pgc2_Bonf_DLPFC_any_feature, pgc2_Bonf_HIPPO_any_feature)))
#                            pgc2_Bonf_HIPPO_any_feature
# pgc2_Bonf_DLPFC_any_feature FALSE TRUE Sum
#                       FALSE    47    5  52
#                       TRUE      9   22  31
#                       Sum      56   27  83

## Start with CLOZUK+PCG
## TWAS FDR <5%
with(gusev_gene, addmargins(table(psycm_FDR_DLPFC_any_feature, psycm_FDR_HIPPO_any_feature)))
#                            psycm_FDR_HIPPO_any_feature
# psycm_FDR_DLPFC_any_feature FALSE TRUE Sum
#                       FALSE    13    4  17
#                       TRUE     15   51  66
#                       Sum      28   55  83
                      
with(gusev_gene, addmargins(table(psycm_Bonf_DLPFC_any_feature, psycm_Bonf_HIPPO_any_feature)))
#                             psycm_Bonf_HIPPO_any_feature
# psycm_Bonf_DLPFC_any_feature FALSE TRUE Sum
#                        FALSE    27    5  32
#                        TRUE     20   31  51
#                        Sum      47   36  83

## Compare by PGC2 and CLOZUK+PGC2 TWAS by merging both brain regions
## TWAS FDR <5%
with(gusev_gene, addmargins(table(
    'PGC2' = pgc2_FDR_DLPFC_any_feature | pgc2_FDR_HIPPO_any_feature,
    'CLOZUK+PGC2' = psycm_FDR_DLPFC_any_feature | psycm_FDR_HIPPO_any_feature
)))
#        CLOZUK+PGC2
# PGC2    FALSE TRUE Sum
#   FALSE    13    3  16
#   TRUE      0   67  67
#   Sum      13   70  83

## TWAS Bonf <5%
with(gusev_gene, addmargins(table(
    'PGC2' = pgc2_Bonf_DLPFC_any_feature | pgc2_Bonf_HIPPO_any_feature,
    'CLOZUK+PGC2' = psycm_Bonf_DLPFC_any_feature | psycm_Bonf_HIPPO_any_feature
)))
#        CLOZUK+PGC2
# PGC2    FALSE TRUE Sum
#   FALSE    26   21  47
#   TRUE      1   35  36
#   Sum      27   56  83


## Compare SCZD t vs TWAS z
pdf('pdf/pgc2_sczd_t_vs_twas_z.pdf', useDingbats = FALSE, width = 21, height = 18)

tt$status <- ifelse(tt$BEST.GWAS.status == 'Other', 'Other', 'Risk Locus')

## Used https://stackoverflow.com/questions/11889625/annotating-text-on-individual-facet-in-ggplot2
dat_text <- map_dfr(features, function(feat) {
    map_dfr(names(outFeat), function(region) {
        map_dfr(c('Other', 'Risk Locus'), function(state) {
            
            i <- which(tt$feature == feat & tt$region == region & tt$status == state)
            data.frame(
                region = region,
                feature = feat,
                status = state,
                cor = signif(cor(tt$TWAS.Z[i], tt$SCZD_t[i]), 3),
                stringsAsFactors = FALSE
            )
        })
    })
})

## Used https://stackoverflow.com/questions/39623636/forcing-r-output-to-be-scientific-notation-with-at-most-two-decimals/39625148
ggplot(tt, aes(x = TWAS.Z, y = SCZD_t, color = BEST.GWAS.P.computed < 5e-08)) +
    geom_point() +
    facet_grid(region * 
        status ~
        factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
    ) +
    theme_bw(base_size = 30) +
    ggtitle('TWAS vs SCZD differential expression') +
    xlab('TWAS Z score') +
    ylab('SCZD vs control t-statistic') +
    labs(caption = 'Risk Loci by BEST GWAS') +
    guides(color=guide_legend(title="BEST GWAS\np < 5e-08")) +
    geom_text(
      data    = dat_text,
      color = 'black',
      size = 7,
      mapping = aes(x = -4, y = 5.5, label = paste0('rho=', formatC(cor, format = "e", digits = 2))),
    )
dev.off()


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
#  date     2019-03-21
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package        * version  date       lib source
#  assertthat       0.2.0    2017-04-11 [2] CRAN (R 3.5.0)
#  BiocGenerics   * 0.28.0   2018-10-30 [1] Bioconductor
#  bitops           1.0-6    2013-08-17 [2] CRAN (R 3.5.0)
#  caTools          1.17.1.2 2019-03-06 [2] CRAN (R 3.5.1)
#  cellranger       1.1.0    2016-07-27 [1] CRAN (R 3.5.0)
#  cli              1.0.1    2018-09-25 [1] CRAN (R 3.5.1)
#  colorout       * 1.2-0    2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace       1.4-0    2019-01-13 [2] CRAN (R 3.5.1)
#  crayon           1.3.4    2017-09-16 [1] CRAN (R 3.5.0)
#  digest           0.6.18   2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr            0.8.0.1  2019-02-15 [1] CRAN (R 3.5.1)
#  fansi            0.4.0    2018-10-05 [1] CRAN (R 3.5.1)
#  formatR          1.6      2019-03-05 [1] CRAN (R 3.5.1)
#  futile.logger  * 1.4.3    2016-07-10 [1] CRAN (R 3.5.0)
#  futile.options   1.0.1    2018-04-20 [2] CRAN (R 3.5.0)
#  gdata            2.18.0   2017-06-06 [2] CRAN (R 3.5.0)
#  ggplot2        * 3.1.0    2018-10-25 [1] CRAN (R 3.5.1)
#  glue             1.3.1    2019-03-12 [1] CRAN (R 3.5.1)
#  gplots         * 3.0.1.1  2019-01-27 [1] CRAN (R 3.5.1)
#  gtable           0.2.0    2016-02-26 [2] CRAN (R 3.5.0)
#  gtools           3.8.1    2018-06-26 [2] CRAN (R 3.5.1)
#  htmltools        0.3.6    2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets      1.3      2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv           1.4.5.1  2018-12-18 [2] CRAN (R 3.5.1)
#  jsonlite         1.6      2018-12-07 [2] CRAN (R 3.5.1)
#  KernSmooth       2.23-15  2015-06-29 [3] CRAN (R 3.5.1)
#  labeling         0.3      2014-08-23 [2] CRAN (R 3.5.0)
#  lambda.r         1.2.3    2018-05-17 [1] CRAN (R 3.5.0)
#  later            0.8.0    2019-02-11 [2] CRAN (R 3.5.1)
#  lattice          0.20-38  2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval         0.2.1    2017-10-29 [2] CRAN (R 3.5.0)
#  magrittr         1.5      2014-11-22 [1] CRAN (R 3.5.0)
#  munsell          0.5.0    2018-06-12 [2] CRAN (R 3.5.1)
#  pillar           1.3.1    2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig        2.0.2    2018-08-16 [1] CRAN (R 3.5.1)
#  plyr             1.8.4    2016-06-08 [2] CRAN (R 3.5.0)
#  png              0.1-7    2013-12-03 [2] CRAN (R 3.5.0)
#  promises         1.0.1    2018-04-13 [2] CRAN (R 3.5.0)
#  purrr          * 0.3.1    2019-03-03 [2] CRAN (R 3.5.1)
#  R6               2.4.0    2019-02-14 [2] CRAN (R 3.5.1)
#  RColorBrewer   * 1.1-2    2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp             1.0.0    2018-11-07 [1] CRAN (R 3.5.1)
#  readxl         * 1.3.1    2019-03-13 [2] CRAN (R 3.5.1)
#  reshape2         1.4.3    2017-12-11 [2] CRAN (R 3.5.0)
#  rlang            0.3.1    2019-01-08 [1] CRAN (R 3.5.1)
#  rmote          * 0.3.4    2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors      * 0.20.1   2018-11-09 [1] Bioconductor
#  scales           1.0.0    2018-08-09 [2] CRAN (R 3.5.1)
#  servr            0.13     2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo    * 1.1.1    2018-11-05 [1] CRAN (R 3.5.1)
#  stringi          1.4.3    2019-03-12 [2] CRAN (R 3.5.1)
#  stringr          1.4.0    2019-02-10 [1] CRAN (R 3.5.1)
#  tibble         * 2.0.1    2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect       0.2.5    2018-10-11 [2] CRAN (R 3.5.1)
#  utf8             1.1.4    2018-05-24 [1] CRAN (R 3.5.0)
#  VennDiagram    * 1.6.20   2018-03-28 [1] CRAN (R 3.5.0)
#  withr            2.1.2    2018-03-15 [2] CRAN (R 3.5.0)
#  xfun             0.5      2019-02-20 [1] CRAN (R 3.5.1)
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
