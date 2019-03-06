## A cleaner and more focused script that explore_twas.R

library('tibble')
library('sessioninfo')
library('purrr')
# library('readr')
library('ggplot2')
# library('gplots')
# library('VennDiagram')
# library('RColorBrewer')
# library('dplyr')

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



## Subset by significant (TWAS FDR < 5%)
ttSig <- map(split(tt, tt$region), ~ .x[.x$TWAS.FDR < 0.05, ])
map_dfr(ttSig, dim)
# # A tibble: 2 x 2
#   DLPFC HIPPO
#   <int> <int>
# 1  3846  2952
# 2    44    44

ttSig_bonf <- map(split(tt, tt$region), ~ .x[.x$TWAS.Bonf < 0.05, ])
map_dfr(ttSig_bonf, dim)
# # A tibble: 2 x 2
#   DLPFC HIPPO
#   <int> <int>
# 1   448   345
# 2    44    44


map_int(ttSig, ~ length(unique(.x$geneid)))
# DLPFC HIPPO
#  1032   910

map_int(ttSig_bonf, ~ length(unique(.x$geneid)))
# DLPFC HIPPO
#   146   130

map_int(ttSig, ~ length(unique(.x$geneid[.x$TWAS.P < 5e-08])))
# DLPFC HIPPO
#    51    47

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
# [1] 0.4178535
with(tt, check_cor(TWAS.P, EQTL.P.computed))
# [1] 0.8177334

map_dbl(ttSig, ~ with(.x, check_cor(TWAS.P, BEST.GWAS.P.computed)))
#     DLPFC     HIPPO
# 0.5658657 0.5437862
map_dbl(ttSig, ~ with(.x, check_cor(TWAS.P, EQTL.P.computed)))
#     DLPFC     HIPPO
# 0.7324814 0.7362852

map_dbl(ttSig_bonf, ~ with(.x, check_cor(TWAS.P, BEST.GWAS.P.computed)))
#     DLPFC     HIPPO
# 0.5388277 0.4828275
map_dbl(ttSig_bonf, ~ with(.x, check_cor(TWAS.P, EQTL.P.computed)))
#     DLPFC     HIPPO
# 0.6384968 0.6563766

map_dbl(split(tt, tt$BEST.GWAS.status), ~ with(.x, check_cor(TWAS.P, BEST.GWAS.P.computed)))
#     Index     Other     Proxy
# 0.1046549 0.3369664 0.2392657
map_dbl(split(tt, tt$BEST.GWAS.status), ~ with(.x, check_cor(TWAS.P, EQTL.P.computed)))
#     Index     Other     Proxy
# 0.7951933 0.7947084 0.8124684


tt_sigonly <- tt[tt$TWAS.FDR < 0.05, ]
with(tt_sigonly, addmargins(table(BEST.GWAS.status, EQTL.status, useNA = 'ifany')))
#                 EQTL.status
# BEST.GWAS.status Index Other Proxy  Sum
#            Index     2   660   113  775
#            Other     0  7354    22 7376
#            Proxy     0  1345   345 1690
#            Sum       2  9359   480 9841

tt_sigonly_bonf <- tt[tt$TWAS.Bonf < 0.05, ]
with(tt_sigonly_bonf, addmargins(table(BEST.GWAS.status, EQTL.status, useNA = 'ifany')))
#                 EQTL.status
# BEST.GWAS.status Index Other Proxy  Sum
#            Index     2   252    60  314
#            Other     0   391    19  410
#            Proxy     0   356   266  622
#            Sum       2   999   345 1346

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
create_gwas_or_eqtl(tt_sigonly, 'pdf/twas_fdr5perc_vs_gwas_or_eqtl.pdf', 'FDR')
create_gwas_or_eqtl(tt_sigonly_bonf, 'pdf/twas_bonf5perc_vs_gwas_or_eqtl.pdf', 'Bonf')


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
create_by_status(tt_sigonly, 'pdf/twas_fdr5perc_by_status.pdf', 'FDR')
create_by_status(tt_sigonly_bonf, 'pdf/twas_bonf5perc_by_status.pdf', 'Bonf')


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
#               FALSE  2100   22 2122
#               TRUE    167  768  935
#               Sum    2267  790 3057
#
# $DLPFC$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   281    5 286
#               TRUE     20  100 120
#               Sum     301  105 406
#
# $DLPFC$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE  1106   19 1125
#               TRUE     61  366  427
#               Sum    1167  385 1552
#
# $DLPFC$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   536    7 543
#               TRUE     35  167 202
#               Sum     571  174 745
#
#
# $HIPPO
# $HIPPO$exon
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE  1354   24 1378
#               TRUE     95  482  577
#               Sum    1449  506 1955
#
# $HIPPO$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   203    3 206
#               TRUE     12   52  64
#               Sum     215   55 270
#
# $HIPPO$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE   895   13  908
#               TRUE     43  298  341
#               Sum     938  311 1249
#
# $HIPPO$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   437    7 444
#               TRUE     31  132 163
#               Sum     468  139 607

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
#               FALSE    38    0  38
#               TRUE     58  275 333
#               Sum      96  275 371
#
# $DLPFC$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    21    2  23
#               TRUE      8   50  58
#               Sum      29   52  81
#
# $DLPFC$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    27    3  30
#               TRUE     24  140 164
#               Sum      51  143 194
#
# $DLPFC$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    31    2  33
#               TRUE     14   75  89
#               Sum      45   77 122
#
#
# $HIPPO
# $HIPPO$exon
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    38    4  42
#               TRUE     35  169 204
#               Sum      73  173 246
#
# $HIPPO$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    17    1  18
#               TRUE     10   28  38
#               Sum      27   29  56
#
# $HIPPO$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    31    2  33
#               TRUE     16  110 126
#               Sum      47  112 159
#
# $HIPPO$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    27    1  28
#               TRUE     15   74  89
#               Sum      42   75 117


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
# 1    44    50 FALSE
# 2    59    66 TRUE
perc_locus(1.1)
#      HIPPO    DLPFC
# 1 57.28155 56.89655

## raggr eQTL locus that have a TWAS FDR <5% result
by_locus(0.05)
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    50    58 FALSE
# 2    53    58 TRUE
perc_locus(0.05)
#      HIPPO DLPFC
# 1 51.45631    50
perc_locus(0.05) / perc_locus(1.1) * 100
#      HIPPO    DLPFC
# 1 89.83051 87.87879

## raggr eQTL locus that have a TWAS Bonf <5% result
by_locus(0.05, var = 'TWAS.Bonf')
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    62    66 FALSE
# 2    41    50 TRUE
perc_locus(0.05, var = 'TWAS.Bonf')
#      HIPPO    DLPFC
# 1 39.80583 43.10345
perc_locus(0.05, var = 'TWAS.Bonf') / perc_locus(1.1) * 100
#      HIPPO    DLPFC
# 1 69.49153 75.75758


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
# 01  13     0     1
# 10   6     1     0
# 11  53     1     1

## Now with the ones that are TWAS FDR < 5%
shared_by_locus(0.05)
#    num HIPPO DLPFC
# 00   0     0     0
# 01  12     0     1
# 10   7     1     0
# 11  46     1     1

## Now with the ones that are TWAS Bonf < 5%
shared_by_locus(0.05, var = 'TWAS.Bonf')
#    num HIPPO DLPFC
# 00   0     0     0
# 01  13     0     1
# 10   4     1     0
# 11  37     1     1

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

pdf('pdf/venn_by_locus.pdf', useDingbats = FALSE)
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
# 3    21    18 FALSE exon
# 4    56    69 TRUE  exon
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

## Example: HIPPO gene
# 34 / (34 + 16) * 100
perc_feature(1.1)
#      HIPPO    DLPFC feature
# 1 68.00000 76.11940    gene
# 2 72.72727 79.31034    exon
# 3 75.00000 75.51020     jxn
# 4 70.76923 77.63158      tx


## Features by locus that have a TWAS FDR<5% result
clean_tabs(by_feature(cut = 0.05))
# # A tibble: 8 x 4
#   HIPPO DLPFC state feature
#   <int> <int> <lgl> <chr>
# 1    25    28 FALSE gene
# 2    25    39 TRUE  gene
# 3    29    28 FALSE exon
# 4    48    59 TRUE  exon
# 5    34    35 FALSE jxn
# 6    54    63 TRUE  jxn
# 7    31    25 FALSE tx
# 8    34    51 TRUE  tx

perc_feature(0.05)
#      HIPPO    DLPFC feature
# 1 50.00000 58.20896    gene
# 2 62.33766 67.81609    exon
# 3 61.36364 64.28571     jxn
# 4 52.30769 67.10526      tx

## Compute the percent using as denominator the number of features considered
cbind(perc_feature(0.05)[, 1:2] / perc_feature(1.1)[, 1:2] * 100, feature = features)
#      HIPPO    DLPFC feature
# 1 73.52941 76.47059    gene
# 2 85.71429 85.50725    exon
# 3 81.81818 85.13514     jxn
# 4 73.91304 86.44068      tx

## Now with Bonf < 5%
clean_tabs(by_feature(cut = 0.05, var = 'TWAS.Bonf'))
# # A tibble: 8 x 4
#   HIPPO DLPFC state feature
#   <int> <int> <lgl> <chr>
# 1    32    38 FALSE gene
# 2    18    29 TRUE  gene
# 3    46    44 FALSE exon
# 4    31    43 TRUE  exon
# 5    51    54 FALSE jxn
# 6    37    44 TRUE  jxn
# 7    35    44 FALSE tx
# 8    30    32 TRUE  tx

perc_feature(0.05, var = 'TWAS.Bonf')
#      HIPPO    DLPFC feature
# 1 36.00000 43.28358    gene
# 2 40.25974 49.42529    exon
# 3 42.04545 44.89796     jxn
# 4 46.15385 42.10526      tx

cbind(perc_feature(0.05, var = 'TWAS.Bonf')[, 1:2] / perc_feature(1.1)[, 1:2] * 100, feature = features)
#      HIPPO    DLPFC feature
# 1 52.94118 56.86275    gene
# 2 55.35714 62.31884    exon
# 3 56.06061 59.45946     jxn
# 4 65.21739 54.23729      tx

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
# 01  20     0     1
# 10   7     1     0
# 11  49     1     1
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
# 01  17     0     1
# 10   3     1     0
# 11  22     1     1
#
# $exon
#    num HIPPO DLPFC
# 00   0     0     0
# 01  17     0     1
# 10   6     1     0
# 11  42     1     1
#
# $jxn
#    num HIPPO DLPFC
# 00   0     0     0
# 01  17     0     1
# 10   8     1     0
# 11  46     1     1
#
# $tx
#    num HIPPO DLPFC
# 00   0     0     0
# 01  21     0     1
# 10   4     1     0
# 11  30     1     1

## TWAS Bonf <5%
shared_by_locus_by_feature(0.05, var = 'TWAS.Bonf')
# $gene
#    num HIPPO DLPFC
# 00   0     0     0
# 01  15     0     1
# 10   4     1     0
# 11  14     1     1
#
# $exon
#    num HIPPO DLPFC
# 00   0     0     0
# 01  17     0     1
# 10   5     1     0
# 11  26     1     1
#
# $jxn
#    num HIPPO DLPFC
# 00   0     0     0
# 01  15     0     1
# 10   8     1     0
# 11  29     1     1
#
# $tx
#    num HIPPO DLPFC
# 00   0     0     0
# 01  11     0     1
# 10   9     1     0
# 11  21     1     1


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

pdf('pdf/venn_by_locus_by_feature.pdf', useDingbats = FALSE)
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
# 1  5482  3977     gene
# 2 39131 25686     exon
# 3 20653 16107      jxn
# 4  9061  7154       tx

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
# 1 23.425348 16.994274     gene
# 2 10.265323  6.738266     exon
# 3  7.670253  5.981928      jxn
# 4 10.184446  8.041003       tx

## Number of features with TWAS FDR<5%
n_feat_twas5perc <- cbind(map_dfr(ttSig, ~ 
    map_dbl(
        split(.x, factor(.x$feature, levels = features)),
        ~ nrow(.x)
    )
), features)
n_feat_twas5perc
#   DLPFC HIPPO features
# 1   406   270     gene
# 2  3057  1955     exon
# 3  1552  1249      jxn
# 4   745   607       tx

cbind(n_feat_twas5perc[, 1:2] / n_feat_considered[, 1:2] * 100, features)
#      DLPFC    HIPPO features
# 1 7.406056 6.789037     gene
# 2 7.812220 7.611150     exon
# 3 7.514647 7.754393      jxn
# 4 8.222051 8.484764       tx

## Now with TWAS Bonf <5%
n_feat_twas5perc_bonf <- cbind(map_dfr(ttSig_bonf, ~ 
    map_dbl(
        split(.x, factor(.x$feature, levels = features)),
        ~ nrow(.x)
    )
), features)
n_feat_twas5perc_bonf
#   DLPFC HIPPO features
# 1    81    56     gene
# 2   371   246     exon
# 3   194   159      jxn
# 4   122   117       tx

cbind(n_feat_twas5perc_bonf[, 1:2] / n_feat_considered[, 1:2] * 100, features)
#       DLPFC     HIPPO features
# 1 1.4775629 1.4080966     gene
# 2 0.9480974 0.9577202     exon
# 3 0.9393308 0.9871484      jxn
# 4 1.3464298 1.6354487       tx


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
# 01 1185     0     1
# 10 2690     1     0
# 11 2792     1     1
#
# $exon
#      num DLPFC HIPPO
# 00     0     0     0
# 01 11906     0     1
# 10 25351     1     0
# 11 13780     1     1
#
# $jxn
#      num DLPFC HIPPO
# 00     0     0     0
# 01  7669     0     1
# 10 12215     1     0
# 11  8438     1     1
#
# $tx
#     num DLPFC HIPPO
# 00    0     0     0
# 01 3040     0     1
# 10 4947     1     0
# 11 4114     1     1

## TWAS FDR <5%
shared_by_feature(0.05)
# $gene
#    num DLPFC HIPPO
# 00   0     0     0
# 01 132     0     1
# 10 268     1     0
# 11 138     1     1
#
# $exon
#     num DLPFC HIPPO
# 00    0     0     0
# 01 1201     0     1
# 10 2303     1     0
# 11  754     1     1
#
# $jxn
#     num DLPFC HIPPO
# 00    0     0     0
# 01  745     0     1
# 10 1048     1     0
# 11  504     1     1
#
# $tx
#    num DLPFC HIPPO
# 00   0     0     0
# 01 347     0     1
# 10 485     1     0
# 11 260     1     1

## TWAS Bonf <5%
shared_by_feature(0.05, var = 'TWAS.Bonf')
# $gene
#    num DLPFC HIPPO
# 00   0     0     0
# 01  29     0     1
# 10  54     1     0
# 11  27     1     1
#
# $exon
#    num DLPFC HIPPO
# 00   0     0     0
# 01 159     0     1
# 10 284     1     0
# 11  87     1     1
#
# $jxn
#    num DLPFC HIPPO
# 00   0     0     0
# 01 108     0     1
# 10 143     1     0
# 11  51     1     1
#
# $tx
#    num DLPFC HIPPO
# 00   0     0     0
# 01  76     0     1
# 10  81     1     0
# 11  41     1     1

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

pdf('pdf/venn_by_feature.pdf', useDingbats = FALSE)
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
# 01 1185     0     1
# 10 2690     1     0
# 11 2792     1     1
#
# $exon
#      num DLPFC HIPPO
# 00     0     0     0
# 01  2660     0     1
# 10 10948     1     0
# 11 28183     1     1
#
# $jxn
#      num DLPFC HIPPO
# 00     0     0     0
# 01  1827     0     1
# 10  3367     1     0
# 11 13177     1     1
#
# $tx
#     num DLPFC HIPPO
# 00    0     0     0
# 01 1662     0     1
# 10 3076     1     0
# 11 5985     1     1

## TWAS FDR <5%
shared_by_geneid(0.05)
# $gene
#    num DLPFC HIPPO
# 00   0     0     0
# 01 132     0     1
# 10 268     1     0
# 11 138     1     1
#
# $exon
#     num DLPFC HIPPO
# 00    0     0     0
# 01  497     0     1
# 10 1317     1     0
# 11 1740     1     1
#
# $jxn
#    num DLPFC HIPPO
# 00   0     0     0
# 01 296     0     1
# 10 463     1     0
# 11 814     1     1
#
# $tx
#    num DLPFC HIPPO
# 00   0     0     0
# 01 248     0     1
# 10 368     1     0
# 11 377     1     1

## TWAS Bonf <5%
shared_by_geneid(0.05, var = 'TWAS.Bonf')
# $gene
#    num DLPFC HIPPO
# 00   0     0     0
# 01  29     0     1
# 10  54     1     0
# 11  27     1     1
#
# $exon
#    num DLPFC HIPPO
# 00   0     0     0
# 01  61     0     1
# 10 138     1     0
# 11 233     1     1
#
# $jxn
#    num DLPFC HIPPO
# 00   0     0     0
# 01  45     0     1
# 10  66     1     0
# 11  97     1     1
#
# $tx
#    num DLPFC HIPPO
# 00   0     0     0
# 01  60     0     1
# 10  65     1     0
# 11  57     1     1

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

pdf('pdf/venn_by_feature_using_geneid.pdf', useDingbats = FALSE)
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

pdf('pdf/venn_by_feature_using_geneid_across_features.pdf', useDingbats = FALSE)
make_pretty_venn_shared_by_geneid2(1.1, 'Features with TWAS weights (by gene ID)')
make_pretty_venn_shared_by_geneid2(0.05, 'Features with TWAS FDR<5% (by gene ID)')
make_pretty_venn_shared_by_geneid2(0.05, 'Features with TWAS Bonf<5% (by gene ID)', var = 'TWAS.Bonf')
dev.off()

system('rm VennDiagram*')

cbind(map_dfr(split(tt, tt$region), ~ map_dfr(split(.x, .x$feature), ~ sum(.x$TWAS.FDR < 0.05))), region = c('DLPFC', 'HIPPO'))
#   exon gene  jxn  tx region
# 1 3057  406 1552 745  DLPFC
# 2 1955  270 1249 607  HIPPO

cbind(map_dfr(split(tt, tt$region), ~ map_dfr(split(.x, .x$feature), ~ sum(.x$TWAS.P < 5e-08))), region = c('DLPFC', 'HIPPO'))
#   exon gene jxn tx region
# 1  205   24  90 47  DLPFC
# 2  128   21  53 37  HIPPO

cbind(map_dfr(split(tt, tt$region), ~ map_dfr(split(.x, .x$feature), ~ sum(p.adjust(.x$TWAS.P, 'bonf') < 0.05))), region = c('DLPFC', 'HIPPO'))
#   exon gene jxn  tx region
# 1  371   81 194 122  DLPFC
# 2  246   56 159 117  HIPPO



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

