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


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

