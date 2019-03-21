## A cleaner and more focused script that explore_twas.R

library('tibble')
library('sessioninfo')
library('purrr')
library('readr')
library('ggplot2')
library('gplots')
library('VennDiagram')
library('RColorBrewer')
source('twas_functions.R')

load('rda/twas_exp.Rdata', verbose = TRUE)

## Andrew's exploration code that focuses on the 'all' part
tt <- twas_exp$all
## Drop TWAS NA p-values
tt <- tt[!is.na(tt$TWAS.P), ]
## Focus on CLOZUK+PGC2 (psycm) GWAS
tt <- tt[which(tt$type == "psycm"),]

## Add GWAS p-value and OR from the original sumstats file
original <- read_tsv('psycm/clozuk_pgc2.meta.sumstats.txt')
original$CHR[original$CHR == 23] <- 'X'
original$hg19_pos <- with(original, paste0(CHR, ':', BP))

snpmap <- read_tsv('psycm/clozuk_pgc2.meta.reformatted.sumstats_hg38_ourname',
    col_types = cols(
      SNP = col_character(),
      A1 = col_character(),
      A2 = col_character(),
      Z = col_double(),
      N = col_double(),
      chr = col_character(),
      basepair = col_double(),
      basepairhg19 = col_double(),
      originalSNP = col_character()
    )
)

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




m_to_map <- match(tt$BEST.GWAS.ID, snpmap$SNP)
table(is.na(tt$BEST.GWAS.ID))
#  FALSE
# 127251
table(is.na(m_to_map))
## Hm... I'm not sure why some are NAs
#  FALSE   TRUE
# 126157   1094
print(tt[head(which(is.na(m_to_map))), ], width = 200)

## Hm....
m_to_map_qtl <- match(tt$EQTL.ID, snpmap$SNP)
table(is.na(m_to_map_qtl))
#  FALSE   TRUE
# 125935   1316

addmargins(table(
    'By BEST.GWAS.ID' = is.na(m_to_map),
    'By EQTL.ID' = is.na(m_to_map_qtl)
))
#                By EQTL.ID
# By BEST.GWAS.ID  FALSE   TRUE    Sum
#           FALSE 124876   1281 126157
#           TRUE    1059     35   1094
#           Sum   125935   1316 127251

## Well, after that it all looks ok
m_to_ori <- match(snpmap$originalSNP[m_to_map], original$SNP)
table(is.na(m_to_ori))
#  FALSE   TRUE
# 126157   1094

## Hm... it's odd that the same number don't match by either name or chr position
m_to_ori2 <- match(tt$BEST.GWAS.pos_hg19, original$hg19_pos)
table(is.na(m_to_ori), is.na(m_to_ori2))
#        FALSE   TRUE
# FALSE 125076   1081
# TRUE    1081     13
## Supplement m_to_ori (by name) with the matching by chr and position in hg19
m_to_ori[is.na(m_to_ori)] <- m_to_ori2[is.na(m_to_ori)]
table(is.na(m_to_ori))
#  FALSE   TRUE
# 127238     13

m_to_map_qtl2 <- match(tt$EQTL.pos_hg19, original$hg19_pos)
table(is.na(m_to_map_qtl), is.na(m_to_map_qtl2))
#        FALSE   TRUE
# FALSE 125935      0
# TRUE       0   1316



## Lets get the originally reported summarized p-values
BEST.GWAS.P <- original$P[m_to_ori]
## This is how the calculate the p-values displayed by FUSION-TWAS
# https://github.com/gusevlab/fusion_twas/blob/master/FUSION.post_process.R#L641
BEST.GWAS.P.computed <- 2*(pnorm( abs(tt$BEST.GWAS.Z ) , lower.tail=F ))


x <- -log10(BEST.GWAS.P[!is.na(m_to_ori)])
y <- -log10(BEST.GWAS.P.computed[!is.na(m_to_ori)])
## The difference seems small, likely due to the number of decimals in the original table
summary(x - y)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.0183092 -0.0033636  0.0002139  0.0002255  0.0040540  0.0166786
summary(abs(x - y))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 4.611e-06 1.705e-03 3.736e-03 4.112e-03 6.098e-03 1.831e-02
head(original$P)
# [1] 8.120e-01 2.585e-01 9.777e-01 7.701e-01 5.431e-01 1.149e-05
head(original$P) * 1e6
# [1] 812000.00 258500.00 977700.00 770100.00 543100.00     11.49

## Ok, assign them to our table
tt$BEST.GWAS.P <- original$P[m_to_ori]
tt$BEST.GWAS.OR <- original$OR[m_to_ori]
tt$BEST.GWAS.SE <- original$SE[m_to_ori]
tt$BEST.GWAS.P.computed <- 2*(pnorm( abs(tt$BEST.GWAS.Z ) , lower.tail=F ))
tt$EQTL.P.computed <- 2*(pnorm( abs(tt$EQTL.GWAS.Z ) , lower.tail=F ))

## Compute FDR/Bonf by region for each feature 4 features
tt <- map_dfr(split(tt, tt$region), function(reg) {
    res <- map_dfr(split(reg, reg$feature), function(reg_feat) {
        reg_feat$TWAS.FDR <- p.adjust(reg_feat$TWAS.P, 'fdr')
        reg_feat$TWAS.Bonf <- p.adjust(reg_feat$TWAS.P, 'bonf')
        reg_feat$BEST.GWAS.FDR <- p.adjust(reg_feat$BEST.GWAS.P, 'fdr')
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


addmargins(table(
    'proxy' = snpMap$pos_hg19 %in% riskLoci$hg19POS2,
    'index' = snpMap$pos_hg19 %in% indexLoci$hg19POS
))
#        index
# proxy     FALSE    TRUE     Sum
#   FALSE 7014124       0 7014124
#   TRUE     9600     135    9735
#   Sum   7023724     135 7023859


snpMap$Status <- 'Other'
snpMap$Status[snpMap$pos_hg19 %in% riskLoci$hg19POS2] <- 'Proxy'
snpMap$Status[snpMap$pos_hg19 %in% indexLoci$hg19POS] <- 'Index'
table(snpMap$Status)
# Index   Other   Proxy
#   135 7014124    9600

## One of the index SNPs has 2 names
length(unique(riskLoci$SNP1_Name))
# [1] 180
which(table(gsub(':.*', '', unique(riskLoci$SNP1_Name))) > 1)
# rs1023497
#         5
## Turns out that it's multi-allelic
unique(riskLoci$SNP1_Name)[grep('rs1023497', unique(riskLoci$SNP1_Name))]
# [1] "rs1023497:42340508:C:G" "rs1023497:42340508:C:A"

tt <- as_tibble(cbind(tt,
    get_proxy_info(tt$BEST.GWAS.pos_hg19, 'BEST.GWAS.'),
    get_proxy_info(tt$EQTL.pos_hg19, 'EQTL.')
))
## BEST.GWAS info
# status
#  Index  Other  Proxy
#   2189 119830   5232
# [1] 113
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -331128  -10256       0   -4290    4113  322056
## EQTL info
# status
#  Index  Other  Proxy
#      3 126738    510
# [1] 50
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -339607  -73401  -16684  -17830   34305  437328

## Label the 'other' with genome significant p-values as proxy
table(tt$BEST.GWAS.P.computed < 5e-8, tt$BEST.GWAS.status)
#        Index  Other  Proxy
# FALSE      0 116979    512
# TRUE    2189   1757   4720
tt$BEST.GWAS.status[tt$BEST.GWAS.P.computed < 5e-8 & tt$BEST.GWAS.status == 'Other'] <- 'Proxy'
table(tt$BEST.GWAS.P.computed < 5e-8, tt$BEST.GWAS.status)
#        Index  Other  Proxy
# FALSE      0 118073    512
# TRUE    2189      0   6477

print(tt, width = 600)
head(as.data.frame(tt))


table(table(tt$ID))
#     1     2
# 69003 29124

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

region_twas_z <- get_variable_by_region('TWAS.Z', NAs_0 = TRUE)

## How I figured out something weird :P
## Some values on the cross with low/high Zs were labeled as "None":
## they were NAs on one of the two regions
head(subset(region_twas_z, FDR.5perc == 'None' & DLPFC < -5))


table(region_twas_z$in_both)
#
# FALSE  TRUE
# 69003 29124
table(region_twas_z$in_both) / nrow(region_twas_z) * 100
#   FALSE    TRUE
# 70.3201 29.6799
addmargins(table(
    'in both' = region_twas_z$in_both,
    'FDR < 0.05' = region_twas_z$FDR.5perc,
    useNA = 'ifany'
))
# in both  None DLPFC HIPPO  Both   Sum
#   FALSE 69003     0     0     0 69003
#   TRUE  26148   707   613  1656 29124
#   Sum   95151   707   613  1656 98127

addmargins(table(
    'in both' = region_twas_z$in_both,
    'Bonf < 0.05' = region_twas_z$Bonf.5perc,
    useNA = 'ifany'
))
# in both  None DLPFC HIPPO  Both   Sum
#   FALSE 68314   436   253     0 69003
#   TRUE  28673   126   119   206 29124
#   Sum   96987   562   372   206 98127

table(region_twas_z$BEST.GWAS.status)
## Initial version
# Other DLPFC HIPPO  Both
# 92402  2675  1354  1696
## Latest
# Other Risk Locus
# 91007       7120



pdf('pdf/twas_z.pdf', useDingbats = FALSE, width = 24, height = 14)
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

## Find the discordant ones
options(width = 200)
subset(region_twas_z, FDR.5perc == 'Both' & sign(DLPFC) != sign(HIPPO))
#                     ID feature             geneid genesymbol    DLPFC    HIPPO in_both TWAS.FDR_DLPFC TWAS.FDR_HIPPO BEST.GWAS.status_DLPFC BEST.GWAS.status_HIPPO FDR.5perc BEST.GWAS.status
# 5269 ENST00000422145.7      tx  ENSG00000236922.9  LINC01378 -2.96921 3.105439    TRUE     0.03937847     0.02796831                  Other                  Other      Both            Other
# 5375           e938460    exon ENSG00000006042.11     TMEM98 -2.94795 3.409439    TRUE     0.04424707     0.01432869                  Other                  Other      Both            Other

# https://www.genecards.org/cgi-bin/carddisp.pl?gene=LINC01378&keywords=ENSG00000236922
# https://www.genecards.org/cgi-bin/carddisp.pl?gene=TMEM98&keywords=ENSG00000006042

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
#   None   3416 2413 5829
#   DLPFC   126   61  187
#   HIPPO    58   51  109
#   Both      0   94   94
#   Sum    3600 2619 6219
#
# $gene$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE Sum
#   None    187  113 300
#   DLPFC    70   11  81
#   HIPPO    18    5  23
#   Both      0   44  44
#   Sum     275  173 448
#
#
# $exon
# $exon$Other
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None  32437 11767 44204
#   DLPFC  1312   288  1600
#   HIPPO   622   232   854
#   Both      0   500   500
#   Sum   34371 12787 47158
#
# $exon$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE  Sum
#   None   1964  611 2575
#   DLPFC   633   70  703
#   HIPPO   289   58  347
#   Both      0  254  254
#   Sum    2886  993 3879
#
#
# $jxn
# $jxn$Other
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None  17539  7201 24740
#   DLPFC   605   155   760
#   HIPPO   416   133   549
#   Both      0   346   346
#   Sum   18560  7835 26395
#
# $jxn$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE  Sum
#   None    908  377 1285
#   DLPFC   259   29  288
#   HIPPO   157   39  196
#   Both      0  158  158
#   Sum    1324  603 1927
#
#
# $tx
# $tx$Other
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None   6964  3481 10445
#   DLPFC   277    76   353
#   HIPPO   169    85   254
#   Both      0   183   183
#   Sum    7410  3825 11235
#
# $tx$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE Sum
#   None    379  185 564
#   DLPFC   115   17 132
#   HIPPO    83   10  93
#   Both      0   77  77
#   Sum     577  289 866



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
#    None   3589 2599 6188
#    DLPFC     8    6   14
#    HIPPO     3    7   10
#    Both      0    7    7
#    Sum    3600 2619 6219
#
# $gene$`Risk Locus`
#         In both
# Bonf <5% FALSE TRUE Sum
#    None    230  139 369
#    DLPFC    31    9  40
#    HIPPO    14    5  19
#    Both      0   20  20
#    Sum     275  173 448
#
#
# $exon
# $exon$Other
#         In both
# Bonf <5% FALSE  TRUE   Sum
#    None  34337 12749 47086
#    DLPFC    18    16    34
#    HIPPO    16    18    34
#    Both      0     4     4
#    Sum   34371 12787 47158
#
# $exon$`Risk Locus`
#         In both
# Bonf <5% FALSE TRUE  Sum
#    None   2590  831 3421
#    DLPFC   203   47  250
#    HIPPO    93   32  125
#    Both      0   83   83
#    Sum    2886  993 3879
#
#
# $jxn
# $jxn$Other
#         In both
# Bonf <5% FALSE  TRUE   Sum
#    None  18534  7806 26340
#    DLPFC    13    11    24
#    HIPPO    13    15    28
#    Both      0     3     3
#    Sum   18560  7835 26395
#
# $jxn$`Risk Locus`
#         In both
# Bonf <5% FALSE TRUE  Sum
#    None   1165  515 1680
#    DLPFC    99   20  119
#    HIPPO    60   20   80
#    Both      0   48   48
#    Sum    1324  603 1927
#
#
# $tx
# $tx$Other
#         In both
# Bonf <5% FALSE  TRUE   Sum
#    None   7378  3806 11184
#    DLPFC    18     6    24
#    HIPPO    14     6    20
#    Both      0     7     7
#    Sum    7410  3825 11235
#
# $tx$`Risk Locus`
#         In both
# Bonf <5% FALSE TRUE Sum
#    None    491  228 719
#    DLPFC    46   11  57
#    HIPPO    40   16  56
#    Both      0   34  34
#    Sum     577  289 866


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
# sapply(tt[, 49:51], function(x) { sum(is.na(x)) })


## Subset by significant (TWAS FDR < 5%)
ttSig <- map(split(tt, tt$region), ~ .x[.x$TWAS.FDR < 0.05, ])
map_dfr(ttSig, dim)
# # A tibble: 2 x 2
#   DLPFC HIPPO
#   <int> <int>
# 1  5760  4081
# 2    51    51

ttSig_bonf <- map(split(tt, tt$region), ~ .x[.x$TWAS.Bonf < 0.05, ])
map_dfr(ttSig_bonf, dim)
# # A tibble: 2 x 2
#   DLPFC HIPPO
#   <int> <int>
# 1   768   578
# 2    51    51


map_int(ttSig, ~ length(unique(.x$geneid)))
# DLPFC HIPPO
#  1514  1255

## Original numbers when computing FDR across all 4 features at the same time:
# DLPFC HIPPO
#  1519  1256

map_int(ttSig_bonf, ~ length(unique(.x$geneid)))
# DLPFC HIPPO
#   240   218

map_int(ttSig, ~ length(unique(.x$geneid[.x$TWAS.P < 5e-08])))
# DLPFC HIPPO
#   115    84

dim(tt)
# [1] 127251     51
dim(region_twas_z)
# [1] 98127    16



## save for later
save(tt, ttSig, ttSig_bonf, get_variable_by_region, ttReg_map, region_twas_z, file = 'rda/tt_objects.Rdata')




## Continue
load('rda/tt_objects.Rdata', verbose = TRUE)


## Check correlations among FDR corrected p-values between TWAS
## and either BEST GWAS FDR or EQTL GWAS FDR
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
# 0.1046549 0.2937820 0.2775410
map_dbl(split(tt, tt$BEST.GWAS.status), ~ with(.x, check_cor(TWAS.P, EQTL.P.computed)))
#     Index     Other     Proxy
# 0.7951933 0.7767584 0.8272976


tt_sigonly <- tt[tt$TWAS.FDR < 0.05, ]
with(tt_sigonly, addmargins(table(BEST.GWAS.status, EQTL.status, useNA = 'ifany')))
#                 EQTL.status
# BEST.GWAS.status Index Other Proxy  Sum
#            Index     2   660   113  775
#            Other     0  6912     0 6912
#            Proxy     0  1787   367 2154
#            Sum       2  9359   480 9841

tt_sigonly_bonf <- tt[tt$TWAS.Bonf < 0.05, ]
with(tt_sigonly_bonf, addmargins(table(BEST.GWAS.status, EQTL.status, useNA = 'ifany')))
#                 EQTL.status
# BEST.GWAS.status Index Other Proxy  Sum
#            Index     2   252    60  314
#            Other     0   230     0  230
#            Proxy     0   517   285  802
#            Sum       2   999   345 1346

create_gwas_or_eqtl(tt_sigonly, 'pdf/twas_fdr5perc_vs_gwas_or_eqtl.pdf', 'FDR')
create_gwas_or_eqtl(tt_sigonly_bonf, 'pdf/twas_bonf5perc_vs_gwas_or_eqtl.pdf', 'Bonf')

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
#               TRUE      0  935  935
#               Sum    2100  957 3057
#
# $DLPFC$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   281    5 286
#               TRUE      0  120 120
#               Sum     281  125 406
#
# $DLPFC$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE  1106   19 1125
#               TRUE      0  427  427
#               Sum    1106  446 1552
#
# $DLPFC$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   536    7 543
#               TRUE      0  202 202
#               Sum     536  209 745
#
#
# $HIPPO
# $HIPPO$exon
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE  1354   24 1378
#               TRUE      0  577  577
#               Sum    1354  601 1955
#
# $HIPPO$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   203    3 206
#               TRUE      0   64  64
#               Sum     203   67 270
#
# $HIPPO$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE  Sum
#               FALSE   895   13  908
#               TRUE      0  341  341
#               Sum     895  354 1249
#
# $HIPPO$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE   437    7 444
#               TRUE      0  163 163
#               Sum     437  170 607

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
#               TRUE      0  333 333
#               Sum      38  333 371
#
# $DLPFC$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    21    2  23
#               TRUE      0   58  58
#               Sum      21   60  81
#
# $DLPFC$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    27    3  30
#               TRUE      0  164 164
#               Sum      27  167 194
#
# $DLPFC$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    31    2  33
#               TRUE      0   89  89
#               Sum      31   91 122
#
#
# $HIPPO
# $HIPPO$exon
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    38    4  42
#               TRUE      0  204 204
#               Sum      38  208 246
#
# $HIPPO$gene
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    17    1  18
#               TRUE      0   38  38
#               Sum      17   39  56
#
# $HIPPO$jxn
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    31    2  33
#               TRUE      0  126 126
#               Sum      31  128 159
#
# $HIPPO$tx
#                    Risk Locus (by BEST GWAS)
# BEST GWAS P < 5e-08 FALSE TRUE Sum
#               FALSE    27    1  28
#               TRUE      0   89  89
#               Sum      27   90 117


## Add locus considered section

## Read in the files that Emily cleaned up at
## https://github.com/LieberInstitute/brainseq_phase2/blob/master/eQTL_GWAS_riskSNPs/create_eqtl_table_indexInfo.R
raggr_clean_files <- c(
    'HIPPO' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/raggr_179_snps_hippo_eqtls_fdr01.csv',
    'DLPFC' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/raggr_179_snps_dlp_eqtls_fdr01.csv'
)
raggr_clean <- map(raggr_clean_files, read.csv, stringsAsFactors = FALSE)
names(raggr_clean) <- names(raggr_clean_files)


clean_by_state(
    map2_dfr(
        raggr_clean,
        split(tt, factor(tt$region, levels = names(raggr_clean))),
        ~ table(unique(.x$IndexSNP) %in% .y$BEST.GWAS.indexSNP)
    )
)
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    26    30 FALSE
# 2    77    86 TRUE

stopifnot(identical(by_locus(1.1), by_locus(1.1, var = 'TWAS.Bonf')))

## raggr eQTL locus considered in the TWAS analysis
by_locus(1.1)
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    26    30 FALSE
# 2    77    86 TRUE
perc_locus(1.1)
#      HIPPO    DLPFC
# 1 74.75728 74.13793

## raggr eQTL locus that have a TWAS FDR <5% result
by_locus(0.05)
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    34    43 FALSE
# 2    69    73 TRUE
perc_locus(0.05)
#      HIPPO    DLPFC
# 1 66.99029 62.93103
perc_locus(0.05) / perc_locus(1.1) * 100
#      HIPPO    DLPFC
# 1 89.61039 84.88372

## raggr eQTL locus that have a TWAS Bonf <5% result
by_locus(0.05, var = 'TWAS.Bonf')
# # A tibble: 2 x 3
#   HIPPO DLPFC state
#   <int> <int> <lgl>
# 1    54    62 FALSE
# 2    49    54 TRUE
perc_locus(0.05, var = 'TWAS.Bonf')
#      HIPPO    DLPFC
# 1 47.57282 46.55172
perc_locus(0.05, var = 'TWAS.Bonf') / perc_locus(1.1) * 100
#      HIPPO   DLPFC
# 1 63.63636 62.7907


## Overlap by region of all loci considered
shared_by_locus(1.1)
#    num HIPPO DLPFC
# 00   0     0     0
# 01  15     0     1
# 10   6     1     0
# 11  71     1     1

## Now with the ones that are TWAS FDR < 5%
shared_by_locus(0.05)
#    num HIPPO DLPFC
# 00   0     0     0
# 01  13     0     1
# 10   9     1     0
# 11  60     1     1

## Now with the ones that are TWAS Bonf < 5%
shared_by_locus(0.05, var = 'TWAS.Bonf')
#    num HIPPO DLPFC
# 00   0     0     0
# 01  12     0     1
# 10   7     1     0
# 11  42     1     1


pdf('pdf/venn_by_locus.pdf', useDingbats = FALSE)
make_pretty_venn(1.1, 'rAggr loci considered in TWAS')
make_pretty_venn(0.05, 'rAggr loci with TWAS FDR<5%')
make_pretty_venn(0.05, 'rAggr loci with TWAS Bonf<5%', var = 'TWAS.Bonf')
dev.off()
system('rm VennDiagram*')

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



pdf('pdf/venn_by_feature.pdf', useDingbats = FALSE)
make_pretty_venn_shared_by_feature(1.1, 'Features with TWAS weights')
make_pretty_venn_shared_by_feature(0.05, 'Features with TWAS FDR<5%')
make_pretty_venn_shared_by_feature(0.05, 'Features with TWAS Bonf<5%', var = 'TWAS.Bonf')
dev.off()
system('rm VennDiagram*')


## Now by gene id
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

pdf('pdf/venn_by_feature_using_geneid.pdf', useDingbats = FALSE)
make_pretty_venn_shared_by_geneid(1.1, 'Features with TWAS weights (by gene ID)')
make_pretty_venn_shared_by_geneid(0.05, 'Features with TWAS FDR<5% (by gene ID)')
make_pretty_venn_shared_by_geneid(0.05, 'Features with TWAS Bonf<5% (by gene ID)', var = 'TWAS.Bonf')
dev.off()
system('rm VennDiagram*')


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



## Compare SCZD t vs TWAS z
pdf('pdf/sczd_t_vs_twas_z.pdf', useDingbats = FALSE, width = 21, height = 18)

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
#  hms              0.4.2    2018-03-10 [2] CRAN (R 3.5.0)
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
#  readr          * 1.3.1    2018-12-21 [1] CRAN (R 3.5.1)
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
