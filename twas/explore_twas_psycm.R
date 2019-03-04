## A cleaner and more focused script that explore_twas.R

library('tibble')
library('sessioninfo')
library('purrr')
library('readr')
library('ggplot2')
# library('dplyr')

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

## Compute FDR by region for each feature 4 features
tt <- map_dfr(split(tt, tt$region), function(reg) {
    res <- map_dfr(split(reg, reg$feature), function(reg_feat) {
        reg_feat$TWAS.FDR <- p.adjust(reg_feat$TWAS.P, 'fdr')
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
length(unique(riskLoci$SNP2_Name))
# [1] 10981
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
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
# hmm <- tt
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
        BEST.GWAS.status_DLPFC = tt$BEST.GWAS.status[ttReg_map$i_DLPFC],
        BEST.GWAS.status_HIPPO = tt$BEST.GWAS.status[ttReg_map$i_HIPPO],
        stringsAsFactors = FALSE
    )

    # print(with(res, addmargins(table(
 #        'DLPFC' = res$TWAS.FDR_DLPFC < 0.05,
 #        'HIPPO' = res$TWAS.FDR_HIPPO < 0.05,
 #        useNA = 'ifany'
 #    ))))
#        HIPPO
# DLPFC   FALSE  TRUE  <NA>   Sum
#   FALSE 26148   613 41806 68567
#   TRUE    707  1656  3397  5760
#   <NA>  21988  1812     0 23800
#   Sum   48843  4081 45203 98127
    
    res$FDR.5perc <- 'None'
    res$FDR.5perc[res$TWAS.FDR_DLPFC < 0.05] <- 'DLPFC'
    res$FDR.5perc[res$TWAS.FDR_HIPPO < 0.05] <- 'HIPPO'
    res$FDR.5perc[res$TWAS.FDR_DLPFC < 0.05 & res$TWAS.FDR_HIPPO < 0.05] <- 'Both'
    
    # print(with(res, table(BEST.GWAS.status_DLPFC, BEST.GWAS.status_HIPPO, useNA = 'ifany')))
#                       BEST.GWAS.status_HIPPO
# BEST.GWAS.status_DLPFC Index Other Proxy  <NA>
#                  Index   524     0     0   723
#                  Other     0 27428     0 42528
#                  Proxy     0     0  1172  1952
#                  <NA>    418 22446   936     0
#
    res$BEST.GWAS.status <- 'Other'
    res$BEST.GWAS.status[c(which(res$BEST.GWAS.status_DLPFC != 'Other'),  which(res$BEST.GWAS.status_HIPPO != 'Other'))] <- 'Risk Locus'
    
    if(NAs_0 == TRUE) {
        res$DLPFC[is.na(res$DLPFC)] <- 0
        res$HIPPO[is.na(res$HIPPO)] <- 0
    }
    
    ## Make the features as factor, so its looks ok when plotting
    res$feature <- factor(res$feature, levels = c('gene', 'exon', 'jxn', 'tx'))
    res$FDR.5perc <- factor(res$FDR.5perc, levels = c('None', 'DLPFC', 'HIPPO', 'Both'))
    res$BEST.GWAS.status <- factor(res$BEST.GWAS.status, levels = c('Other', 'Risk Locus'))
    
    return(res)
}


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

table(region_twas_z$BEST.GWAS.status)
## Initial version
# Other DLPFC HIPPO  Both
# 92402  2675  1354  1696
## Latest
# Other Risk Locus
# 92402       5725



pdf('pdf/twas_z.pdf', useDingbats = FALSE, width = 24, height = 14)
ggplot(region_twas_z,
    aes(x = DLPFC, y = HIPPO, color = FDR.5perc, shape = in_both)) +
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
#   None   3453 2436 5889
#   DLPFC   136   64  200
#   HIPPO    63   51  114
#   Both      0  101  101
#   Sum    3652 2652 6304
#
# $gene$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE Sum
#   None    150   90 240
#   DLPFC    60    8  68
#   HIPPO    13    5  18
#   Both      0   37  37
#   Sum     223  140 363
#
#
# $exon
# $exon$Other
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None  32859 11911 44770
#   DLPFC  1440   292  1732
#   HIPPO   676   238   914
#   Both      0   535   535
#   Sum   34975 12976 47951
#
# $exon$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE  Sum
#   None   1542  467 2009
#   DLPFC   505   66  571
#   HIPPO   235   52  287
#   Both      0  219  219
#   Sum    2282  804 3086
#
#
# $jxn
# $jxn$Other
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None  17741  7262 25003
#   DLPFC   650   158   808
#   HIPPO   441   138   579
#   Both      0   359   359
#   Sum   18832  7917 26749
#
# $jxn$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE  Sum
#   None    706  316 1022
#   DLPFC   214   26  240
#   HIPPO   132   34  166
#   Both      0  145  145
#   Sum    1052  521 1573
#
#
# $tx
# $tx$Other
#        In both
# FDR <5% FALSE  TRUE   Sum
#   None   7039  3520 10559
#   DLPFC   295    76   371
#   HIPPO   181    87   268
#   Both      0   200   200
#   Sum    7515  3883 11398
#
# $tx$`Risk Locus`
#        In both
# FDR <5% FALSE TRUE Sum
#   None    304  146 450
#   DLPFC    97   17 114
#   HIPPO    71    8  79
#   Both      0   60  60
#   Sum     472  231 703
  

## Subset by significant (TWAS FDR < 5%)
ttSig <- map(split(tt, tt$region), ~ .x[.x$TWAS.FDR < 0.05, ])
map(ttSig, dim)
# $DLPFC
# [1] 5760   27
#
# $HIPPO
# [1] 4081   27

map_int(ttSig, ~ length(unique(.x$geneid)))
# DLPFC HIPPO
#  1514  1255

## Original numbers when computing FDR across all 4 features at the same time:
# DLPFC HIPPO
#  1519  1256

map_int(ttSig, ~ length(unique(.x$geneid[.x$TWAS.P < 5e-08])))
# DLPFC HIPPO
#   115    84
## Just checking...
map_int(split(tt, tt$region), ~ length(unique(.x$geneid[.x$TWAS.P < 5e-08])))
# DLPFC HIPPO
#   115    84

## save for later
save(tt, ttSig, get_variable_by_region, ttReg_map, region_twas_z, file = 'rda/tt_objects.Rdata')







## Continue


## Check correlations among FDR corrected p-values between TWAS
## and either BEST GWAS FDR or EQTL GWAS FDR
check_cor <- function(x, y) {
    cor(-log10(x), -log10(y))
}
with(tt, check_cor(TWAS.FDR, BEST.GWAS.FDR.computed))
# [1] 0.425936
with(tt, check_cor(TWAS.FDR, EQTL.FDR.computed))
# [1] 0.8243939

map_dbl(ttSig, ~ with(.x, check_cor(TWAS.FDR, BEST.GWAS.FDR.computed)))
#     DLPFC     HIPPO
# 0.5629601 0.5427811
map_dbl(ttSig, ~ with(.x, check_cor(TWAS.FDR, EQTL.FDR.computed)))
#     DLPFC     HIPPO
# 0.7460110 0.7498989


map_dbl(split(tt, tt$BEST.GWAS.status), ~ with(.x, check_cor(TWAS.FDR, BEST.GWAS.FDR.computed)))
#     Index     Other     Proxy
# 0.1169518 0.3485689 0.2396344
map_dbl(split(tt, tt$BEST.GWAS.status), ~ with(.x, check_cor(TWAS.FDR, EQTL.FDR.computed)))
#     Index     Other     Proxy
# 0.7901372 0.8054720 0.8109420

map_dbl(split(tt, tt$BEST.GWAS.status), ~ with(.x, check_cor(TWAS.P, BEST.GWAS.P.computed)))
#     Index     Other     Proxy
# 0.1046549 0.3369664 0.2392657
map_dbl(split(tt, tt$BEST.GWAS.status), ~ with(.x, check_cor(TWAS.P, EQTL.P.computed)))
#     Index     Other     Proxy
# 0.7951933 0.7947084 0.8124684

tt_sigonly <- tt[tt$TWAS.FDR < 0.05, ]
ggplot(tt_sigonly, aes(
    x = -log10(TWAS.FDR),
    y = -log10(BEST.GWAS.FDR.computed),
    color = TWAS.P < 5e-08
)) + geom_point() + facet_grid(region ~ feature)

ggplot(tt_sigonly, aes(
    x = -log10(TWAS.P),
    y = -log10(BEST.GWAS.P.computed),
    color = TWAS.P < 5e-08
)) + geom_point() + facet_grid(region ~ feature)

ggplot(tt_sigonly, aes(
    x = TWAS.Z,
    y = BEST.GWAS.Z,
    color = TWAS.P < 5e-08
)) + geom_point() + facet_grid(region ~ feature)

ggplot(tt_sigonly, aes(
    x = -log10(TWAS.FDR),
    y = -log10(EQTL.FDR.computed),
    color = TWAS.P < 5e-08
)) + geom_point() + facet_grid(region ~ feature)

ggplot(tt_sigonly, aes(
    x = -log10(TWAS.P),
    y = -log10(EQTL.P.computed),
    color = TWAS.P < 5e-08
)) + geom_point() + facet_grid(region ~ feature)

ggplot(tt_sigonly, aes(
    x = TWAS.Z,
    y = EQTL.GWAS.Z,
    color = TWAS.P < 5e-08
)) + geom_point() + facet_grid(region ~ feature)

map(ttSig, ~ addmargins(table(
    'TWAS P < 5e-08' = .x$TWAS.P < 5e-08,
    'BEST GWAS P < 5e-08' = .x$BEST.GWAS.P.computed < 5e-08
)))
# $DLPFC
#               BEST GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE  Sum
#          FALSE  4067 1327 5394
#          TRUE      9  357  366
#          Sum    4076 1684 5760
#
# $HIPPO
#               BEST GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE  Sum
#          FALSE  2933  909 3842
#          TRUE      3  236  239
#          Sum    2936 1145 4081

map(ttSig, ~ addmargins(table(
    'TWAS P < 5e-08' = .x$TWAS.P < 5e-08,
    'EQTL GWAS P < 5e-08' = .x$EQTL.P.computed < 5e-08
)))
# $DLPFC
#               EQTL GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE  Sum
#          FALSE  5218  176 5394
#          TRUE    128  238  366
#          Sum    5346  414 5760
#
# $HIPPO
#               EQTL GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE  Sum
#          FALSE  3707  135 3842
#          TRUE     62  177  239
#          Sum    3769  312 4081

map(ttSig, ~ map(split(.x, .x$feature), ~ addmargins(table(
    'TWAS P < 5e-08' = .x$TWAS.P < 5e-08,
    'BEST GWAS P < 5e-08' = .x$BEST.GWAS.P.computed < 5e-08
))))
# $DLPFC
# $DLPFC$exon
#               BEST GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE  Sum
#          FALSE  2117  735 2852
#          TRUE      5  200  205
#          Sum    2122  935 3057
#
# $DLPFC$gene
#               BEST GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE Sum
#          FALSE   285   97 382
#          TRUE      1   23  24
#          Sum     286  120 406
#
# $DLPFC$jxn
#               BEST GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE  Sum
#          FALSE  1122  340 1462
#          TRUE      3   87   90
#          Sum    1125  427 1552
#
# $DLPFC$tx
#               BEST GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE Sum
#          FALSE   543  155 698
#          TRUE      0   47  47
#          Sum     543  202 745
#
#
# $HIPPO
# $HIPPO$exon
#               BEST GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE  Sum
#          FALSE  1377  450 1827
#          TRUE      1  127  128
#          Sum    1378  577 1955
#
# $HIPPO$gene
#               BEST GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE Sum
#          FALSE   206   43 249
#          TRUE      0   21  21
#          Sum     206   64 270
#
# $HIPPO$jxn
#               BEST GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE  Sum
#          FALSE   906  290 1196
#          TRUE      2   51   53
#          Sum     908  341 1249
#
# $HIPPO$tx
#               BEST GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE Sum
#          FALSE   444  126 570
#          TRUE      0   37  37
#          Sum     444  163 607


map(ttSig, ~ map(split(.x, .x$feature), ~ addmargins(table(
    'TWAS P < 5e-08' = .x$TWAS.P < 5e-08,
    'EQTL GWAS P < 5e-08' = .x$EQTL.P.computed < 5e-08
))))
# $DLPFC
# $DLPFC$exon
#               EQTL GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE  Sum
#          FALSE  2763   89 2852
#          TRUE     87  118  205
#          Sum    2850  207 3057
#
# $DLPFC$gene
#               EQTL GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE Sum
#          FALSE   371   11 382
#          TRUE      9   15  24
#          Sum     380   26 406
#
# $DLPFC$jxn
#               EQTL GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE  Sum
#          FALSE  1409   53 1462
#          TRUE     20   70   90
#          Sum    1429  123 1552
#
# $DLPFC$tx
#               EQTL GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE Sum
#          FALSE   675   23 698
#          TRUE     12   35  47
#          Sum     687   58 745
#
#
# $HIPPO
# $HIPPO$exon
#               EQTL GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE  Sum
#          FALSE  1751   76 1827
#          TRUE     36   92  128
#          Sum    1787  168 1955
#
# $HIPPO$gene
#               EQTL GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE Sum
#          FALSE   244    5 249
#          TRUE      8   13  21
#          Sum     252   18 270
#
# $HIPPO$jxn
#               EQTL GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE  Sum
#          FALSE  1160   36 1196
#          TRUE     12   41   53
#          Sum    1172   77 1249
#
# $HIPPO$tx
#               EQTL GWAS P < 5e-08
# TWAS P < 5e-08 FALSE TRUE Sum
#          FALSE   552   18 570
#          TRUE      6   31  37
#          Sum     558   49 607

## Venn diagrams of features by region, then joint (grouped by gene id)

## Compare TWAS Z-scores across DLPFC and HIPPO

## Read in the 179 CLOZUK+PGC2 snps

## Use the raggr output to find the proxy snps


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
