## A cleaner and more focused script that explore_twas.R

library('tibble')
library('sessioninfo')
library('purrr')
library('readr')
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
tt$EQTL.GWAS.P.computed <- 2*(pnorm( abs(tt$EQTL.GWAS.Z ) , lower.tail=F ))

## Compute FDR by region for each feature 4 features
tt <- map_dfr(split(tt, tt$region), function(reg) {
    res <- map_dfr(split(reg, reg$feature), function(reg_feat) {
        reg_feat$TWAS.FDR <- p.adjust(reg_feat$TWAS.P, 'fdr')
        reg_feat$BEST.GWAS.FDR <- p.adjust(reg_feat$BEST.GWAS.P, 'fdr')
        reg_feat$BEST.GWAS.FDR.computed <- p.adjust(reg_feat$BEST.GWAS.P.computed, 'fdr')
        reg_feat$EQTL.GWAS.FDR.computed <- p.adjust(reg_feat$EQTL.GWAS.P.computed, 'fdr')
        return(reg_feat)
    })
    return(res[order(res$TWAS.P), ])
})
print(tt, width = 400)

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

map_int(ttSig, ~ length(unique(.x$geneid[.x$TWAS.P < 5e-08])))
# DLPFC HIPPO
#   115    84
## Just checking...
map_int(split(tt, tt$region), ~ length(unique(.x$geneid[.x$TWAS.P < 5e-08])))
# DLPFC HIPPO
#   115    84
 
## save for later
save(tt, ttSig, file = 'rda/tt_objects.Rdata')


## Continue


## Check correlations among FDR corrected p-values between TWAS
## and either BEST GWAS FDR or EQTL GWAS FDR
check_cor <- function(x, y) {
    cor(-log10(x), -log10(y))
}
with(tt, check_cor(TWAS.FDR, BEST.GWAS.FDR.computed))
# [1] 0.425936
with(tt, check_cor(TWAS.FDR, EQTL.GWAS.FDR.computed))
# [1] 0.8243939

map_dbl(ttSig, ~ with(.x, check_cor(TWAS.FDR, BEST.GWAS.FDR.computed)))
#     DLPFC     HIPPO
# 0.5629601 0.5427811
map_dbl(ttSig, ~ with(.x, check_cor(TWAS.FDR, EQTL.GWAS.FDR.computed)))
#     DLPFC     HIPPO
# 0.7460110 0.7498989


library('ggplot2')

plot_cor <- function(df, xvar, yvar) {
    xvarquo <- enquo(-log10(xvar))
    yvarquo <- enquo(-log10(yvar))
    ggplot(df, aes(!!xvarquo, !!yvarquo)) + geom_point()
}

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
    y = -log10(EQTL.GWAS.FDR.computed),
    color = TWAS.P < 5e-08
)) + geom_point() + facet_grid(region ~ feature)

ggplot(tt_sigonly, aes(
    x = -log10(TWAS.P),
    y = -log10(EQTL.GWAS.P.computed),
    color = TWAS.P < 5e-08
)) + geom_point() + facet_grid(region ~ feature)

ggplot(tt_sigonly, aes(
    x = TWAS.Z,
    y = EQTL.GWAS.Z,
    color = TWAS.P < 5e-08
)) + geom_point() + facet_grid(region ~ feature)


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
