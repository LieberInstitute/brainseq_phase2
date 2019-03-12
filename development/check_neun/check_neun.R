library('sessioninfo')
library('purrr')
library('jaffelab')
library('ggplot2')

message(paste(Sys.time(), 'loading ../rda/pcheck_both.Rdata'))
load('../rda/pcheck_both.Rdata', verbose = TRUE)

get_de <- function(x) {
    sign(x$F) == sign(x$span_F) & x$span_P.Value < 0.05 & x$P.Bonf < 0.01
}
pcheck_both$de <- get_de(pcheck_both)

## Rename for simplicity
dev <- pcheck_both

## For now:
dev <- dev[dev$type != 'exon', ]

# features <- c('gene', 'exon', 'jxn', 'tx')
features <- c('gene', 'jxn', 'tx')
top <- lapply(features, function(type) {
    f <- paste0('../rda/limma_dev_interaction_adjNeunProp_', type, '.Rdata')
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    top$type <- type
    top$P.Bonf <- p.adjust(top$P.Value, 'bonf')
    return(top)
})
names(top) <- features


neun <- do.call(rbind, map(top, function(x) {
    colnames(x) <- paste0('neun_', colnames(x))
    x
}))

stopifnot(identical(nrow(dev), nrow(neun)))

m <- match(rownames(dev), rownames(neun))
stopifnot(!any(is.na(m)))
neun <- neun[m, ]

dev <- cbind(dev, neun)
head(dev)


dev_type <- split(dev, dev$type)
tab_pbonf <- map(
    dev_type,
    ~ with(.x, table('Original Bonf<1%' = P.Bonf < 0.01, 'NeuN Bonf<1%' = neun_P.Bonf < 0.01))
)
map(tab_pbonf, addmargins)
# $gene
#                 NeuN Bonf<1%
# Original Bonf<1% FALSE  TRUE   Sum
#            FALSE  6905     2  6907
#            TRUE  17061   684 17745
#            Sum   23966   686 24652
#
# $jxn
#                 NeuN Bonf<1%
# Original Bonf<1%  FALSE   TRUE    Sum
#            FALSE 104732     62 104794
#            TRUE  191304   1083 192387
#            Sum   296036   1145 297181
#
# $tx
#                 NeuN Bonf<1%
# Original Bonf<1% FALSE  TRUE   Sum
#            FALSE 88902    15 88917
#            TRUE   3345   470  3815
#            Sum   92247   485 92732
map_dbl(tab_pbonf, getOR)
#       gene        jxn         tx
# 138.415685   9.562955 832.764126
map_dbl(tab_pbonf, ~ chisq.test(.x)$p.value)
#         gene          jxn           tx
# 3.859747e-60 2.838502e-99 0.000000e+00


tab_pbonf_span <- map(
    dev_type,
    ~ with(.x, table('Original Bonf<1% & Rep BrainSpan' = P.Bonf < 0.01 & span_P.Value < 0.05, 'NeuN Bonf<1% & Rep BrainSpan' = neun_P.Bonf < 0.01  & span_P.Value < 0.05))
)
map(tab_pbonf_span , addmargins)
# $gene
#                                 NeuN Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan FALSE  TRUE   Sum
#                            FALSE 13813     0 13813
#                            TRUE  10293   546 10839
#                            Sum   24106   546 24652
#
# $jxn
#                                 NeuN Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan  FALSE   TRUE    Sum
#                            FALSE 153227     59 153286
#                            TRUE  142957    938 143895
#                            Sum   296184    997 297181
#
# $tx
#                                 NeuN Bonf<1% & Rep BrainSpan
# Original Bonf<1% & Rep BrainSpan FALSE  TRUE   Sum
#                            FALSE 91012     5 91017
#                            TRUE   1415   300  1715
#                            Sum   92427   305 92732
map_dbl(tab_pbonf_span , getOR)
# gene        jxn         tx
#  Inf   17.04044 3859.16608
map_dbl(tab_pbonf_span , ~ chisq.test(.x)$p.value)
#          gene           jxn            tx
# 2.916268e-156 3.084074e-183  0.000000e+00

make_table <- function(ov) {
    
    res <- map_dfr(ov,
        ~ as.data.frame(matrix(as.vector(.x[1:2, 1:2]), nrow = 1, dimnames = list(1, c('Null_both', 'Original_only', 'NeuN_only', 'Both'))))
    )
    res$feature <- names(ov)
    res$OR <- map_dbl(ov, ~ getOR(.x[1:2, 1:2]))
    res$pval <- map_dbl(ov, ~ chisq.test(.x[1:2, 1:2])$p.value)
    res$pval_bonf <- p.adjust(res$pval, 'bonf')
    return(res)
}

options(width = 120)
make_table(tab_pbonf)
#   Null_both Original_only NeuN_only Both feature         OR         pval    pval_bonf
# 1      6905         17061         2  684    gene 138.415685 3.859747e-60 1.157924e-59
# 2    104732        191304        62 1083     jxn   9.562955 2.838502e-99 8.515505e-99
# 3     88902          3345        15  470      tx 832.764126 0.000000e+00 0.000000e+00
make_table(tab_pbonf_span)
#   Null_both Original_only NeuN_only Both feature         OR          pval     pval_bonf
# 1     13813         10293         0  546    gene        Inf 2.916268e-156 8.748805e-156
# 2    153227        142957        59  938     jxn   17.04044 3.084074e-183 9.252221e-183
# 3     91012          1415         5  300      tx 3859.16608  0.000000e+00  0.000000e+00

map_dbl(dev_type, ~ cor(.x$F, .x$neun_F))
#       gene        jxn         tx
# 0.08688325 0.16308335 0.70752450

## Compute the correlation on the scale that I'm actually plotting below
map_dbl(dev_type, ~ cor(log10(.x$F), log10(.x$neun_F)))
#      gene       jxn        tx
# 0.1507852 0.2738099 0.6006053

corrs <- cbind(map_dfr(dev_type, ~ map_dbl(split(.x, .x$de), ~ cor(log10(.x$F), log10(.x$neun_F)))), DE = c('FALSE', 'TRUE'))
corrs
#         gene       jxn        tx    DE
# 1 0.16026235 0.3079739 0.5696692 FALSE
# 2 0.06215385 0.1700112 0.5925908  TRUE


# ggplot(dev_type$gene, aes(x = F, y = neun_F, color = de)) + geom_point() + scale_x_log10() + scale_y_log10() + facet_grid( ~ de)

pdf('f_original_vs_f_adjNeuN_by_feature.pdf', width = 12, useDingbats = FALSE)
map2(dev_type, names(dev_type), function(df, type) {
    print(
        ggplot(df, aes(x = F, y = neun_F)) +
            geom_hex(aes(fill=..density..), bins = 100) +
            scale_x_log10() +
            scale_y_log10() +
            facet_grid( ~ de) +
            theme_bw(base_size = 30) +
            xlab('F-statistic: original') + 
            ylab('F-statistic: adj. NeuN prop') +
            labs(caption = 'Separated by DE status', title =  paste(type, 'corr =', paste(signif(corrs[, type], 3), collapse = ', ')))
    )
    return(invisible(NULL))
})
dev.off()

pdf('f_original_vs_f_adjNeuN.pdf', width = 12, useDingbats = FALSE, height = 18)
ggplot(dev, aes(x = F, y = neun_F)) +
    geom_hex(aes(fill=..density..), bins = 100) +
    scale_x_log10() +
    scale_y_log10() +
    facet_grid(type ~ de) +
    theme_bw(base_size = 30) +
    xlab('F-statistic: original') + 
    ylab('F-statistic: adj. NeuN prop') +
    labs(caption = 'Separated by DE status')
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
