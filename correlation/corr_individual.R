## Usage:
# qrsh -l mem_free=100G,h_vmem=100G,h_fsize=100G
# Rscript corr_individual.R  > logs/corr_individual.txt 2>&1

library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('scales')
library('derfinder')
library('clusterProfiler')

dir.create('pdf', showWarnings = FALSE)
dir.create('rda', showWarnings = FALSE)

## Load expression data (cleaned and uncleaned) + rse objects for the pheno data
load('rda/expr_and_cleaned.Rdata', verbose = TRUE)
load('rda/rse_and_modQsva.Rdata', verbose = TRUE)

## Function for computing paired cors
paircor <- function(x, y) {
    stopifnot(nrow(x) == nrow(y))
    stopifnot(ncol(x) == ncol(y))
    sapply(seq_len(ncol(x)), function(i) {
        cor(x[, i], y[, i])
    })
}

## Test paircor
test <- cor(cleaned[['DLPFC']][['geneRpkm']][1:100, ], cleaned[['HIPPO']][['geneRpkm']][1:100, ])
test2 <- paircor(cleaned[['DLPFC']][['geneRpkm']][1:100, ], cleaned[['HIPPO']][['geneRpkm']][1:100, ])
stopifnot(all(diag(test) - test2 == 0))

computecor <- function(exp) {
    sets <- names(exp[[1]])
    res <- lapply(sets, function(feature) {
        message(paste(Sys.time(), 'processing feature', feature))
        paircor(exp[['DLPFC']][[feature]], exp[['HIPPO']][[feature]])
    })
    names(res) <- sets
    return(res)
}

## Compute individual level correlations across features
indv_expr <- computecor(expr)
indv_cleaned <- computecor(cleaned)
save(indv_expr, indv_cleaned, file = 'rda/indv_corr.Rdata')


pdf('pdf/indv_hist.pdf', useDingbats = FALSE)
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = 'light blue')
}, indv_expr, names(indv_expr), 'expr')
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = "#009E73")
}, indv_cleaned, names(indv_cleaned), 'cleaned expr (keeping Dx)')
dev.off()

pdf('pdf/indv_expr_vs_cleaned.pdf', useDingbats = FALSE)
mapply(function(x, y, set) {
    m <- !is.na(x) & !is.na(y)
    plot(x[m], y[m], xlab = 'expr', ylab = 'cleaned', main = set, pch = 19, col = alpha("#D55E00", 1))
    abline(a = 0, b = 1, col = 'blue')
}, indv_expr, indv_cleaned, names(indv_expr))
dev.off()

pdf('pdf/indv_expr_vs_cleaned_sameLims.pdf', useDingbats = FALSE)
mapply(function(x, y, set) {
    m <- !is.na(x) & !is.na(y)
    lims <- c(min(c(x[m], y[m])), max(c(x[m], y[m])))
    plot(x[m], y[m], xlab = 'expr', ylab = 'cleaned', main = set, pch = 19, col = alpha("#D55E00", 1),
        xlim = lims, ylim = lims
    )
    abline(a = 0, b = 1, col = 'blue')
}, indv_expr, indv_cleaned, names(indv_expr))
dev.off()


pdf('pdf/indv_box_sczd.pdf', useDingbats = FALSE)
mapply(function(x, y, set, type) {
    m <- !is.na(x)
    dx <- colData(y)$Dx[m]
    dx <- factor(ifelse(dx == 'Schizo', 'SCZD', ifelse(dx == 'Control', 'Control', 'hmm')))
    f <- lm(x[m] ~ dx)
    p <- summary(f)$coef[2, 4]
    ylim <- range(x[m])
    boxplot(x[m] ~ dx, main = paste(type, '-', set, '\n p-value:', signif(p, 3)),
        xlab = 'SCZD diagnosis', outline = FALSE, yab = 'Correlation', ylim = ylim)
    points(x[m] ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = 'light blue')
}, indv_expr, simple_rse[[1]], names(indv_expr), 'expr')
mapply(function(x, y, set, type) {
    m <- !is.na(x)
    dx <- colData(y)$Dx[m]
    f <- lm(x[m] ~ dx)
    p <- summary(f)$coef[2, 4]
    ylim <- range(x[m])
    dx <- factor(ifelse(dx == 'Schizo', 'SCZD', ifelse(dx == 'Control', 'Control', 'hmm')))
    boxplot(x[m] ~ dx, main = paste(type, '-', set, '\n p-value:', signif(p, 3)),
        xlab = 'SCZD diagnosis', outline = FALSE, ylab = 'Correlation', ylim = ylim)
    points(x[m] ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = '#009E73')
}, indv_cleaned, simple_rse[[1]], names(indv_cleaned), 'cleaned expr (keeping Dx)')
dev.off()



## Explore covariates to see if low correlation is related to any main covariate
pdf('pdf/indv_corr_vs_covariates.pdf', useDingbats = FALSE)
lapply(c('Age'), function(var) {
    message(paste(Sys.time(), 'processing', var))
        pd <- colData(simple_rse[[1]][[1]])
        y <- pd[, var]
        if(is(y, 'CompressedNumericList') | is(y, 'CompressedIntegerList')) {
            y <- mean(y)
        } else if (is.character(y)) {
            y <- as.factor(y)
        }
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - expr', pch = 19, col = 'deepskyblue3')
        }, indv_expr, names(indv_expr))
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - cleaned', pch = 19, col = '#009E73')
        }, indv_cleaned, names(indv_cleaned))
        return(NULL)
})
lapply(c('Sex', 'Race'), function(var) {
    message(paste(Sys.time(), 'processing', var))
        pd <- colData(simple_rse[[1]][[1]])
        y <- pd[, var]
        if(is(y, 'CompressedNumericList') | is(y, 'CompressedIntegerList')) {
            y <- mean(y)
        } else if (is.character(y)) {
            y <- as.factor(y)
        }
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - expr')
        }, indv_expr, names(indv_expr))
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - cleaned')
        }, indv_cleaned, names(indv_cleaned))
        return(NULL)
})

lapply(c('snpPC1', 'snpPC2', 'RIN', 'rRNA_rate', 'totalAssignedGene', 'mitoRate', 'totalMapped', 'mitoMapped', 'overallMapRate'), function(var) {
    message(paste(Sys.time(), 'processing', var))
    lapply(names(simple_rse), function(region) {
        pd <- colData(simple_rse[[region]][[1]])
        y <- pd[, var]
        if(is(y, 'CompressedNumericList') | is(y, 'CompressedIntegerList')) {
            y <- mean(y)
        } else if (is.character(y)) {
            y <- as.factor(y)
        }
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = paste0(region, ': ', var), xlab = 'Corr - expr', pch = 19, col = 'deepskyblue3')
        }, indv_expr, names(indv_expr))
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = paste0(region, ': ', var), xlab = 'Corr - cleaned', pch = 19, col = '#009E73')
        }, indv_cleaned, names(indv_cleaned))
        return(NULL)
    })
})
dev.off()



## Explore expression across all genes for a few samples with either very low or very high correlation (expr, not cleaned)
pdf('pdf/indv_expr_high_low_corr.pdf', useDingbats = FALSE)
lapply(names(indv_expr)[1], function(set) {
    i <- match(1:5, rank(indv_expr[[set]]))
    j <- match(261:265, rank(indv_expr[[set]]))
    
    mapply(function(k, rk) {
        plot(expr[[1]][[set]][, k], expr[[2]][[set]][, k], xlab = 'DLPFC', ylab = 'HIPPO', main = paste(set, 'rank:', rk, '- expr'), pch = 19, col = alpha('deepskyblue3', 1/10))
        abline(a = 0, b = 1, col = 'grey80')
        
        plot(cleaned[[1]][[set]][, k], cleaned[[2]][[set]][, k], xlab = 'DLPFC', ylab = 'HIPPO', main = paste(set, 'rank:', rk, '- cleaned expr'), pch = 19, col = alpha('#009E73', 1/10))
        abline(a = 0, b = 1, col = 'grey80')
        
    }, c(i, j), c(1:5, 261:265))
    return(NULL)
})
dev.off()


## Get KEGG sets
## code by Andrew Jaffe
KEGG_DATA <- clusterProfiler:::get_data_from_KEGG_db("hsa")
keggList = get("PATHID2EXTID", envir = KEGG_DATA)
keggTerms <- DOSE:::TERM2NAME(names(keggList), KEGG_DATA)

## Split genes by KEGG pathway
kegg_df <- data.frame(
    'entrez' = unlist(keggList),
    'kegg_path' = rep(names(keggList), sapply(keggList, length)),
    'kegg_term' = rep(keggTerms, sapply(keggList, length)),
    stringsAsFactors = FALSE
)

m <- match(rowRanges(simple_rse[['DLPFC']][['gene']])$EntrezID, as.integer(kegg_df$entrez))
m1 <- !is.na(m)
table(m1)
m2 <- m[m1]

## 1,420 KEGG entrez ids are missing
length(unique(kegg_df$entrez))
m_rev <- match(as.integer(unique(kegg_df$entrez)), rowRanges(simple_rse[['DLPFC']][['gene']])$EntrezID)
m1_rev <- !is.na(m_rev)
table(m1_rev)

kegg_df$m_gene <- match(as.integer(kegg_df$entrez), rowRanges(simple_rse[['DLPFC']][['gene']])$EntrezID)
kegg_df_comp <- kegg_df[!is.na(kegg_df$m_gene), ]

kegg_genes <- split(kegg_df_comp$m_gene, kegg_df_comp$kegg_path)
length(kegg_genes)
kegg_genes <- kegg_genes[sapply(kegg_genes, length) > 1]
length(kegg_genes)


computecor_kegg <- function(exp) {
    pathways <- names(kegg_genes)
    res <- lapply(pathways, function(pathway) {
        message(paste(Sys.time(), 'processing pathway', pathway))
        paircor(exp[['DLPFC']][['geneRpkm']][kegg_genes[[pathway]], ], exp[['HIPPO']][['geneRpkm']][kegg_genes[[pathway]], ])
    })
    names(res) <- pathways
    return(res)
}

## Compute individual level correlations across features
indv_kegg_expr <- computecor_kegg(expr)
indv_kegg_cleaned <- computecor_kegg(cleaned)
save(indv_kegg_expr, indv_kegg_cleaned, kegg_df, kegg_df_comp, kegg_genes, file = 'rda/indv_kegg_corr.Rdata')


pdf('pdf/indv_mean_kegg.pdf', useDingbats = FALSE)
plot(sapply(indv_kegg_expr, mean), sapply(indv_kegg_cleaned, mean), xlab = 'expr', ylab = 'cleaned', main = paste('mean correlation for', length(kegg_genes), 'KEGG pathways'), pch = 19, col = 'darkviolet')
hist(sapply(indv_kegg_expr, mean), main = 'mean KEGG corr - expr', col = 'deepskyblue3')
hist(sapply(indv_kegg_cleaned, mean), main = 'mean KEGG corr - cleaned', col = '#009E73')
dev.off()

top_kegg <- function(kegg_cor, n = 10) {
    top <- head(sort(sapply(kegg_cor, mean), decreasing = TRUE), n = n)
    res <- cbind(kegg_df_comp[match(names(top), kegg_df_comp$kegg_path), ], mean_corr = top)
    res$n_genes <- sapply(res$kegg_path, function(x) { sum(kegg_df_comp$kegg_path == x) })
    rownames(res) <- NULL
    return(res[, c('kegg_path', 'kegg_term', 'mean_corr', 'n_genes')])
}
options(width = 120)
print('expr')
top_kegg(indv_kegg_expr)
print('cleaned')
top_kegg(indv_kegg_cleaned)

# > print('expr')
# [1] "expr"
# > top_kegg(indv_kegg_expr)
#    kegg_path                                           kegg_term mean_corr n_genes
# 1   hsa00400 Phenylalanine, tyrosine and tryptophan biosynthesis 0.9949677       3
# 2   hsa03450                          Non-homologous end-joining 0.9868043      11
# 3   hsa00785                              Lipoic acid metabolism 0.9846689       3
# 4   hsa04744                                   Phototransduction 0.9827878      12
# 5   hsa00740                               Riboflavin metabolism 0.9815958       7
# 6   hsa00524                 Butirosin and neomycin biosynthesis 0.9808939       4
# 7   hsa00471              D-Glutamine and D-glutamate metabolism 0.9808905       4
# 8   hsa04114                                      Oocyte meiosis 0.9800333      95
# 9   hsa03010                                            Ribosome 0.9795631      86
# 10  hsa01040             Biosynthesis of unsaturated fatty acids 0.9795592      21
# > print('cleaned')
# [1] "cleaned"
# > top_kegg(indv_kegg_cleaned)
#    kegg_path                                           kegg_term mean_corr n_genes
# 1   hsa00780                                   Biotin metabolism 1.0000000       2
# 2   hsa00785                              Lipoic acid metabolism 0.9813169       3
# 3   hsa00400 Phenylalanine, tyrosine and tryptophan biosynthesis 0.9783005       3
# 4   hsa00524                 Butirosin and neomycin biosynthesis 0.9763295       4
# 5   hsa00740                               Riboflavin metabolism 0.9720471       7
# 6   hsa04614                            Renin-angiotensin system 0.9648249       9
# 7   hsa03450                          Non-homologous end-joining 0.9602165      11
# 8   hsa00860                Porphyrin and chlorophyll metabolism 0.9554085      23
# 9   hsa00190                           Oxidative phosphorylation 0.9496137     122
# 10  hsa00062                               Fatty acid elongation 0.9489531       8


kegg_sczd <- function(kegg_cor, type) {
    dx <- colData(simple_rse[['DLPFC']][['gene']])$Dx
    dx <- factor(ifelse(dx == 'Schizo', 'SCZD', ifelse(dx == 'Control', 'Control', 'hmm')))
    # > levels(dx)
    # [1] "Control" "SCZD"
    
    res <- do.call(rbind, mapply(function(pathway, pathname) {
        f <- lm(pathway ~ dx)
        p <- summary(f)$coef[2, 4]
        ylim <- range(pathway)
        
        boxplot(pathway ~ dx, main = paste0(pathname, ': ',
            kegg_df_comp$kegg_term[which(kegg_df_comp$kegg_path == pathname)[1]],
            '\np-value: ', signif(p, 3)),
            ylim = ylim,
            xlab = 'SCZD diagnosis', outline = FALSE, ylab = paste('Correlation -', type))
        points(pathway ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = ifelse(type == 'expr', 'deepskyblue3', '#009E73'))
                
        c(summary(f)$coef[2, ], mean_sczd = mean(pathway[dx == 'SCZD']), mean_control = mean(pathway[dx == 'Control']))
    }, kegg_cor, names(kegg_cor), SIMPLIFY = FALSE))
    m <- match(rownames(res), kegg_df_comp$kegg_path)
    res <- cbind(res, kegg_df_comp[m, c('kegg_path', 'kegg_term')])
    res$n_genes <- sapply(res$kegg_path, function(x) { sum(kegg_df_comp$kegg_path == x) })
    res$padj <- p.adjust(res[, 'Pr(>|t|)'], method = 'fdr')
    res <- res[order(res$padj, decreasing = FALSE), ]
    rownames(res) <- NULL
    
    return(res)
}

pdf('pdf/indv_box_sczd_kegg.pdf', useDingbats = FALSE)
indv_kegg_expr_sczd <- kegg_sczd(indv_kegg_expr, 'expr')
indv_kegg_cleaned_sczd <- kegg_sczd(indv_kegg_cleaned, 'cleaned')
dev.off()

nrow(indv_kegg_expr_sczd)
# [1] 227
table(indv_kegg_expr_sczd$padj < 0.05)
# FALSE
#   227
table(indv_kegg_cleaned_sczd$padj < 0.05)
# FALSE  TRUE
#   214    13
options(width = 200)
subset(indv_kegg_cleaned_sczd, padj < 0.05)
#        Estimate   Std. Error   t value     Pr(>|t|) mean_sczd mean_control kegg_path                               kegg_term n_genes         padj
# 1  -0.011095422 0.0023079374 -4.807506 2.573804e-06 0.8441098    0.8552052  hsa03320                  PPAR signaling pathway      51 0.0005842534
# 2  -0.004573640 0.0009932781 -4.604591 6.440982e-06 0.8916175    0.8961911  hsa04142                                Lysosome     109 0.0007310515
# 3  -0.014778589 0.0037914968 -3.897825 1.232383e-04 0.7534047    0.7681833  hsa04640              Hematopoietic cell lineage      42 0.0093250285
# 4  -0.009263283 0.0024614478 -3.763347 2.067723e-04 0.8204990    0.8297623  hsa05145                           Toxoplasmosis     105 0.0117343254
# 5  -0.004902039 0.0013356126 -3.670255 2.933653e-04 0.8842804    0.8891824  hsa04210                               Apoptosis      76 0.0133187844
# 6  -0.011181312 0.0032033927 -3.490459 5.652704e-04 0.7469369    0.7581182  hsa00511                Other glycan degradation      16 0.0169275367
# 7   0.009304367 0.0026565159  3.502470 5.414783e-04 0.8257343    0.8164300  hsa04740                  Olfactory transduction      41 0.0169275367
# 8  -0.015689972 0.0045146292 -3.475362 5.965652e-04 0.7248851    0.7405750  hsa05140                           Leishmaniasis      58 0.0169275367
# 9  -0.010150273 0.0029939956 -3.390210 8.056484e-04 0.7362409    0.7463911  hsa00071                  Fatty acid degradation      33 0.0203202440
# 10  0.009464357 0.0028201798  3.355941 9.076858e-04 0.9244389    0.9149745  hsa04744                       Phototransduction      12 0.0206044679
# 11 -0.002924630 0.0008956438 -3.265394 1.238149e-03 0.9067666    0.9096912  hsa00970             Aminoacyl-tRNA biosynthesis      41 0.0239204922
# 12 -0.004382762 0.0013447454 -3.259176 1.264519e-03 0.9349081    0.9392909  hsa01040 Biosynthesis of unsaturated fatty acids      21 0.0239204922
# 13 -0.015603131 0.0049429940 -3.156615 1.781964e-03 0.7893232    0.8049264  hsa00460              Cyanoamino acid metabolism       5 0.0311158305
save(indv_kegg_expr_sczd, indv_kegg_cleaned_sczd, file = 'rda/indv_kegg_sczd.Rdata')

## Re-loading if necessary
if(FALSE) {
    f <- dir('rda', full.names = TRUE)
    f <- f[!grepl('limma', f)]
    for(ff in f) load(ff, verbose = TRUE)
    rm(ff, f)
}

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.5.0 Patched (2018-04-30 r74679)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       US/Eastern
#  date     2018-08-23
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  acepack                1.4.1     2016-10-29 CRAN (R 3.5.0)
#  AnnotationDbi          1.42.1    2018-05-17 Bioconductor
#  assertthat             0.2.0     2017-04-11 CRAN (R 3.5.0)
#  backports              1.1.2     2017-12-13 CRAN (R 3.5.0)
#  base                 * 3.5.0     2018-05-02 local
#  base64enc              0.1-3     2015-07-28 CRAN (R 3.5.0)
#  bibtex                 0.4.2     2017-06-30 CRAN (R 3.5.0)
#  bindr                  0.1.1     2018-03-13 CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 CRAN (R 3.5.0)
#  Biobase              * 2.40.0    2018-05-02 Bioconductor
#  BiocGenerics         * 0.26.0    2018-05-03 Bioconductor
#  BiocParallel         * 1.14.2    2018-07-08 Bioconductor
#  biomaRt                2.36.1    2018-06-13 Bioconductor
#  Biostrings             2.48.0    2018-05-03 Bioconductor
#  bit                    1.1-14    2018-05-29 CRAN (R 3.5.0)
#  bit64                  0.9-7     2017-05-08 CRAN (R 3.5.0)
#  bitops                 1.0-6     2013-08-17 CRAN (R 3.5.0)
#  blob                   1.1.1     2018-03-25 CRAN (R 3.5.0)
#  BSgenome               1.48.0    2018-05-03 Bioconductor
#  bumphunter             1.22.0    2018-05-03 Bioconductor
#  checkmate              1.8.5     2017-10-24 CRAN (R 3.5.0)
#  cluster                2.0.7-1   2018-04-13 CRAN (R 3.5.0)
#  clusterProfiler      * 3.8.1     2018-06-13 Bioconductor
#  codetools              0.2-15    2016-10-05 CRAN (R 3.5.0)
#  colorout             * 1.2-0     2018-05-02 Github (jalvesaq/colorout@c42088d)
#  colorspace             1.3-2     2016-12-14 CRAN (R 3.5.0)
#  compiler               3.5.0     2018-05-02 local
#  cowplot                0.9.3     2018-07-15 CRAN (R 3.5.0)
#  crayon                 1.3.4     2017-09-16 CRAN (R 3.5.0)
#  data.table             1.11.4    2018-05-27 cran (@1.11.4)
#  datasets             * 3.5.0     2018-05-02 local
#  DBI                    1.0.0     2018-05-02 CRAN (R 3.5.0)
#  DelayedArray         * 0.6.2     2018-07-23 Bioconductor
#  derfinder            * 1.14.0    2018-05-03 Bioconductor
#  derfinderHelper        1.14.0    2018-05-03 Bioconductor
#  devtools             * 1.13.6    2018-06-27 CRAN (R 3.5.0)
#  digest                 0.6.15    2018-01-28 CRAN (R 3.5.0)
#  DO.db                  2.9       2018-05-03 Bioconductor
#  doRNG                  1.7.1     2018-06-22 CRAN (R 3.5.0)
#  DOSE                   3.6.1     2018-06-20 Bioconductor
#  dplyr                  0.7.6     2018-06-29 CRAN (R 3.5.0)
#  enrichplot             1.0.2     2018-06-13 Bioconductor
#  fastmatch              1.1-0     2017-01-28 CRAN (R 3.5.0)
#  fgsea                  1.6.0     2018-05-03 Bioconductor
#  foreach                1.4.4     2017-12-12 CRAN (R 3.5.0)
#  foreign                0.8-70    2017-11-28 CRAN (R 3.5.0)
#  Formula                1.2-3     2018-05-03 CRAN (R 3.5.0)
#  GenomeInfoDb         * 1.16.0    2018-05-03 Bioconductor
#  GenomeInfoDbData       1.1.0     2018-04-17 Bioconductor
#  GenomicAlignments      1.16.0    2018-05-03 Bioconductor
#  GenomicFeatures        1.32.0    2018-05-03 Bioconductor
#  GenomicFiles           1.16.0    2018-05-03 Bioconductor
#  GenomicRanges        * 1.32.6    2018-07-20 Bioconductor
#  ggforce                0.1.3     2018-07-07 CRAN (R 3.5.0)
#  ggplot2                3.0.0     2018-07-03 CRAN (R 3.5.0)
#  ggraph                 1.0.2     2018-07-07 CRAN (R 3.5.0)
#  ggrepel                0.8.0     2018-05-09 CRAN (R 3.5.0)
#  ggridges               0.5.0     2018-04-05 CRAN (R 3.5.0)
#  glue                   1.3.0     2018-07-17 CRAN (R 3.5.0)
#  GO.db                  3.6.0     2018-05-03 Bioconductor
#  GOSemSim               2.6.0     2018-05-03 Bioconductor
#  graphics             * 3.5.0     2018-05-02 local
#  grDevices            * 3.5.0     2018-05-02 local
#  grid                   3.5.0     2018-05-02 local
#  gridExtra              2.3       2017-09-09 CRAN (R 3.5.0)
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.5.0)
#  Hmisc                  4.1-1     2018-01-03 CRAN (R 3.5.0)
#  hms                    0.4.2     2018-03-10 CRAN (R 3.5.0)
#  htmlTable              1.12      2018-05-26 CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 CRAN (R 3.5.0)
#  htmlwidgets            1.2       2018-04-19 CRAN (R 3.5.0)
#  httpuv                 1.4.5     2018-07-19 CRAN (R 3.5.0)
#  httr                   1.3.1     2017-08-20 CRAN (R 3.5.0)
#  igraph                 1.2.1     2018-03-10 CRAN (R 3.5.0)
#  IRanges              * 2.14.10   2018-05-17 Bioconductor
#  iterators              1.0.10    2018-07-13 CRAN (R 3.5.0)
#  jaffelab             * 0.99.21   2018-05-03 Github (LieberInstitute/jaffelab@7ed0ab7)
#  knitr                  1.20      2018-02-20 CRAN (R 3.5.0)
#  later                  0.7.3     2018-06-08 CRAN (R 3.5.0)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.5.0)
#  latticeExtra           0.6-28    2016-02-09 CRAN (R 3.5.0)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.5.0)
#  limma                  3.36.2    2018-06-21 Bioconductor
#  locfit                 1.5-9.1   2013-04-20 CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 CRAN (R 3.5.0)
#  MASS                   7.3-50    2018-04-30 CRAN (R 3.5.0)
#  Matrix                 1.2-14    2018-04-13 CRAN (R 3.5.0)
#  matrixStats          * 0.54.0    2018-07-23 CRAN (R 3.5.0)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.5.0)
#  methods              * 3.5.0     2018-05-02 local
#  munsell                0.5.0     2018-06-12 CRAN (R 3.5.0)
#  nnet                   7.3-12    2016-02-02 CRAN (R 3.5.0)
#  parallel             * 3.5.0     2018-05-02 local
#  pillar                 1.3.0     2018-07-14 CRAN (R 3.5.0)
#  pkgconfig              2.0.1     2017-03-21 CRAN (R 3.5.0)
#  pkgmaker               0.27      2018-05-25 CRAN (R 3.5.0)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 CRAN (R 3.5.0)
#  prettyunits            1.0.2     2015-07-13 CRAN (R 3.5.0)
#  progress               1.2.0     2018-06-14 CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 CRAN (R 3.5.0)
#  purrr                  0.2.5     2018-05-29 CRAN (R 3.5.0)
#  qvalue                 2.12.0    2018-05-03 Bioconductor
#  R6                     2.2.2     2017-06-17 CRAN (R 3.5.0)
#  rafalib              * 1.0.0     2015-08-09 CRAN (R 3.5.0)
#  RColorBrewer           1.1-2     2014-12-07 CRAN (R 3.5.0)
#  Rcpp                   0.12.18   2018-07-23 CRAN (R 3.5.0)
#  RCurl                  1.95-4.11 2018-07-15 CRAN (R 3.5.0)
#  registry               0.5       2017-12-03 CRAN (R 3.5.0)
#  reshape2               1.4.3     2017-12-11 CRAN (R 3.5.0)
#  rlang                  0.2.1     2018-05-30 cran (@0.2.1)
#  rmote                * 0.3.4     2018-05-02 deltarho (R 3.5.0)
#  rngtools               1.3.1     2018-05-15 CRAN (R 3.5.0)
#  rpart                  4.1-13    2018-02-23 CRAN (R 3.5.0)
#  Rsamtools              1.32.2    2018-07-03 Bioconductor
#  RSQLite                2.1.1     2018-05-06 CRAN (R 3.5.0)
#  rstudioapi             0.7       2017-09-07 CRAN (R 3.5.0)
#  rtracklayer            1.40.3    2018-06-13 Bioconductor
#  rvcheck                0.1.0     2018-05-23 CRAN (R 3.5.0)
#  S4Vectors            * 0.18.3    2018-06-13 Bioconductor
#  scales               * 0.5.0     2017-08-24 CRAN (R 3.5.0)
#  segmented              0.5-3.0   2017-11-30 CRAN (R 3.5.0)
#  servr                  0.10      2018-05-30 CRAN (R 3.5.0)
#  splines                3.5.0     2018-05-02 local
#  stats                * 3.5.0     2018-05-02 local
#  stats4               * 3.5.0     2018-05-02 local
#  stringi                1.2.4     2018-07-20 CRAN (R 3.5.0)
#  stringr                1.3.1     2018-05-10 CRAN (R 3.5.0)
#  SummarizedExperiment * 1.10.1    2018-05-17 Bioconductor
#  survival               2.42-3    2018-04-16 CRAN (R 3.5.0)
#  tibble                 1.4.2     2018-01-22 CRAN (R 3.5.0)
#  tidyr                  0.8.1     2018-05-18 CRAN (R 3.5.0)
#  tidyselect             0.2.4     2018-02-26 CRAN (R 3.5.0)
#  tools                  3.5.0     2018-05-02 local
#  tweenr                 0.1.5     2016-10-10 CRAN (R 3.5.0)
#  units                  0.6-0     2018-06-09 CRAN (R 3.5.0)
#  UpSetR                 1.3.3     2017-03-21 CRAN (R 3.5.0)
#  utils                * 3.5.0     2018-05-02 local
#  VariantAnnotation      1.26.1    2018-07-04 Bioconductor
#  viridis                0.5.1     2018-03-29 CRAN (R 3.5.0)
#  viridisLite            0.3.0     2018-02-01 CRAN (R 3.5.0)
#  withr                  2.1.2     2018-03-15 CRAN (R 3.5.0)
#  xfun                   0.3       2018-07-06 CRAN (R 3.5.0)
#  XML                    3.98-1.12 2018-07-15 CRAN (R 3.5.0)
#  xtable                 1.8-2     2016-02-05 CRAN (R 3.5.0)
#  XVector                0.20.0    2018-05-03 Bioconductor
#  zlibbioc               1.26.0    2018-05-02 Bioconductor
