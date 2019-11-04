## Usage:
# qrsh -l mem_free=100G,h_vmem=100G,h_fsize=100G
# Rscript corr_individual_projectors.R  > logs/corr_individual_projectors.txt 2>&1

library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('scales')
library('derfinder')
library('clusterProfiler')
library('biomaRt')
library('grr') # https://stackoverflow.com/questions/21947086/how-do-i-match-all-occurrences-in-r


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
rm(test, test2)

computecor <- function(exp) {
    sets <- names(exp[[1]])
    res <- lapply(sets, function(feature) {
        message(paste(Sys.time(), 'processing feature', feature))
        paircor(exp[['DLPFC']][[feature]], exp[['HIPPO']][[feature]])
    })
    names(res) <- sets
    return(res)
}

## get projectors
load("/dcl01/lieber/ajaffe/Keri/SoLo_Flow/Syn_vs_Retro_Cre/rdas/exprsGenes_voom_pooledCounts_RetroVsSyn.rda", verbose = TRUE)

## get homologs
ensembl = useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),
               mart = ensembl)

## add mouse homologs
# res$ensemblID = rowData(rse_gene[rownames(res),])$ensemblID
outGene$hsapien_homolog <- MMtoHG$hsapiens_homolog_ensembl_gene[
	match(outGene$ensemblID, MMtoHG$ensembl_gene_id)]

## Keep only those with human homologs
outGene_human <- subset(outGene, hsapien_homolog != '')
nrow(outGene)
# [1] 17287
sum(outGene$adj.P.Val < 0.05)
# [1] 1343
nrow(outGene_human)
# [1] 13114
sum(outGene_human$adj.P.Val < 0.05)
# 1] 1273

## Match to BSP2
bsp2_matches <- lapply(names(simple_rse$DLPFC), function(feature) {
    row <- rowRanges(simple_rse$DLPFC[[feature]])
    if (feature == 'tx') {
        row$ensemblID <- gsub('\\..*', '', row$gene_id)
    }
    IntegerList(matches(outGene_human$hsapien_homolog, row$ensemblID, list = TRUE))
})
names(bsp2_matches) <- names(simple_rse$DLPFC)
bsp2_matches <- DataFrame(bsp2_matches)

## Check which are missing        
lapply(bsp2_matches, function(x) { addmargins(table(
    'Missing' = elementNROWS(x) == 0,
    'DE' = outGene_human$adj.P.Val < 0.05
)) })
# $gene
#        DE
# Missing FALSE  TRUE   Sum
#   FALSE 11001  1223 12224
#   TRUE    840    50   890
#   Sum   11841  1273 13114
#
# $exon
#        DE
# Missing FALSE  TRUE   Sum
#   FALSE 11198  1235 12433
#   TRUE    643    38   681
#   Sum   11841  1273 13114
#
# $jxn
#        DE
# Missing FALSE  TRUE   Sum
#   FALSE 10945  1213 12158
#   TRUE    896    60   956
#   Sum   11841  1273 13114
#
# $tx
#        DE
# Missing FALSE  TRUE   Sum
#   FALSE 10976  1197 12173
#   TRUE    865    76   941
#   Sum   11841  1273 13114

## Check which are missing at any level
bsp2_matches_any <- do.call(cbind, lapply(bsp2_matches,
    function(x) elementNROWS(x) > 0))
addmargins(table('Missing (any)' = rowSums(bsp2_matches_any) == 0,
    'DE' = outGene_human$adj.P.Val < 0.05))
#              DE
# Missing (any) FALSE  TRUE   Sum
#         FALSE 11484  1251 12735
#         TRUE    357    22   379
#         Sum   11841  1273 13114

## Keep those that have at least one feature match
outGene_bsp2 <- outGene_human[rowSums(bsp2_matches_any) > 0, ]
outGene_bsp2 <- cbind(outGene_bsp2, bsp2_matches[rowSums(bsp2_matches_any) > 0, ])

nrow(outGene_bsp2)
# [1] 12735
sum(outGene_bsp2$adj.P.Val < 0.05)
# [1] 1251

## Subset to only genes
expr_gene <- lapply(expr, '[', 'geneRpkm')
cleaned_gene <- lapply(cleaned, '[', 'geneRpkm')
simple_rse <- lapply(simple_rse, '[', 'gene')

## Separate into genes from DE in mouse vs non vs all
mouse_sets <- function(exp) {
    res <- lapply(1:4, function(i) {
        feature <- names(simple_rse$DLPFC)[i]
        all <- exp[[i]]
        
        ## Subset to those in the mouse data (or absent)
        mouse <- all[unlist(outGene_bsp2[[feature]]), ]
        notmouse <- all[-unlist(outGene_bsp2[[feature]]), ]
        
        mouse_up <- all[unlist(outGene_bsp2[[feature]][outGene_bsp2$logFC > 0]), ]
        mouse_down <- all[unlist(outGene_bsp2[[feature]][outGene_bsp2$logFC < 0]), ]
        
        mouse_up_de <- all[unlist(outGene_bsp2[[feature]][outGene_bsp2$logFC > 0 & outGene_bsp2$adj.P.Val < 0.05]), ]
        mouse_down_de <- all[unlist(outGene_bsp2[[feature]][outGene_bsp2$logFC < 0 & outGene_bsp2$adj.P.Val < 0.05]), ]
        
        result <- list(
            'mouse' = mouse,
            'human_only' = notmouse,
            'mouse_up' = mouse_up,
            'mouse_down' = mouse_down,
            'mouse_up_de' = mouse_up_de,
            'mouse_down_de' = mouse_down_de
        )
        names(result) <- paste0(feature, '_', names(result))
        return(result)
        
    })
    res <- do.call(c, res)
    return(res)
}

## Save the original data
expr_ori <- expr
cleaned_ori <- cleaned

## Now create the mouse subsets
expr <- lapply(expr_ori, mouse_sets)
cleaned <- lapply(cleaned_ori, mouse_sets)

sapply(expr[[1]], nrow)
#       gene_mouse    gene_human_only      gene_mouse_up    gene_mouse_down
#            12224              12577               6069               6155
# gene_mouse_up_de gene_mouse_down_de         exon_mouse    exon_human_only
#              590                633             320657              78370
#    exon_mouse_up    exon_mouse_down   exon_mouse_up_de exon_mouse_down_de
#           166872             153785              17503              15852
#        jxn_mouse     jxn_human_only       jxn_mouse_up     jxn_mouse_down
#           162968             135215              86210              76758
#  jxn_mouse_up_de  jxn_mouse_down_de           tx_mouse      tx_human_only
#             9924               8629              63928              29368
#      tx_mouse_up      tx_mouse_down     tx_mouse_up_de   tx_mouse_down_de
#            33070              30858               3643               3268

save(outGene_bsp2, file = 'rda/mouse_input_data.Rdata')

## Compute individual level correlations across features
indv_expr <- computecor(expr)
indv_cleaned <- computecor(cleaned)
save(indv_expr, indv_cleaned, file = 'rda/mouse_indv_corr.Rdata')


pdf('pdf/mouse_indv_hist.pdf', useDingbats = FALSE)
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = 'light blue')
}, indv_expr, names(indv_expr), 'expr')
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = "#009E73")
}, indv_cleaned, names(indv_cleaned), 'cleaned expr (keeping Dx)')
dev.off()

pdf('pdf/mouse_indv_expr_vs_cleaned.pdf', useDingbats = FALSE)
mapply(function(x, y, set) {
    m <- !is.na(x) & !is.na(y)
    plot(x[m], y[m], xlab = 'expr', ylab = 'cleaned', main = set, pch = 19, col = alpha("#D55E00", 1))
    abline(a = 0, b = 1, col = 'blue')
}, indv_expr, indv_cleaned, names(indv_expr))
dev.off()

pdf('pdf/mouse_indv_expr_vs_cleaned_sameLims.pdf', useDingbats = FALSE)
mapply(function(x, y, set) {
    m <- !is.na(x) & !is.na(y)
    lims <- c(min(c(x[m], y[m])), max(c(x[m], y[m])))
    plot(x[m], y[m], xlab = 'expr', ylab = 'cleaned', main = set, pch = 19, col = alpha("#D55E00", 1),
        xlim = lims, ylim = lims
    )
    abline(a = 0, b = 1, col = 'blue')
}, indv_expr, indv_cleaned, names(indv_expr))
dev.off()


pdf('pdf/mouse_indv_box_sczd.pdf', useDingbats = FALSE)
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
    points(x[m] ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = ifelse(dx == 'SCZD', 'aquamarine4', 'orchid4'))
}, indv_cleaned, simple_rse[[1]], names(indv_cleaned), 'cleaned expr (keeping Dx)')
dev.off()



## Explore covariates to see if low correlation is related to any main covariate
pdf('pdf/mouse_indv_corr_vs_covariates.pdf', useDingbats = FALSE)
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
pdf('pdf/mouse_indv_expr_high_low_corr.pdf', useDingbats = FALSE)
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

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.6.1 Patched (2019-09-06 r77160)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-10-08
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version    date       lib source
#  acepack                1.4.1      2016-10-29 [2] CRAN (R 3.6.1)
#  AnnotationDbi          1.46.1     2019-08-20 [1] Bioconductor
#  assertthat             0.2.1      2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.5      2019-10-02 [1] CRAN (R 3.6.1)
#  base64enc              0.1-3      2015-07-28 [2] CRAN (R 3.6.1)
#  bibtex                 0.4.2      2017-06-30 [1] CRAN (R 3.6.1)
#  Biobase              * 2.44.0     2019-05-02 [2] Bioconductor
#  BiocGenerics         * 0.30.0     2019-05-02 [1] Bioconductor
#  BiocManager            1.30.4     2018-11-13 [1] CRAN (R 3.6.1)
#  BiocParallel         * 1.18.1     2019-08-06 [1] Bioconductor
#  biomaRt              * 2.40.5     2019-10-01 [1] Bioconductor
#  Biostrings             2.52.0     2019-05-02 [1] Bioconductor
#  bit                    1.1-14     2018-05-29 [2] CRAN (R 3.6.1)
#  bit64                  0.9-7      2017-05-08 [2] CRAN (R 3.6.1)
#  bitops                 1.0-6      2013-08-17 [2] CRAN (R 3.6.1)
#  blob                   1.2.0      2019-07-09 [2] CRAN (R 3.6.1)
#  BSgenome               1.52.0     2019-05-02 [1] Bioconductor
#  bumphunter             1.27.1     2019-10-03 [1] Github (lcolladotor/bumphunter@308ea44)
#  callr                  3.3.2      2019-09-22 [1] CRAN (R 3.6.1)
#  checkmate              1.9.4      2019-07-04 [1] CRAN (R 3.6.1)
#  cli                    1.1.0      2019-03-19 [1] CRAN (R 3.6.1)
#  cluster                2.1.0      2019-06-19 [3] CRAN (R 3.6.1)
#  clusterProfiler      * 3.12.0     2019-05-02 [1] Bioconductor
#  codetools              0.2-16     2018-12-24 [3] CRAN (R 3.6.1)
#  colorout             * 1.2-2      2019-09-26 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1      2019-03-18 [2] CRAN (R 3.6.1)
#  cowplot                1.0.0      2019-07-11 [1] CRAN (R 3.6.1)
#  crayon                 1.3.4      2017-09-16 [1] CRAN (R 3.6.1)
#  curl                   4.2        2019-09-24 [1] CRAN (R 3.6.1)
#  data.table             1.12.4     2019-10-03 [1] CRAN (R 3.6.1)
#  DBI                    1.0.0      2018-05-02 [2] CRAN (R 3.6.1)
#  DelayedArray         * 0.10.0     2019-05-02 [2] Bioconductor
#  derfinder            * 1.18.9     2019-09-20 [1] Bioconductor
#  derfinderHelper        1.18.1     2019-05-22 [1] Bioconductor
#  desc                   1.2.0      2018-05-01 [2] CRAN (R 3.6.1)
#  devtools             * 2.2.1      2019-09-24 [1] CRAN (R 3.6.1)
#  digest                 0.6.21     2019-09-20 [1] CRAN (R 3.6.1)
#  DO.db                  2.9        2019-09-26 [1] Bioconductor
#  doRNG                  1.7.1      2018-06-22 [2] CRAN (R 3.6.1)
#  DOSE                   3.10.2     2019-06-24 [1] Bioconductor
#  dplyr                  0.8.3      2019-07-04 [1] CRAN (R 3.6.1)
#  ellipsis               0.3.0      2019-09-20 [1] CRAN (R 3.6.1)
#  enrichplot             1.4.0      2019-05-02 [1] Bioconductor
#  europepmc              0.3        2018-04-20 [1] CRAN (R 3.6.1)
#  farver                 1.1.0      2018-11-20 [1] CRAN (R 3.6.1)
#  fastmatch              1.1-0      2017-01-28 [1] CRAN (R 3.6.1)
#  fgsea                  1.10.1     2019-08-21 [1] Bioconductor
#  foreach                1.4.7      2019-07-27 [2] CRAN (R 3.6.1)
#  foreign                0.8-71     2018-07-20 [3] CRAN (R 3.6.1)
#  Formula                1.2-3      2018-05-03 [2] CRAN (R 3.6.1)
#  fs                     1.3.1      2019-05-06 [1] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.20.0     2019-05-02 [1] Bioconductor
#  GenomeInfoDbData       1.2.1      2019-09-09 [2] Bioconductor
#  GenomicAlignments      1.20.1     2019-06-18 [1] Bioconductor
#  GenomicFeatures        1.36.4     2019-07-10 [1] Bioconductor
#  GenomicFiles           1.20.0     2019-05-02 [1] Bioconductor
#  GenomicRanges        * 1.36.1     2019-09-06 [1] Bioconductor
#  ggforce                0.3.1      2019-08-20 [2] CRAN (R 3.6.1)
#  ggplot2                3.2.1      2019-08-10 [1] CRAN (R 3.6.1)
#  ggplotify              0.0.4      2019-08-06 [1] CRAN (R 3.6.1)
#  ggraph                 2.0.0      2019-09-02 [2] CRAN (R 3.6.1)
#  ggrepel                0.8.1      2019-05-07 [1] CRAN (R 3.6.1)
#  ggridges               0.5.1      2018-09-27 [1] CRAN (R 3.6.1)
#  glue                   1.3.1      2019-03-12 [1] CRAN (R 3.6.1)
#  GO.db                  3.8.2      2019-09-26 [1] Bioconductor
#  googledrive            1.0.0      2019-08-19 [1] CRAN (R 3.6.1)
#  GOSemSim               2.10.0     2019-05-02 [1] Bioconductor
#  graphlayouts           0.5.0      2019-08-20 [2] CRAN (R 3.6.1)
#  gridExtra              2.3        2017-09-09 [2] CRAN (R 3.6.1)
#  gridGraphics           0.4-1      2019-05-20 [1] CRAN (R 3.6.1)
#  grr                  * 0.9.5      2016-08-26 [1] CRAN (R 3.6.1)
#  gtable                 0.3.0      2019-03-25 [2] CRAN (R 3.6.1)
#  Hmisc                  4.2-0      2019-01-26 [1] CRAN (R 3.6.1)
#  hms                    0.5.1      2019-08-23 [2] CRAN (R 3.6.1)
#  htmlTable              1.13.2     2019-09-22 [1] CRAN (R 3.6.1)
#  htmltools              0.4.0      2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5        2019-10-04 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2      2019-09-11 [1] CRAN (R 3.6.1)
#  httr                   1.4.1      2019-08-05 [1] CRAN (R 3.6.1)
#  igraph                 1.2.4.1    2019-04-22 [2] CRAN (R 3.6.1)
#  IRanges              * 2.18.3     2019-09-24 [1] Bioconductor
#  iterators              1.0.12     2019-07-26 [2] CRAN (R 3.6.1)
#  jaffelab             * 0.99.29    2019-09-26 [1] Github (LieberInstitute/jaffelab@a7d87cb)
#  jsonlite               1.6        2018-12-07 [2] CRAN (R 3.6.1)
#  knitr                  1.25       2019-09-18 [1] CRAN (R 3.6.1)
#  later                  1.0.0      2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38    2018-11-04 [3] CRAN (R 3.6.1)
#  latticeExtra           0.6-28     2016-02-09 [2] CRAN (R 3.6.1)
#  lazyeval               0.2.2      2019-03-15 [2] CRAN (R 3.6.1)
#  lifecycle              0.1.0      2019-08-01 [1] CRAN (R 3.6.1)
#  limma                  3.40.6     2019-07-26 [1] Bioconductor
#  locfit                 1.5-9.1    2013-04-20 [2] CRAN (R 3.6.1)
#  magrittr               1.5        2014-11-22 [1] CRAN (R 3.6.1)
#  MASS                   7.3-51.4   2019-03-31 [3] CRAN (R 3.6.1)
#  Matrix                 1.2-17     2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0     2019-09-07 [1] CRAN (R 3.6.1)
#  memoise                1.1.0      2017-04-21 [2] CRAN (R 3.6.1)
#  munsell                0.5.0      2018-06-12 [2] CRAN (R 3.6.1)
#  nnet                   7.3-12     2016-02-02 [3] CRAN (R 3.6.1)
#  pillar                 1.4.2      2019-06-29 [1] CRAN (R 3.6.1)
#  pkgbuild               1.0.5      2019-08-26 [2] CRAN (R 3.6.1)
#  pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 3.6.1)
#  pkgload                1.0.2      2018-10-29 [2] CRAN (R 3.6.1)
#  pkgmaker               0.27       2018-05-25 [2] CRAN (R 3.6.1)
#  plyr                   1.8.4      2016-06-08 [2] CRAN (R 3.6.1)
#  png                    0.1-7      2013-12-03 [2] CRAN (R 3.6.1)
#  polyclip               1.10-0     2019-03-14 [2] CRAN (R 3.6.1)
#  prettyunits            1.0.2      2015-07-13 [1] CRAN (R 3.6.1)
#  processx               3.4.1      2019-07-18 [1] CRAN (R 3.6.1)
#  progress               1.2.2      2019-05-16 [1] CRAN (R 3.6.1)
#  promises               1.1.0      2019-10-04 [1] CRAN (R 3.6.1)
#  ps                     1.3.0      2018-12-21 [2] CRAN (R 3.6.1)
#  purrr                  0.3.2      2019-03-15 [2] CRAN (R 3.6.1)
#  qvalue                 2.16.0     2019-05-02 [1] Bioconductor
#  R6                     2.4.0      2019-02-14 [2] CRAN (R 3.6.1)
#  rafalib              * 1.0.0      2015-08-09 [1] CRAN (R 3.6.1)
#  RColorBrewer           1.1-2      2014-12-07 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.2      2019-07-25 [1] CRAN (R 3.6.1)
#  RCurl                  1.95-4.12  2019-03-04 [2] CRAN (R 3.6.1)
#  registry               0.5-1      2019-03-05 [2] CRAN (R 3.6.1)
#  remotes                2.1.0.9000 2019-09-26 [1] Github (MangoTheCat/remotes@2f1a040)
#  reshape2               1.4.3      2017-12-11 [2] CRAN (R 3.6.1)
#  rlang                  0.4.0      2019-06-25 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4      2019-09-26 [1] Github (cloudyr/rmote@fbce611)
#  rngtools               1.4        2019-07-01 [2] CRAN (R 3.6.1)
#  rpart                  4.1-15     2019-04-12 [3] CRAN (R 3.6.1)
#  rprojroot              1.3-2      2018-01-03 [2] CRAN (R 3.6.1)
#  Rsamtools              2.0.2      2019-10-02 [1] Bioconductor
#  RSQLite                2.1.2      2019-07-24 [2] CRAN (R 3.6.1)
#  rstudioapi             0.10       2019-03-19 [2] CRAN (R 3.6.1)
#  rtracklayer            1.44.4     2019-09-06 [1] Bioconductor
#  rvcheck                0.1.5      2019-10-01 [1] CRAN (R 3.6.1)
#  S4Vectors            * 0.22.1     2019-09-09 [1] Bioconductor
#  scales               * 1.0.0      2018-08-09 [2] CRAN (R 3.6.1)
#  segmented              1.0-0      2019-06-17 [2] CRAN (R 3.6.1)
#  servr                  0.15       2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo            1.1.1      2018-11-05 [1] CRAN (R 3.6.1)
#  stringi                1.4.3      2019-03-12 [2] CRAN (R 3.6.1)
#  stringr                1.4.0      2019-02-10 [1] CRAN (R 3.6.1)
#  SummarizedExperiment * 1.14.1     2019-07-31 [1] Bioconductor
#  survival               2.44-1.1   2019-04-01 [3] CRAN (R 3.6.1)
#  testthat               2.2.1      2019-07-25 [1] CRAN (R 3.6.1)
#  tibble                 2.1.3      2019-06-06 [1] CRAN (R 3.6.1)
#  tidygraph              1.1.2      2019-02-18 [1] CRAN (R 3.6.1)
#  tidyr                  1.0.0      2019-09-11 [1] CRAN (R 3.6.1)
#  tidyselect             0.2.5      2018-10-11 [2] CRAN (R 3.6.1)
#  triebeard              0.3.0      2016-08-04 [1] CRAN (R 3.6.1)
#  tweenr                 1.0.1      2018-12-14 [1] CRAN (R 3.6.1)
#  UpSetR                 1.4.0      2019-05-22 [1] CRAN (R 3.6.1)
#  urltools               1.7.3      2019-04-14 [1] CRAN (R 3.6.1)
#  usethis              * 1.5.1      2019-07-04 [1] CRAN (R 3.6.1)
#  VariantAnnotation      1.30.1     2019-05-19 [1] Bioconductor
#  vctrs                  0.2.0      2019-07-05 [1] CRAN (R 3.6.1)
#  viridis                0.5.1      2018-03-29 [2] CRAN (R 3.6.1)
#  viridisLite            0.3.0      2018-02-01 [2] CRAN (R 3.6.1)
#  withr                  2.1.2      2018-03-15 [2] CRAN (R 3.6.1)
#  xfun                   0.10       2019-10-01 [1] CRAN (R 3.6.1)
#  XML                    3.98-1.20  2019-06-06 [2] CRAN (R 3.6.1)
#  xml2                   1.2.2      2019-08-09 [2] CRAN (R 3.6.1)
#  xtable                 1.8-4      2019-04-21 [2] CRAN (R 3.6.1)
#  XVector                0.24.0     2019-05-02 [1] Bioconductor
#  zeallot                0.1.0      2018-01-28 [1] CRAN (R 3.6.1)
#  zlibbioc               1.30.0     2019-05-02 [2] Bioconductor
#
# [1] /users/lcollado/R/3.6
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6/R/3.6/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6/R/3.6/lib64/R/library
