## Usage:
# qrsh -l mem_free=100G,h_vmem=100G,h_fsize=100G
# Rscript corr_analysis.R  > logs/corr_analysis.txt 2>&1

library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('scales')
library('derfinder')
library('clusterProfiler')

dir.create('pdf', showWarnings = FALSE)
dir.create('rda', showWarnings = FALSE)

## Load expr data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata", verbose = TRUE)

regions <- c('DLPFC', 'HIPPO')
qinfo <- lapply(regions, function(region) {
    f <- paste0('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_', region, '.Rdata')
    message(paste(Sys.time(), 'loading file', f))
    load(f, verbose = TRUE)
    res <- list(
        'qsvBonf' = qsvBonf,
        'qSVs' = qSVs,
        'mod' = mod,
        'modQsva' = modQsva,
        'keepIndex' = keepIndex
    )
    return(res)
})
names(qinfo) <- regions

## Keep only samples that were observed twice
brains <- unlist(lapply(regions, function(region) {
    colData(rse_gene)$BrNum[qinfo[[region]]$keepIndex]
}))
both <- brains[which(duplicated(brains))]
stopifnot(length(both) == 265)


## Subset the rses for each brain region and re-order
simple_rse <- lapply(regions, function(region) {
    
    res <- list(
        'gene' = rse_gene[, qinfo[[region]]$keepIndex],
        'exon' = rse_exon[, qinfo[[region]]$keepIndex],
        'jxn' = rse_jxn[, qinfo[[region]]$keepIndex],
        'tx' = rse_tx[, qinfo[[region]]$keepIndex]
    )
    
    m <- match(both, colData(res$gene)$BrNum)
    stopifnot(all(!is.na(m)))
    ## Re-order by brain id
    lapply(res, function(x) {
        x[, m]
    })
})
names(simple_rse) <- regions

## Keep only dup brains on the mods and re-order
modQsva <- lapply(regions, function(region) {
    m <- match(both, colData(rse_gene[, qinfo[[region]]$keepIndex])$BrNum)
    stopifnot(all(!is.na(m)))
    
    qinfo[[region]]$modQsva[m, ]
})
names(modQsva) <- regions
stopifnot(all(sapply(modQsva, nrow) == 265))
stopifnot(all(sapply(modQsva, function(x) { colnames(x)[2] }) == 'DxSchizo'))

## Extract expr
expr <- lapply(simple_rse, function(rses) {
    res <- list(
        "geneRpkm" = assays(rses$gene)$rpkm,
        "exonRpkm" = assays(rses$exon)$rpkm,
        "jxnRp10m" = assays(rses$jxn)$rp10m,
        "txTpm" = assays(rses$tx)$tpm
    )
    lapply(res, function(x) { log2(x + 0.5) })
})

## Create cleaned expr versions protecting Dx
cleaned <- lapply(regions, function(region) {
    lapply(expr[[region]], cleaningY, modQsva[[region]], P = 2)
})
names(cleaned) <- regions
sapply(cleaned, function(x) sapply(x, dim) )


## Function for computing paired cors
paircor <- function(x, y) {
    stopifnot(nrow(x) == nrow(y))
    stopifnot(ncol(x) == ncol(y))
    sapply(seq_len(nrow(x)), function(i) {
        cor(x[i, ], y[i, ])
    })
}

## Test paircor
test <- cor(t(cleaned[['DLPFC']][['geneRpkm']][1:100, ]), t(cleaned[['HIPPO']][['geneRpkm']][1:100, ]))
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


corr_expr <- computecor(expr)
corr_cleaned <- computecor(cleaned)

message(paste(Sys.time(), 'saving rse and modQsva info'))
save(simple_rse, modQsva, file = 'rda/rse_and_modQsva.Rdata')
message(paste(Sys.time(), 'saving expression info'))
save(expr, cleaned, file = 'rda/expr_and_cleaned.Rdata')
message(paste(Sys.time(), 'saving corr results info'))
save(corr_expr, corr_cleaned, file = 'rda/corrs.Rdata')


pdf('pdf/hist.pdf', useDingbats = FALSE)
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = 'light blue')
}, corr_expr, names(corr_expr), 'expr')
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = "#009E73")
}, corr_cleaned, names(corr_cleaned), 'cleaned expr (keeping Dx)')
dev.off()

pdf('pdf/expr_vs_cleaned.pdf', useDingbats = FALSE)
mapply(function(x, y, set) {
    m <- !is.na(x) & !is.na(y)
    plot(x[m], y[m], xlab = 'expr', ylab = 'cleaned', main = set, pch = 19, col = alpha("#D55E00", 1/10))
}, corr_expr, corr_cleaned, names(corr_expr))
dev.off()



computecor_perm <- function(exp, n = 100, start_seed = 20180816) {
    stopifnot(n > 1)
    sets <- names(exp[[1]])[1] ## only gene level for now
    res <- lapply(sets, function(feature) {
        message(paste(Sys.time(), 'processing feature', feature))
        res_perm <- lapply(seq_len(n), function(i) {
            message(paste(Sys.time(), 'permutation', i))
            j <- nrow(exp[['HIPPO']][[feature]])
            set.seed(start_seed + i)
            paircor(exp[['DLPFC']][[feature]], exp[['HIPPO']][[feature]][sample(seq_len(j), j), ])
        })
        names(res_perm) <- paste0('permutation_', seq_len(n))
        res_perm <- do.call(cbind, res_perm)
        return(res_perm)
    })
    names(res) <- sets
    return(res)
}

corr_expr_perm <- computecor_perm(expr, n = 1000)
corr_cleaned_perm <- computecor_perm(cleaned, n = 1000)

## Hmmmm? Used for a previous version of computecor_perm that didn't cbind
# corr_expr_perm[[1]] <- do.call(cbind, corr_expr_perm[[1]])
# corr_cleaned_perm[[1]] <- do.call(cbind, corr_cleaned_perm[[1]])

message(paste(Sys.time(), 'saving corr perm results info'))
save(corr_expr_perm, corr_cleaned_perm, file = 'rda/corrs_perm.Rdata')

pdf('pdf/hist_perm.pdf', useDingbats = FALSE)
mapply(function(dat, set, type) {
    hist(as.vector(dat), main = paste(type, '-', set), col = 'light blue')
}, corr_expr_perm, names(corr_expr_perm), 'expr')
mapply(function(dat, set, type) {
    hist(as.vector(dat), main = paste(type, '-', set), col = "#009E73")
}, corr_cleaned_perm, names(corr_cleaned_perm), 'cleaned expr (keeping Dx)')
dev.off()


expr_pval <- mapply(edge.pvalue, corr_expr[1], corr_expr_perm, SIMPLIFY = FALSE)
cleaned_pval <- mapply(edge.pvalue, corr_cleaned[1], corr_cleaned_perm, SIMPLIFY = FALSE)

## FDR
expr_fdr <- lapply(expr_pval, p.adjust, method = 'fdr')
cleaned_fdr <- lapply(cleaned_pval, p.adjust, method = 'fdr')

## FWER
calc_fwer <- function(obs, perm) {
    derfinder:::.calcPval(obs, apply(perm, 2, max, na.rm = TRUE))
}
expr_fwer <- mapply(calc_fwer, corr_expr[1], corr_expr_perm, SIMPLIFY = FALSE)
cleaned_fwer <- mapply(calc_fwer, corr_cleaned[1], corr_cleaned_perm, SIMPLIFY = FALSE)

pdf('pdf/hist_pval.pdf', useDingbats = FALSE)
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = 'light blue')
}, expr_pval, names(expr_pval), 'expr')
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = "#009E73")
}, cleaned_pval, names(cleaned_pval), 'cleaned expr (keeping Dx)')
dev.off()

pdf('pdf/hist_fdr.pdf', useDingbats = FALSE)
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = 'light blue')
}, expr_fdr, names(expr_fdr), 'expr')
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = "#009E73")
}, cleaned_fdr, names(cleaned_fdr), 'cleaned expr (keeping Dx)')
dev.off()

pdf('pdf/hist_fwer.pdf', useDingbats = FALSE)
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = 'light blue')
}, expr_fwer, names(expr_fwer), 'expr')
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = "#009E73")
}, cleaned_fwer, names(cleaned_fwer), 'cleaned expr (keeping Dx)')
dev.off()



pdf('pdf/expr_vs_cleaned_pval.pdf', useDingbats = FALSE)
mapply(function(x, y, set) {
    m <- !is.na(x) & !is.na(y)
    plot(-log10(x[m]), -log10(y[m]), xlab = 'expr', ylab = 'cleaned', main = set, pch = 19, col = alpha("#D55E00", 1/10))
    abline(a = 0, b = 1, col = 'red')
}, expr_pval, cleaned_pval, names(expr_pval))
dev.off()

pdf('pdf/expr_vs_cleaned_fdr.pdf', useDingbats = FALSE)
mapply(function(x, y, set) {
    m <- !is.na(x) & !is.na(y)
    plot(-log10(x[m]), -log10(y[m]), xlab = 'expr', ylab = 'cleaned', main = set, pch = 19, col = alpha("#D55E00", 1/10))
    abline(a = 0, b = 1, col = 'red')
}, expr_fdr, cleaned_fdr, names(expr_fdr))
dev.off()

pdf('pdf/expr_vs_cleaned_fwer.pdf', useDingbats = FALSE)
mapply(function(x, y, set) {
    m <- !is.na(x) & !is.na(y)
    plot(-log10(x[m]), -log10(y[m]), xlab = 'expr', ylab = 'cleaned', main = set, pch = 19, col = alpha("#D55E00", 1/10))
    abline(a = 0, b = 1, col = 'red')
}, expr_fwer, cleaned_fwer, names(expr_fwer))
dev.off()


gene_pinfo <- data.frame(
    expr.pval = expr_pval[[1]],
    expr.fdr = expr_fdr[[1]],
    expr.fwer = expr_fwer[[1]],
    cleaned.pval = cleaned_pval[[1]],
    cleaned.fdr = cleaned_fdr[[1]],
    cleaned.fwer = cleaned_fwer[[1]],
    corr.expr = corr_expr[[1]],
    corr.cleaned = corr_cleaned[[1]]
)
rownames(gene_pinfo) <- rownames(rse_gene)

save(gene_pinfo, file = 'rda/gene_pinfo.Rdata')
save(expr_pval, expr_fdr, expr_fwer, cleaned_pval, cleaned_fdr, cleaned_fwer, file = 'rda/all_pinfo.Rdata')

options(width = 120)
print('Alpha = 0.1')
summary(gene_pinfo[, 1:6] < 0.1)
print('Alpha = 0.05')
summary(gene_pinfo[, 1:6] < 0.05)
print('Alpha = 0.01')
summary(gene_pinfo[, 1:6] < 0.01)
# > print('Alpha = 0.1')
# [1] "Alpha = 0.1"
# > summary(gene_pinfo[, 1:6] < 0.1)
#  expr.pval        expr.fdr       expr.fwer       cleaned.pval    cleaned.fdr     cleaned.fwer
#  Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical
#  FALSE:6610      FALSE:8403      FALSE:24168     FALSE:5041      FALSE:5565      FALSE:23060
#  TRUE :18042     TRUE :16249     TRUE :483       TRUE :19611     TRUE :19087     TRUE :1592
#                                  NA's :1
# > print('Alpha = 0.05')
# [1] "Alpha = 0.05"
# > summary(gene_pinfo[, 1:6] < 0.05)
#  expr.pval        expr.fdr       expr.fwer       cleaned.pval    cleaned.fdr     cleaned.fwer
#  Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical
#  FALSE:9665      FALSE:13037     FALSE:24312     FALSE:6453      FALSE:7162      FALSE:23722
#  TRUE :14987     TRUE :11615     TRUE :339       TRUE :18199     TRUE :17490     TRUE :930
#                                  NA's :1
# > print('Alpha = 0.01')
# [1] "Alpha = 0.01"
# > summary(gene_pinfo[, 1:6] < 0.01)
#  expr.pval        expr.fdr       expr.fwer       cleaned.pval    cleaned.fdr     cleaned.fwer
#  Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical   Mode :logical
#  FALSE:16268     FALSE:20540     FALSE:24627     FALSE:9489      FALSE:10485     FALSE:24462
#  TRUE :8384      TRUE :4112      TRUE :24        TRUE :15163     TRUE :14167     TRUE :190
#                                  NA's :1


## Load DE gene results
files <- c(
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV_matchHIPPO.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV_matchDLPFC.rda'
)

outGene <- lapply(files, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(outGene)
})
names(outGene) <- c('HIPPO_matchQSV', 'DLPFC_matchQSV')

## Load BrainSeq Phase 1 and Common Mind results
load('/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/caseControl/rdas/expressed_de_features.rda', verbose = TRUE)

prev <- list(
    'BSP1' = data.frame(
        ensemblID = names(outStatsExprs$Gene),
        adj.P.Val = outStatsExprs$Gene$fdr_qsva,
        logFC = outStatsExprs$Gene$log2FC_qsva,
        t = outStatsExprs$Gene$tstat_qsva
        ),
    'CMC' = data.frame(
        ensemblID = names(outStatsExprs$Gene),
        adj.P.Val = p.adjust(outStatsExprs$Gene$CMC_pval_qsva, method = 'fdr'),
        logFC = outStatsExprs$Gene$CMC_log2FC_qsva,
        t = outStatsExprs$Gene$CMC_tstat_qsva
        )
)
outGene <- c(outGene, prev)

## Add corr results
outGene <- lapply(outGene, function(out) {
    if('gencodeID' %in% colnames(out)) {
        m <- match(out$gencodeID, rownames(gene_pinfo))
    } else {
        m <- match(out$ensemblID, gsub('\\..*', '', rownames(gene_pinfo)))
    }
    cbind(out, gene_pinfo[m, ])
})

## Load exon/jx/tx results
outFeat <- lapply(c('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda', '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda'), function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    outTx$ensemblID <- gsub('\\..*', '', outTx$gene_id)
    return(list('exon' = outExon, 'jxn' = outJxn, 'tx' = outTx))
})
names(outFeat) <- c('DLPFC', 'HIPPO')

## Add corr results to features
outFeat <- lapply(outFeat, function(region) {
    lapply(region, function(out) {
        if('gencodeID' %in% colnames(out)) {
            m <- match(out$gencodeID, rownames(gene_pinfo))
        } else if('gencodeGeneID' %in% colnames(out)) {
            m <- match(out$gencodeGeneID, rownames(gene_pinfo))
        } else if('gene_id' %in% colnames(out)) {
            m <- match(out$gene_id, rownames(gene_pinfo))
        }
        cbind(out, gene_pinfo[m, ])
    })
})

save(outGene, outFeat, file = 'rda/out_info.Rdata')

compare_corr_de <- function(out, cutde = 0.05, cutvar = cutde) {
    out_comp <- out[ !is.na(out$adj.P.Val) & !is.na(out$expr.pval), ]
    res <- lapply(colnames(gene_pinfo)[1:6], function(var) {
        table(DE = out$adj.P.Val < cutde, Corr = out[, var] < cutvar)
    })
    names(res) <- colnames(gene_pinfo)[1:6]
    return(res)
}
calc_chi_pval <- function(x) {
    sapply(lapply(x, chisq.test), '[[', 'p.value')
}



corr_de_gene <- lapply(outGene, compare_corr_de)
print('Gene level tables for chisq')
corr_de_gene
# > print('Gene level tables for chisq')
# [1] "Gene level tables for chisq"
# > corr_de_gene
# $HIPPO_matchQSV
# $HIPPO_matchQSV$expr.pval
#        Corr
# DE      FALSE  TRUE
#   FALSE  9655 14949
#   TRUE     10    38
#
# $HIPPO_matchQSV$expr.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE 13021 11583
#   TRUE     16    32
#
# $HIPPO_matchQSV$expr.fwer
#        Corr
# DE      FALSE  TRUE
#   FALSE 24265   338
#   TRUE     47     1
#
# $HIPPO_matchQSV$cleaned.pval
#        Corr
# DE      FALSE  TRUE
#   FALSE  6452 18152
#   TRUE      1    47
#
# $HIPPO_matchQSV$cleaned.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE  7160 17444
#   TRUE      2    46
#
# $HIPPO_matchQSV$cleaned.fwer
#        Corr
# DE      FALSE  TRUE
#   FALSE 23675   929
#   TRUE     47     1
#
#
# $DLPFC_matchQSV
# $DLPFC_matchQSV$expr.pval
#        Corr
# DE      FALSE  TRUE
#   FALSE  9603 14804
#   TRUE     62   183
#
# $DLPFC_matchQSV$expr.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE 12942 11465
#   TRUE     95   150
#
# $DLPFC_matchQSV$expr.fwer
#        Corr
# DE      FALSE  TRUE
#   FALSE 24071   335
#   TRUE    241     4
#
# $DLPFC_matchQSV$cleaned.pval
#        Corr
# DE      FALSE  TRUE
#   FALSE  6417 17990
#   TRUE     36   209
#
# $DLPFC_matchQSV$cleaned.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE  7119 17288
#   TRUE     43   202
#
# $DLPFC_matchQSV$cleaned.fwer
#        Corr
# DE      FALSE  TRUE
#   FALSE 23480   927
#   TRUE    242     3
#
#
# $BSP1
# $BSP1$expr.pval
#        Corr
# DE      FALSE  TRUE
#   FALSE  6283 12843
#   TRUE     48   123
#
# $BSP1$expr.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE  8971 10155
#   TRUE     78    93
#
# $BSP1$expr.fwer
#        Corr
# DE      FALSE  TRUE
#   FALSE 18847   279
#   TRUE    171     0
#
# $BSP1$cleaned.pval
#        Corr
# DE      FALSE  TRUE
#   FALSE  3382 15744
#   TRUE     19   152
#
# $BSP1$cleaned.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE  3897 15229
#   TRUE     24   147
#
# $BSP1$cleaned.fwer
#        Corr
# DE      FALSE  TRUE
#   FALSE 18318   808
#   TRUE    169     2
#
#
# $CMC
# $CMC$expr.pval
#        Corr
# DE      FALSE  TRUE
#   FALSE  6214 12744
#   TRUE    117   222
#
# $CMC$expr.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE  8888 10070
#   TRUE    161   178
#
# $CMC$expr.fwer
#        Corr
# DE      FALSE  TRUE
#   FALSE 18682   276
#   TRUE    336     3
#
# $CMC$cleaned.pval
#        Corr
# DE      FALSE  TRUE
#   FALSE  3341 15617
#   TRUE     60   279
#
# $CMC$cleaned.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE  3847 15111
#   TRUE     74   265
#
# $CMC$cleaned.fwer
#        Corr
# DE      FALSE  TRUE
#   FALSE 18154   804
#   TRUE    333     6

print('Gene level chisq test pvals')
do.call(cbind, lapply(corr_de_gene, calc_chi_pval))
#              HIPPO_matchQSV DLPFC_matchQSV       BSP1        CMC
# expr.pval      0.0138237025   1.019910e-05 0.21361029 0.53773829
# expr.fdr       0.0101263445   1.176472e-05 0.79505488 0.86643488
# expr.fwer      1.0000000000   9.425265e-01 0.20437530 0.52004724
# cleaned.pval   0.0002763259   5.437789e-05 0.03199006 1.00000000
# cleaned.fdr    0.0002703122   9.062332e-05 0.05047347 0.52943731
# cleaned.fwer   0.8136789291   5.296141e-02 0.07316440 0.03467391

do.call(cbind, lapply(corr_de_gene, calc_chi_pval)) < 0.05
#              HIPPO_matchQSV DLPFC_matchQSV  BSP1   CMC
# expr.pval              TRUE           TRUE FALSE FALSE
# expr.fdr               TRUE           TRUE FALSE FALSE
# expr.fwer             FALSE          FALSE FALSE FALSE
# cleaned.pval           TRUE           TRUE  TRUE FALSE
# cleaned.fdr            TRUE           TRUE FALSE FALSE
# cleaned.fwer          FALSE          FALSE FALSE  TRUE

corr_de_feat <- lapply(outFeat, function(region) { lapply(region, compare_corr_de) })
print('Exon/jxn/tx level chisq table info')
corr_de_feat
# > print('Exon/jxn/tx level chisq table info')
# [1] "Exon/jxn/tx level chisq table info"
# > corr_de_feat
# $DLPFC
# $DLPFC$exon
# $DLPFC$exon$expr.pval
#        Corr
# DE       FALSE   TRUE
#   FALSE 127706 260019
#   TRUE     138    301
#
# $DLPFC$exon$expr.fdr
#        Corr
# DE       FALSE   TRUE
#   FALSE 185760 201965
#   TRUE     207    232
#
# $DLPFC$exon$expr.fwer
#        Corr
# DE       FALSE   TRUE
#   FALSE 383380   4344
#   TRUE     428     11
#
# $DLPFC$exon$cleaned.pval
#        Corr
# DE       FALSE   TRUE
#   FALSE  42826 344899
#   TRUE      92    347
#
# $DLPFC$exon$cleaned.fdr
#        Corr
# DE       FALSE   TRUE
#   FALSE  51620 336105
#   TRUE     101    338
#
# $DLPFC$exon$cleaned.fwer
#        Corr
# DE       FALSE   TRUE
#   FALSE 373360  14365
#   TRUE     429     10
#
#
# $DLPFC$jxn
# $DLPFC$jxn$expr.pval
#        Corr
# DE       FALSE   TRUE
#   FALSE  65383 128445
#   TRUE       5     16
#
# $DLPFC$jxn$expr.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE 94981 98847
#   TRUE      9    12
#
# $DLPFC$jxn$expr.fwer
#        Corr
# DE       FALSE   TRUE
#   FALSE 191568   2260
#   TRUE      21      0
#
# $DLPFC$jxn$cleaned.pval
#        Corr
# DE       FALSE   TRUE
#   FALSE  21041 172787
#   TRUE       2     19
#
# $DLPFC$jxn$cleaned.fdr
#        Corr
# DE       FALSE   TRUE
#   FALSE  25689 168139
#   TRUE       4     17
#
# $DLPFC$jxn$cleaned.fwer
#        Corr
# DE       FALSE   TRUE
#   FALSE 186715   7113
#   TRUE      21      0
#
#
# $DLPFC$tx
# $DLPFC$tx$expr.pval
#        Corr
# DE      FALSE  TRUE
#   FALSE 27472 54470
#   TRUE      1     5
#
# $DLPFC$tx$expr.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE 39463 42479
#   TRUE      1     5
#
# $DLPFC$tx$expr.fwer
#        Corr
# DE      FALSE  TRUE
#   FALSE 81124   817
#   TRUE      6     0
#
# $DLPFC$tx$cleaned.pval
#        Corr
# DE      FALSE  TRUE
#   FALSE 11374 70568
#   TRUE      1     5
#
# $DLPFC$tx$cleaned.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE 13363 68579
#   TRUE      2     4
#
# $DLPFC$tx$cleaned.fwer
#        Corr
# DE      FALSE  TRUE
#   FALSE 79002  2940
#   TRUE      6     0
#
#
#
# $HIPPO
# $HIPPO$exon
# $HIPPO$exon$expr.pval
#        Corr
# DE       FALSE   TRUE
#   FALSE  65382 128441
#   TRUE       6     20
#
# $HIPPO$exon$expr.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE 94978 98845
#   TRUE     12    14
#
# $HIPPO$exon$expr.fwer
#        Corr
# DE       FALSE   TRUE
#   FALSE 191563   2260
#   TRUE      26      0
#
# $HIPPO$exon$cleaned.pval
#        Corr
# DE       FALSE   TRUE
#   FALSE  21039 172784
#   TRUE       4     22
#
# $HIPPO$exon$cleaned.fdr
#        Corr
# DE       FALSE   TRUE
#   FALSE  25689 168134
#   TRUE       4     22
#
# $HIPPO$exon$cleaned.fwer
#        Corr
# DE       FALSE   TRUE
#   FALSE 186711   7112
#   TRUE      25      1
#
#
# $HIPPO$jxn
# $HIPPO$jxn$expr.pval
#        Corr
# DE       FALSE   TRUE
#   FALSE  65382 128441
#   TRUE       6     20
#
# $HIPPO$jxn$expr.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE 94978 98845
#   TRUE     12    14
#
# $HIPPO$jxn$expr.fwer
#        Corr
# DE       FALSE   TRUE
#   FALSE 191563   2260
#   TRUE      26      0
#
# $HIPPO$jxn$cleaned.pval
#        Corr
# DE       FALSE   TRUE
#   FALSE  21039 172784
#   TRUE       4     22
#
# $HIPPO$jxn$cleaned.fdr
#        Corr
# DE       FALSE   TRUE
#   FALSE  25689 168134
#   TRUE       4     22
#
# $HIPPO$jxn$cleaned.fwer
#        Corr
# DE       FALSE   TRUE
#   FALSE 186711   7112
#   TRUE      25      1
#
#
# $HIPPO$tx
# $HIPPO$tx$expr.pval
#        Corr
# DE      FALSE  TRUE
#   FALSE 27473 54475
#
# $HIPPO$tx$expr.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE 39464 42484
#
# $HIPPO$tx$expr.fwer
#        Corr
# DE      FALSE  TRUE
#   FALSE 81130   817
#
# $HIPPO$tx$cleaned.pval
#        Corr
# DE      FALSE  TRUE
#   FALSE 11375 70573
#
# $HIPPO$tx$cleaned.fdr
#        Corr
# DE      FALSE  TRUE
#   FALSE 13365 68583
#
# $HIPPO$tx$cleaned.fwer
#        Corr
# DE      FALSE  TRUE
#   FALSE 79008  2940



print('Exon/jxn/tx level chisq test pvals')
lapply(corr_de_feat, function(region) { do.call(cbind, lapply(region, calc_chi_pval)) } )
# $DLPFC
#                      exon       jxn        tx
# expr.pval    5.362391e-01 0.4648102 0.6582342
# expr.fdr     7.873276e-01 0.7300490 0.2562549
# expr.fwer    1.148746e-02 1.0000000 1.0000000
# cleaned.pval 6.063518e-11 1.0000000 1.0000000
# cleaned.fdr  3.577957e-09 0.6446391 0.5644553
# cleaned.fwer 1.453986e-01 0.7534788 1.0000000
#
# $HIPPO
#                   exon       jxn           tx
# expr.pval    0.3463280 0.3463280 0.000000e+00
# expr.fdr     0.9248146 0.9248146 5.099046e-26
# expr.fwer    1.0000000 1.0000000 0.000000e+00
# cleaned.pval 0.6692196 0.6692196 0.000000e+00
# cleaned.fdr  0.9751166 0.9751166 0.000000e+00
# cleaned.fwer 1.0000000 1.0000000 0.000000e+00

lapply(corr_de_feat, function(region) { do.call(cbind, lapply(region, calc_chi_pval)) < 0.05 } )
# $DLPFC
#               exon   jxn    tx
# expr.pval    FALSE FALSE FALSE
# expr.fdr     FALSE FALSE FALSE
# expr.fwer     TRUE FALSE FALSE
# cleaned.pval  TRUE FALSE FALSE
# cleaned.fdr   TRUE FALSE FALSE
# cleaned.fwer FALSE FALSE FALSE
#
# $HIPPO
#               exon   jxn   tx
# expr.pval    FALSE FALSE TRUE
# expr.fdr     FALSE FALSE TRUE
# expr.fwer    FALSE FALSE TRUE
# cleaned.pval FALSE FALSE TRUE
# cleaned.fdr  FALSE FALSE TRUE
# cleaned.fwer FALSE FALSE TRUE

save(corr_de_gene, corr_de_feat, file = 'rda/corr_de.Rdata')



## Plot t-stat DE vs correlation
pdf('pdf/corr_vs_tstat_fwer.pdf', useDingbats = FALSE)
mapply(function(out, set) {
    plot(out$t, out$corr.expr, col = alpha(ifelse(out$expr.fwer < 0.05, 'red', 'black'), 1/10), pch = 19, xlab = 'SCZD t-statistic', ylab = 'Corr - expr', main = set)
    plot(out$t, out$corr.cleaned, col = alpha(ifelse(out$cleaned.fwer < 0.05, 'red', 'black'), 1/10), pch = 19, xlab = 'SCZD t-statistic', ylab = 'Corr - cleaned', main = set)
}, outGene, names(outGene))
dev.off()

pdf('pdf/corrpvals_vs_sczd_fwer.pdf', useDingbats = FALSE)
mapply(function(out, set) {
    plot(-log10(out$adj.P.Val), -log10(out$expr.fwer), col = alpha(ifelse(out$expr.fwer < 0.05, 'red', 'black'), 1/10), pch = 19, xlab = 'SCZD FDR', ylab = 'Corr FWER - expr', main = set)
    plot(-log10(out$adj.P.Val), -log10(out$cleaned.fwer), col = alpha(ifelse(out$cleaned.fwer < 0.05, 'red', 'black'), 1/10), pch = 19, xlab = 'SCZD FDR', ylab = 'Corr FWER - cleaned', main = set)
}, outGene, names(outGene))
dev.off()






run_go <- function(genes, ont = c('BP', 'MF', 'CC')) {
    ## Change to ENSEMBL ids and remove NAs
    genes_ens <- lapply(lapply(genes, function(x) { gsub('\\..*', '', x) }), function(y) y[!is.na(y)])

    #genes_venn <- venn(genes_ens, show.plot = FALSE)

    ## Run GO analysis
    go_cluster <- lapply(ont, function(bp) {
        message(paste(Sys.time(), 'running GO analysis for', bp))
        tryCatch(compareCluster(genes_ens, fun = "enrichGO",
            universe = uni, OrgDb = 'org.Hs.eg.db',
            ont = bp, pAdjustMethod = "BH",
            pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
            readable = TRUE, keyType = 'ENSEMBL'),
            error = function(e) { return(NULL) })
    })
    names(go_cluster) <- ont
    
    message(paste(Sys.time(), 'running GO analysis for KEGG'))
    genes_ncbi <- lapply(lapply(genes_ens, bitr, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db'), function(x) x$ENTREZID)
    
    uni_ncbi <- bitr(uni, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID
    
    go_cluster$KEGG <- tryCatch(compareCluster(genes_ncbi, fun = 'enrichKEGG',
        universe = uni_ncbi, organism = 'hsa', pAdjustMethod = 'BH',
        pvalueCutoff = 0.1, qvalueCutoff = 0.05, keyType = 'ncbi-geneid'),
        error = function(e) { return(NULL) })
    
    return(go_cluster)
}


corr_genes <- lapply(c('expr', 'cleaned'), function(set) {
    fwer <- gene_pinfo[, paste0(set, '.fwer')]
    rownames(gene_pinfo[fwer < 0.05, ])
})
names(corr_genes) <- c('expr', 'cleaned')

## Use all genes as the universe
uni <- gsub('\\..*', '', rownames(gene_pinfo))
length(uni)
# [1] 24652
sapply(corr_genes, length)
# expr cleaned
#  340     930

system.time( go_corr_genes <- run_go(corr_genes) )
message(paste(Sys.time(), 'saving rda/go_corr_genes.Rdata'))
save(go_corr_genes, file = 'rda/go_corr_genes.Rdata')


# simplify_go <- function(x) {
#     #gsub('QSV|IPPO|LPFC', '', x)
#     #gsub('_matchQSV', '', x)
#     gsub('IPPO|LPFC|ontrol|chizo|CZD|xon_|xn_|x_|ene_|\\.', '', x)
# }

plot_go <- function(go_cluster, cat = 10) {
    lapply(names(go_cluster), function(bp) {
        go <- go_cluster[[bp]]
        if(is.null(go)) {
            message(paste(Sys.time(), 'found no results for', bp))
            return(NULL)
        }

        ## Simplify names
        # go@compareClusterResult$Cluster <- simplify_go(go@compareClusterResult$Cluster)
        # names(go@geneClusters) <- simplify_go(names(go@geneClusters))

        print(plot(go, title = paste('ontology:', bp), font.size = 18, showCategory = cat, includeAll = TRUE))
        return(NULL)
    })
}

## Plot using conda_R/3.4.x
pdf('pdf/go_corr_genes.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_go(go_corr_genes)
dev.off()

pdf('pdf/go_all_corr_genes.pdf', width = 16, height = 70, useDingbats = FALSE)
plot_go(go_corr_genes, cat = NULL)
dev.off()




# n_de_sign <- do.call(rbind, lapply(c(0.05, 0.1, 0.15, 0.2), function(cut) {
#     xx <- lapply(c(outGene, do.call(c, outFeat)), function(x) {
#         y <- table(factor(x$adj.P.Val < cut, levels = c('FALSE', 'TRUE')), factor(sign(x$logFC), levels = c(-1, 0, 1)))
#         data.frame(de_status = rownames(y), n = as.vector(y), sign = rep(colnames(y), each = 2), cutoff = cut)
#     })
#     for(i in seq_len(length(xx))) { xx[[i]]$model = names(xx)[i] }
#     names(xx) <- NULL
#     do.call(rbind, xx)
# }))
# n_de_sign$group <- ifelse(n_de_sign$sign == 0, 'none', ifelse(n_de_sign$sign == -1, 'Control', 'SCZD'))
# n_de_sign


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
#  date     2018-08-21
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  acepack                1.4.1     2016-10-29 CRAN (R 3.5.0)
#  AnnotationDbi        * 1.42.1    2018-05-17 Bioconductor
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
#  mime                   0.5       2016-07-07 CRAN (R 3.5.0)
#  munsell                0.5.0     2018-06-12 CRAN (R 3.5.0)
#  nnet                   7.3-12    2016-02-02 CRAN (R 3.5.0)
#  org.Hs.eg.db         * 3.6.0     2018-05-03 Bioconductor
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
