## Usage:
# qrsh -l mem_free=100G,h_vmem=100G,h_fsize=100G
# Rscript corr_analysis.R  > logs/corr_analysis.txt 2>&1

library('SummarizedExperiment')
library('jaffelab')
library('devtools')

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

dir.create('rda', showWarnings = FALSE)
message(paste(Sys.time(), 'saving rse and modQsva info'))
save(simple_rse, modQsva, file = 'rda/rse_and_modQsva.Rdata')
message(paste(Sys.time(), 'saving expression info'))
save(expr, cleaned, file = 'rda/expr_and_cleaned.Rdata')
message(paste(Sys.time(), 'saving corr results info'))
save(corr_expr, corr_cleaned, file = 'rda/corrs.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
