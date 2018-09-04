library('SummarizedExperiment')
library('xlsx')
library('devtools')

## For Table S3
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/correlation/rda/indv_kegg_sczd.Rdata', verbose = TRUE)
info <- subset(indv_kegg_cleaned_sczd, padj < 0.05)
colnames(info)[colnames(info) == 'padj'] <- 'FDR'
write.xlsx(info, file = 'SupplementaryTable3.xlsx', sheetName = 'KEGG', append = FALSE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/correlation/rda/indv_go_sczd.Rdata', verbose = TRUE)
info <- lapply(lapply(indv_go_cleaned_sczd, subset, padj < 0.05), function(x) { x[, -c(11, 12)]})
info <- lapply(info, function(x) { colnames(x)[colnames(x) == 'padj'] <- 'FDR'; return(x) })
for(i in 1:3) {
    write.xlsx(info[[i]], file = 'SupplementaryTable3.xlsx', sheetName = paste0('GO-', names(info)[i]), append = TRUE)
}

## For table S4
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/degradation_rse_phase2_usingJoint_justFirst.rda", verbose = TRUE)
info <- as.data.frame(rowRanges(cov_rse))[, -c(6:7)]
rownames(info) <- NULL
write.xlsx(info, file = 'SupplementaryTable4.xlsx', sheetName = 'Top1kDegradationExpressedRegions', append = FALSE)

## For table S2
## Development model first
source('/dcl01/lieber/ajaffe/lab/brainseq_phase2/development/load_funs.R')
lvls <- c('gene', 'exon', 'jxn', 'tx')
rses <- lapply(lvls, load_foo)
names(rses) <- lvls

load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/development/rda/pcheck_both.Rdata', verbose = TRUE)
stopifnot(sum(sapply(rses, nrow)) == nrow(pcheck_both))
info <- split(pcheck_both[, -which(colnames(pcheck_both) %in% c('global_fdr', 'global_bonf'))], pcheck_both$type)[names(rses)]
stopifnot(all(sapply(rses, nrow) - sapply(info, nrow) == 0))
for(i in 1:4) {
    info[[i]] <- cbind(info[[i]], as.data.frame(rowRanges(rses[[i]])))
    ## Replication: Bonf < 1% and BrainSpan p-val < 5%
    info[[i]]$replicates_in_BrainSpan <- info[[i]]$P.Bonf < 0.01 & info[[i]]$span_P.Value < 0.05
}
sapply(info, function(x) sum(x$replicates_in_BrainSpan) )
#  gene   exon    jxn     tx
# 10839 169253 143895   1715

## Save for grouping across the 3 models
deres <- vector('list', 3)
names(deres) <- c('development', 'region', 'sczd')
deres$development <- info

## Move on to the region-specific (adult, prenatal) results
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/region_specific/rda/pcheck_both.Rdata', verbose = TRUE)
stopifnot(sum(sapply(rses, nrow)) * 2 == nrow(pcheck_both))
info <- split(pcheck_both[, -which(colnames(pcheck_both) %in% c('global_fdr', 'global_bonf'))], pcheck_both$type)[names(rses)]
stopifnot(all(sapply(rses, nrow) * 2 - sapply(info, nrow) == 0))
for(i in 1:4) {
    info[[i]] <- cbind(info[[i]], rbind(as.data.frame(rowRanges(rses[[i]])), as.data.frame(rowRanges(rses[[i]]))))
    ## Replication: log2FC has to have the same sing, Bonf < 1% and BrainSpan p-val < 5%
    info[[i]]$replicates_in_BrainSpan <- sign(info[[i]]$logFC) == sign(info[[i]]$span_logFC) & info[[i]]$P.Bonf < 0.01 & info[[i]]$span_P.Value < 0.05
    info[[i]]$age[info[[i]]$age == 'fetal'] <- 'prenatal'
    info[[i]]$span_age[info[[i]]$span_age == 'fetal'] <- 'prenatal'
    rownames(info[[i]]) <- gsub('fetal', 'prenatal', rownames(info[[i]]))
}
sapply(info, function(x) sum(x$replicates_in_BrainSpan[x$age == 'adult']) )
# gene  exon   jxn    tx
#  1612 15442  5561  1739
sapply(info, function(x) sum(x$replicates_in_BrainSpan[x$age == 'prenatal']) )
# gene exon  jxn   tx
#   32   71   18    3
deres$region <- info

## Now the SCZD case-control results
outFeat <- lapply(c('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda', '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda'), function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    outTx$ensemblID <- gsub('\\..*', '', outTx$gene_id)
    return(list('gene' = outGene, 'exon' = outExon, 'jxn' = outJxn, 'tx' = outTx))
})
names(outFeat) <- c('DLPFC', 'HIPPO')
info <- mapply(rbind, outFeat[[1]], outFeat[[2]])

for(i in 1:4) {
    info[[i]]$type <- names(info)[i]
}
deres$sczd <- info

## Save full info
save(deres, file = 'deres.Rdata')

## Subset to DE features
rep_in_span <- function(x) {
    x[x$replicates_in_BrainSpan, ]
}
deres_sig <- deres
deres_sig[[1]] <- lapply(deres_sig[[1]], rep_in_span)
deres_sig[[2]] <- lapply(deres_sig[[2]], rep_in_span)
deres_sig[[3]] <- lapply(deres_sig[[3]], function(x) { x[x$adj.P.Val < 0.1, ] })

sapply(deres_sig, lapply, nrow)
#      development region sczd
# gene 10839       1644   733
# exon 169253      15513  2384
# jxn  143895      5579   155
# tx   1715        1742   21
save(deres_sig, file = 'deres_sig.Rdata')

## Export to Excel TableS2
# load('deres_sig.Rdata', verbose = TRUE)
# system('rm SupplementaryTable2.xlsx')
for(i in 1:3) {
    for(j in 1:4) {
        message(paste(Sys.time(), 'processing', names(deres_sig)[i], 'at the', names(deres_sig[[i]])[j], 'level'))
        write.csv(deres_sig[[i]][[j]], file = paste0('SupplementaryTable2_' names(deres_sig)[i], '_', names(deres_sig[[i]])[j], '.csv'))
    }
}



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()