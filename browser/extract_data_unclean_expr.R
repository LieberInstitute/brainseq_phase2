## Adapted from /users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/filter_replicated_eqtls.R

library('SummarizedExperiment')
library("data.table")
library("sessioninfo")

## Load expr data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata", verbose = TRUE)

## Fix exon labels
load('rda/exon_name_map.Rdata', verbose = TRUE)

## Replace names
rownames(rse_exon) <- exon_name_map$gencode
rowRanges(rse_exon)$exon_libdID <- exon_name_map$libd_bsp2
rowRanges(rse_exon)$exon_gencodeID <- exon_name_map$gencode
rowRanges(rse_exon)$exon_libdID_gtex <- exon_name_map$libd_gtex


## Load pheno data
load("rda/pd.Rdata", verbose = TRUE)


## Cleaned eQTL expression
cleaned <- lapply(c('hippo', 'dlpfc', 'interaction'), function(modtype) {
    message(paste(Sys.time(), 'processing', modtype))
    if(modtype == 'hippo') {
        load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/rdas/pcs_hippo_4features_filtered_over13.rda', verbose = TRUE)
    } else if (modtype == 'dlpfc') {
        load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/rdas/pcs_dlpfc_4features_filtered_over13.rda', verbose = TRUE)
    } else if (modtype == 'interaction') {
        load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/rdas/pcs_4features_combined_regions_filtered_over13.rda', verbose = TRUE)
    }
    
    keepInd <- pd[[paste0('analysis_eqtl_', modtype)]]

    ## extract pd and rpkms
    exprs <- list(
        'gene' = assays(rse_gene)$rpkm[,keepInd],
        'exon' = assays(rse_exon)$rpkm[,keepInd],
        'jxn' = assays(rse_jxn)$rp10m[,keepInd],
        'tx' = assays(rse_tx)$tpm[,keepInd]
    )
    exprs <- lapply(exprs, function(x) { log2(x + 1) })
    return(exprs)
})
names(cleaned) <- c('hippo', 'dlpfc', 'interaction')

## Export
mapply(function(exprs, type) {
    message(paste(Sys.time(), 'processing', type))
    mapply(function(expr, exprtype) {
        message(paste(Sys.time(), 'processing', exprtype))
        fwrite(as.data.frame(expr), sep = '\t', row.names = TRUE, file = paste0('BrainSeqPhaseII_unclean_expression_eqtl_', type, '_', exprtype, '.txt'))
    }, exprs, names(exprs))
    return(NULL)
}, cleaned, names(cleaned))
cleaned_eqtl <- cleaned
rm(cleaned_eqtl, cleaned)

## Cleaned development expression
de_analyses <- c('development', 'sczd_casecontrol_interaction', 'sczd_casecontrol_hippo', 'sczd_casecontrol_dlpfc', 'regionspecific_adult', 'regionspecific_prenatal')
cleaned <- lapply(de_analyses, function(modtype) {
    message(paste(Sys.time(), 'processing', modtype))
    
    if(modtype == 'development') {
        keepInd <- pd$analysis_development
    } else if(modtype == 'sczd_casecontrol_interaction') {
        load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold.Rdata', verbose = TRUE)
        keepInd <- keepIndex
    } else if(modtype == 'sczd_casecontrol_dlpfc') {
        load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_DLPFC.Rdata', verbose = TRUE)
        keepInd <- keepIndex
    } else if(modtype == 'sczd_casecontrol_hippo') {
        load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata', verbose = TRUE)
        keepInd <- keepIndex
    } else if(modtype == 'regionspecific_adult') {
        keepInd <- pd$analysis_regionspecific_adult
    } else if(modtype == 'regionspecific_prenatal') {
        keepInd <- pd$analysis_regionspecific_prenatal
    }

    ## extract pd and rpkms
    exprs <- list(
        'gene' = assays(rse_gene)$rpkm[,keepInd],
        'exon' = assays(rse_exon)$rpkm[,keepInd],
        'jxn' = assays(rse_jxn)$rp10m[,keepInd],
        'tx' = assays(rse_tx)$tpm[,keepInd]
    )
    exprs <- lapply(exprs, function(x) { log2(x + 1) })
    
    return(exprs)
})
names(cleaned) <- de_analyses

## Export
mapply(function(exprs, type) {
    message(paste(Sys.time(), 'processing', type))
    mapply(function(expr, exprtype) {
        message(paste(Sys.time(), 'processing', exprtype))
        fwrite(as.data.frame(expr), sep = '\t', row.names = TRUE, file = paste0('BrainSeqPhaseII_unclean_expression_', type, '_', exprtype, '.txt'))
    }, exprs, names(exprs))
    return(NULL)
}, cleaned, names(cleaned))
cleaned_de_analyses <- cleaned
rm(cleaned_de_analyses, cleaned)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
