library('SummarizedExperiment')
library('jaffelab')
library('devtools')

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata", verbose = TRUE)


load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda", verbose = TRUE)
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

# ## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
# snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
# snpMap = snpMap[-snpInd,]
# snp = snp[-snpInd,]

# ## drop SNPs not mapping to hg38
# keepIndex = which(!is.na(snpMap$chr_hg38))
# snpMap = snpMap[keepIndex,]
# snp = snp[keepIndex,]

pd = colData(rse_gene)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum

dim(mds)
dim(pd)
dim(snp)


eQTLfiles <- c(
    'HIPPOfull' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_hippo_4features.rda',
    'DLPFCfull' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda',
    'INTfull' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_hippo_4features.rda',
    'HIPPOraggr' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/eqtl_tables/mergedEqtl_output_hippo_raggr_4features.rda',
    'DLPFCraggr' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/eqtl_tables/mergedEqtl_output_dlpfc_raggr_4features.rda')
    
stopifnot(all(file.exists(eQTLfiles)))



cleaned <- lapply(c('HIPPO', 'DLPFC', 'INT'), function(modtype) {
    if(modtype == 'HIPPO') {
        keepInd = which(colData(rse_gene)$Age > 13 & colData(rse_gene)$Region == "HIPPO")
        load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/rdas/pcs_hippo_4features_filtered_over13.rda', verbose = TRUE)
    } else if (modtype == 'DLPFC') {
        keepInd = which(colData(rse_gene)$Age > 13 & colData(rse_gene)$Region == "DLPFC")
        load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/rdas/pcs_dlpfc_4features_filtered_over13.rda', verbose = TRUE)
    } else if (modtype == 'INT') {
        keepInd = which(colData(rse_gene)$Age > 13)
        load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/rdas/pcs_4features_combined_regions_filtered_over13.rda', verbose = TRUE)
    }

    ## extract pd and rpkms
    pd = colData(rse_gene[,keepInd])
    
    exprs <- list(
        'gene' = assays(rse_gene)$rpkm[,keepInd],
        'exon' = assays(rse_exon)$rpkm[,keepInd],
        'jxn' = assays(rse_jxn)$rp10m[,keepInd],
        'tx' = assays(rse_tx)$tpm[,keepInd]
    )
    
    pd$Dx = factor(pd$Dx, levels = c("Control", "SCZD"))

    if(modtype == 'INT') {
        mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]) + Region,
        	data = pd)
    } else {
        mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]),
        	data = pd)
    }
    colnames(mod)[4:8] = colnames(mds)[1:5]
    
    mods <- lapply(list('gene' = genePCs, 'exon' = exonPCs, 'jxn' = jxnPCs, 'tx' = txPCs), function(pc) {
        cbind(mod, pc)
    })
    
    res <- mapply(function(expr, model) {
        cleaningY(expr, model, P = 1)
    }, exprs, mods)
    
    return(res)
})


eqtls <- lapply(eQTLfiles, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(allEqtl)
})

#################################
### export database files #######

e_stats <- function(allEqtls, name) {
    eqtlStatsOut = allEqtls[,c(1:7,18)]
    eqtlStatsOut = as.data.frame(eqtlStatsOut)
    write.table(eqtlStatsOut, row.names=FALSE, sep = "\t",
    	file=gzfile(paste0("eqtl_stats_", name, ".txt.gz")))
}
mapply(e_stats, eqtls, names(eqtls))



## expression anno
e_ann <- function(allEqtls, name) {
    exprsAnnoOut = allEqtls[,c(2,9:11,13:17)]
    exprsAnnoOut = exprsAnnoOut[!duplicated(exprsAnnoOut$Feature),]
    exprsAnnoOut$chr = ss(exprsAnnoOut$Coordinates, ":")
    exprsAnnoOut$start = as.integer(ss(ss(exprsAnnoOut$Coordinates, ":",2),"-"))
    exprsAnnoOut$end = as.integer(ss(ss(ss(exprsAnnoOut$Coordinates, ":",2),"-",2),"\\("))
    exprsAnnoOut = as.data.frame(exprsAnnoOut)
    exprsAnnoOut$WhichTx = sapply(exprsAnnoOut$WhichTx, paste, collapse=";")
    write.table(exprsAnnoOut, row.names=FALSE, sep = "\t",
    	file=gzfile(paste0("eqtl_exprsAnno_", name, ".txt.gz")))    
}
mapply(e_ann, eqtls, names(eqtls))

## SNP anno	
snpAnnoOut = snpMap
colnames(snpAnnoOut)[c(1,3,4:5)] = c("chr","pos","Counted","Ref")
snpAnnoOut$chr = paste0("chr", snpAnnoOut$chr)
snpAnnoOut = snpAnnoOut[,c(2,6,1,3:5,7)]

snpAnnoOut = as.data.frame(snpAnnoOut)
write.table(snpAnnoOut, row.names=FALSE, sep = "\t",
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/AllEqtls/snpAnno_allConsidered.txt.gz"))

## expression data
cleanExprs = rbind(cleanGeneSub, cleanExonSub, cleanJxnSub,	cleanTxSub, cleanErSub)
write.table(cleanExprs, sep = "\t",
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/AllEqtls/eqtl_log2cleanExprsData_allStats.txt.gz"))

## snp data
write.table(snpSub, sep = "\t",
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/AllEqtls/eqtl_snpData_allStats.txt.gz"))

## replication stats
repStatsOut = allEqtls[,c(1:2, 19:34)]
write.table(repStatsOut, sep = "\t",row.names=FALSE,
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/AllEqtls/replication_stats_allStats.txt.gz"))
## Reproducibility information


print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
