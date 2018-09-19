## Adapted from /users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/filter_replicated_eqtls.R

library('SummarizedExperiment')
library('data.table')
library('jaffelab')
library('devtools')
library('readxl')

dir.create('rda', showWarnings = FALSE)

## Load expr data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata", verbose = TRUE)

## Fix exon labels
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_both/rse_gtex_exon.Rdata', verbose = TRUE)
ov <- findOverlaps(rowRanges(rse_exon), rowRanges(rse_gtex_exon), type = 'equal', ignore.strand = FALSE)
length(unique(queryHits(ov)))
length(unique(subjectHits(ov)))
exon_name_map <- data.frame(
    'libd_bsp2' = rownames(rse_exon),
    'gencode' = NA,
    'libd_gtex' = NA,
    stringsAsFactors = FALSE
)
exon_name_map$gencode[queryHits(ov)] <- rowRanges(rse_gtex_exon)$exon_gencodeID[subjectHits(ov)]
exon_name_map$libd_gtex[queryHits(ov)] <- rownames(rse_gtex_exon)[subjectHits(ov)]

## libd_bsp2 and libd_gtex are not the same!
with(exon_name_map, identical(libd_bsp2[!is.na(libd_gtex)], libd_gtex[!is.na(libd_gtex)]))
# [1] FALSE

## 2 are missing at this point
subset(exon_name_map, is.na(gencode))
#        libd_bsp2 gencode libd_gtex
# 171891   e514496    <NA>      <NA>
# 290930   e862727    <NA>      <NA>
options(width = 200)
rowRanges(rse_exon[!seq_len(nrow(rse_exon)) %in% queryHits(ov)])
# GRanges object with 2 ranges and 11 metadata columns:
#           seqnames              ranges strand |    Length         gencodeID       ensemblID      gene_type      Symbol  EntrezID       Class         meanExprs     NumTx         gencodeTx passExprsCut
#              <Rle>           <IRanges>  <Rle> | <integer>       <character>     <character>    <character> <character> <integer> <character>         <numeric> <integer>   <CharacterList>    <logical>
#   e514496     chr8 123348034-123348130      - |        97 ENSG00000283172.1 ENSG00000283172          miRNA                  <NA>       InGen 0.453936954830773         1 ENST00000636914.1         TRUE
#   e862727    chr16     3384459-3384941      + |       483 ENSG00000262621.4 ENSG00000262621 protein_coding                  <NA>       InGen 0.301426795813603         1 ENST00000618352.1         TRUE
#   -------
#   seqinfo: 25 sequences from an unspecified genome; no seqlengths

## Resolve manually since both genes are 1 transcript genes and the coordinates match only 1 exon
# http://useast.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000283172;r=8:123348034-123348130;t=ENST00000636914
exon_name_map$gencode[exon_name_map$libd_bsp2 == 'e514496'] <- 'ENSE00003792174'
# http://useast.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000262621;r=16:3382113-3397745;t=ENST00000618352
exon_name_map$gencode[exon_name_map$libd_bsp2 == 'e862727'] <- 'ENSE00003747793'

## Replace names
rownames(rse_exon) <- exon_name_map$gencode
rowRanges(rse_exon)$exon_libdID <- exon_name_map$libd_bsp2
rowRanges(rse_exon)$exon_gencodeID <- exon_name_map$gencode
rowRanges(rse_exon)$exon_libdID_gtex <- exon_name_map$libd_gtex

## Export
exon_name_map <- data.table(exon_name_map)
save(exon_name_map, file = 'rda/exon_name_map.Rdata')
## Don't really need it as a separate file, since it's part of the rse_exon (see below)
fwrite(exon_name_map, row.names = FALSE, sep = '\t', file = 'rda/BrainSeqPhaseII_exon_name_map.txt')

## Export expr annotation
export_ann <- function(x, f) {
    y <- as.data.frame(rowRanges(x))
    y$feature_id <- rownames(y)
    fwrite(y, row.names = FALSE, sep = '\t', file = f)
}
export_ann(rse_gene, 'BrainSeqPhaseII_feature_annotation_gene.txt')
export_ann(rse_exon, 'BrainSeqPhaseII_feature_annotation_exon.txt')
export_ann(rse_jxn, 'BrainSeqPhaseII_feature_annotation_jxn.txt')
export_ann(rse_tx, 'BrainSeqPhaseII_feature_annotation_tx.txt')

## Load pheno data
pd = colData(rse_gene)

## Add correlation results
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/correlation/rda/indv_corr.Rdata', verbose = TRUE)
rm(indv_expr)
indv_cleaned <- do.call(cbind, indv_cleaned)
colnames(indv_cleaned) <- paste0('DLFPC_HIPPO_correlation_', colnames(indv_cleaned))
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/correlation/rda/expr_and_cleaned.Rdata', verbose = TRUE)
rm(expr)
rownames(indv_cleaned) <- colnames(cleaned$DLPFC$geneRpkm)
tmp <- matrix(NA, nrow = nrow(pd), ncol = ncol(indv_cleaned))
rownames(tmp) <- rownames(pd)
m <- match(rownames(tmp), rownames(indv_cleaned))
tmp[!is.na(m), ] <- indv_cleaned[m[!is.na(m)], ]
colnames(tmp) <- colnames(indv_cleaned)
pd <- cbind(pd, tmp)
## The correlation values get all assigned to DLPFC
# > table(pd$Region[!is.na(pd$DLFPC_HIPPO_correlation_geneRpkm)])
#
# DLPFC
#   265
rm(tmp, cleaned, indv_cleaned, m)

## Fix Dx labels
pd$Dx <- factor(ifelse(pd$Dx == 'Control', 'Control', 'SCZD'))

## Relevel some variables
pd$Region <- relevel(factor(pd$Region), 'DLPFC')
pd$Race <- relevel(factor(pd$Race), ref = 'CAUC')
pd$Sex <- relevel(factor(pd$Sex), ref = 'F')

## Add means
pd$mean_mitoRate <- mean(pd$mitoRate)
pd$mean_totalAssignedGene <- mean(pd$totalAssignedGene)
pd$mean_rRNA_rate <- mean(pd$rRNA_rate)
pd$mean_RIN <- mean(pd$RIN)

## Add age linear splines
fetal <- ifelse(pd$Age < 0, 1,0)
birth <- pd$Age
birth[birth < 0] <- 0 # linear spline
infant <- pd$Age - 1
infant[infant < 0] <- 0 # linear spline
child <- pd$Age - 10
child[child < 0] <- 0 # linear spline
teen <- pd$Age - 20
teen[teen < 0] <- 0 # linear spline
adult <- pd$Age - 50
adult[adult < 0] <- 0 # linear spline
pd$agespline_fetal <- fetal
pd$agespline_birth <- birth
pd$agespline_infant <- infant
pd$agespline_child <- child
pd$agespline_teen <- teen
pd$agespline_adult <- adult
rm(fetal, birth, infant, child, teen, adult)

## Add prep protocol
pd$protocol <- 'RiboZeroGold'
hipxl <- read_excel('/dcl01/lieber/ajaffe/lab/brainseq_phase2/LIBD_PhaseII_HIPPO_RiboZero_sample_list_01_28_2015.xlsx')
hipHMR <- as.integer(gsub('R', '', rownames(pd))) %in% hipxl$RNum[hipxl$Protocol == 'RiboZeroHMR']
pd$protocol[hipHMR] <- 'RiboZeroHMR'
rm(hipHMR, hipxl)

## Add labels for each analysis
pd$analysis_regionspecific_adult <- with(pd, Dx == 'Control' & Age >= 18)
pd$analysis_regionspecific_prenatal <- with(pd, Dx == 'Control' & Age <= 0)
pd$analysis_development <- with(pd, Dx == 'Control')
pd$analysis_sczd_casecontrol_dlpfc <- with(pd, Age > 17 & Region == 'DLPFC')
pd$analysis_sczd_casecontrol_hippo <- with(pd, Age > 17 & protocol == 'RiboZeroHMR' & Region == 'HIPPO')
pd$analysis_sczd_casecontrol_interaction <- with(pd, Age > 17 & (protocol == 'RiboZeroHMR' | Region == 'DLPFC'))
pd$analysis_eqtl_dlpfc <- with(pd, Age > 13 & Region == 'DLPFC')
pd$analysis_eqtl_hippo <- with(pd, Age > 13 & Region == 'HIPPO')
pd$analysis_eqtl_interaction <- with(pd, Age > 13)

## Print the sample sizes for each analysis
t(summary(as.data.frame(pd[, grep('analysis_', colnames(pd))])))[, -1]
# analysis_regionspecific_adult         FALSE:440       TRUE :460
# analysis_regionspecific_prenatal      FALSE:844       TRUE :56
# analysis_development                  FALSE:286       TRUE :614
# analysis_sczd_casecontrol_dlpfc       FALSE:521       TRUE :379
# analysis_sczd_casecontrol_hippo       FALSE:567       TRUE :333
# analysis_sczd_casecontrol_interaction FALSE:188       TRUE :712
# analysis_eqtl_dlpfc                   FALSE:503       TRUE :397
# analysis_eqtl_hippo                   FALSE:505       TRUE :395
# analysis_eqtl_interaction             FALSE:108       TRUE :792

## Export
dim(pd)
# [1] 900  78
write.table(pd, sep = "\t", row.names=FALSE, file = 'BrainSeqPhaseII_sample_metadata.txt', quote = FALSE)
system('wc -l BrainSeqPhaseII_sample_metadata.txt')
## Save in R format too, since reading from the text would be a bit of a pain (due to NumericList variables)
save(pd, file = 'rda/pd.Rdata')


load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda", verbose = TRUE)
## No need for the actual snp info or mds info (first 10 snpPCs are already on pd)
rm(snp, mds)

## No need to change from 23 to X
# > table(snpMap$CHR == "23")
#
#   FALSE
# 7023860

## Don't need this since the 2 columns are already there
# snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

## No need to drop this once, since it won't show up later in the stats
# ## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
# snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
# snpMap = snpMap[-snpInd,]
# snp = snp[-snpInd,]

# ## drop SNPs not mapping to hg38
# keepIndex = which(!is.na(snpMap$chr_hg38))
# snpMap = snpMap[keepIndex,]
# snp = snp[keepIndex,]

## SNP anno	
snpAnnoOut = snpMap
colnames(snpAnnoOut) <- tolower(colnames(snpAnnoOut))
snpAnnoOut$chr = paste0("chr", snpAnnoOut$chr)
colnames(snpAnnoOut)[colnames(snpAnnoOut) == 'chr'] <- 'chr_hg19'
colnames(snpAnnoOut)[colnames(snpAnnoOut) == 'pos'] <- 'pos_hg19'
## Just export everything (the column indexes don't match this table anyway)
# snpAnnoOut = snpAnnoOut[,c(2,6,1,3:5,7)]
#
# snpAnnoOut = as.data.frame(snpAnnoOut)

## Re-order a bit: keep the snp as the first column
snpAnnoOut <- snpAnnoOut[, c('snp', 'chr_hg38', 'pos_hg38', 'chr_hg19', 'pos_hg19', 'cm', 'counted', 'alt', 'type', 'newref', 'newcount', 'name', 'rsnumguess')]

## Export
dim(snpAnnoOut)
# [1] 7023860      13
fwrite(snpAnnoOut, row.names=FALSE, sep = "\t",	file = "BrainSeqPhaseII_snp_annotation.txt")
system('wc -l BrainSeqPhaseII_snp_annotation.txt')
system('head BrainSeqPhaseII_snp_annotation.txt')
rm(snpAnnoOut, snpMap)

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
    pd2 = pd[keepInd, ]
    
    exprs <- list(
        'gene' = assays(rse_gene)$rpkm[,keepInd],
        'exon' = assays(rse_exon)$rpkm[,keepInd],
        'jxn' = assays(rse_jxn)$rp10m[,keepInd],
        'tx' = assays(rse_tx)$tpm[,keepInd]
    )

    if(modtype == 'interaction') {
        mod = model.matrix(~Dx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + Region,
        	data = pd2)
    } else {
        mod = model.matrix(~Dx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5,
        	data = pd2)
    }
    
    mods <- lapply(list('gene' = genePCs, 'exon' = exonPCs, 'jxn' = jxnPCs, 'tx' = txPCs), function(pc) {
        cbind(mod, pc)
    })
    
    res <- mapply(function(expr, model) {
        as.data.frame(cleaningY(expr, model, P = 1))
    }, exprs, mods)
    
    return(res)
})
names(cleaned) <- c('hippo', 'dlpfc', 'interaction')

## Export
mapply(function(exprs, type) {
    message(paste(Sys.time(), 'processing', type))
    mapply(function(expr, exprtype) {
        message(paste(Sys.time(), 'processing', exprtype))
        fwrite(expr, sep = '\t', row.names = FALSE, file = paste0('BrainSeqPhaseII_clean_expression_eqtl_', type, '_', exprtype, '.txt'))
    }, exprs, names(exprs))
    return(NULL)
}, cleaned, names(cleaned))
cleaned_eqtl <- cleaned
message(paste(Sys.time(), 'saving clean_expr_eqtl.Rdata'))
save(cleaned_eqtl, file = 'rda/clean_expr_eqtl.Rdata')
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
    pd2 = pd[keepInd, ]
    
    exprs <- list(
        'gene' = assays(rse_gene)$rpkm[,keepInd],
        'exon' = assays(rse_exon)$rpkm[,keepInd],
        'jxn' = assays(rse_jxn)$rp10m[,keepInd],
        'tx' = assays(rse_tx)$tpm[,keepInd]
    )
    
    if(modtype == 'development') {
        design = model.matrix(~Age *Region + agespline_fetal * Region + agespline_birth *Region + agespline_infant *Region + agespline_child * Region + agespline_teen * Region + agespline_adult * Region + Sex + Region + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN,
        	data = pd2)
            
        ## Re-organize to protect variables of interest
        protect <- grepl(':RegionHIPPO|RegionHIPPO:|Intercept', colnames(design))
        design <- cbind(design[, protect], design[, !protect])
        design_p <- sum(protect)
    } else if(modtype == 'sczd_casecontrol_interaction') {
        design <- cbind(with(pd2, model.matrix(~ Dx * Region)),
            modQsva[, -grep('Dx|RegionHIPPO|Intercept', colnames(modQsva))]
        )
        design_p <- 4
    } else if (modtype %in% c('sczd_casecontrol_dlpfc', 'sczd_casecontrol_hippo')) {
        design <- modQsva
        design_p <- 2
    } else if (grepl('regionspecific', modtype)) {
        design <- model.matrix(~Region + Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN, data = pd2)
        design_p <- 2
    }
    
    res <- lapply(exprs, function(expr) { as.data.frame(cleaningY(expr, design, P = design_p)) })
    return(res)
})
names(cleaned) <- de_analyses

## Export
mapply(function(exprs, type) {
    message(paste(Sys.time(), 'processing', type))
    mapply(function(expr, exprtype) {
        message(paste(Sys.time(), 'processing', exprtype))
        fwrite(expr, sep = '\t', row.names = FALSE, file = paste0('BrainSeqPhaseII_clean_expression_', type, '_', exprtype, '.txt'))
    }, exprs, names(exprs))
    return(NULL)
}, cleaned, names(cleaned))
cleaned_de_analyses <- cleaned
message(paste(Sys.time(), 'saving clean_expr_de_analysis.Rdata'))
save(cleaned_de_analyses, file = 'rda/clean_expr_de_analyses.Rdata')
rm(cleaned_de_analyses, cleaned)


## Process main eQTL files
eQTLfiles <- c(
    'hippo_full' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_hippo_4features.rda',
    'dlpfc_full' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda',
    'interaction_full' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_hippo_4features.rda',
    'hippo_raggr' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/eqtl_tables/mergedEqtl_output_hippo_raggr_4features.rda',
    'dlpfc_raggr' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/eqtl_tables/mergedEqtl_output_dlpfc_raggr_4features.rda',
    'hippo_replication_GTEx' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/eqtl_tables/mergedEqtl_output_hippo_4features.rda',
    'dlpfc_replication_GTEx' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda',
    'interaction_replication_GTEx' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/eqtl_tables/mergedEqtl_output_interaction_4features.rda'
)
stopifnot(all(file.exists(eQTLfiles)))


## Looks like the id used is alwasy the libd one (even for the GTEx replication ones)
## which makes this easier
setkey(exon_name_map, libd_bsp2)

for(i in seq_len(length(eQTLfiles))) {
    f <- eQTLfiles[i]
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    
    ## Keep only some columns: leave the annotation to the other tables
    DT <- data.table(as.data.frame(allEqtl[, c('snps', 'gene', 'statistic', 'pvalue', 'FDR', 'beta', 'Type')]))
    
    ## So the column matches the snpAnno column name
    colnames(DT)[1] <- 'snp'
    
    ## So the column matches the feature annotation column
    colnames(DT)[2] <- 'feature_id'

    message(paste(Sys.time(), 'fixing exon ids'))
    e <- DT$feature_id[DT$Type == 'Exon']

    ## Use data.table syntax ^^
    DT$feature_id[DT$Type == 'Exon'] <- exon_name_map[.(e), gencode]
    
    ## Make Type lowercase to match file names from other tables
    DT$Type <- tolower(DT$Type)
    
    f_new <- paste0('BrainSeqPhaseII_eQTL_', names(eQTLfiles[i]), '.txt')
    message(paste(Sys.time(), 'writing', f_new))
    fwrite(DT, file = f_new, sep = '\t', row.names = FALSE)
    
    ## Clean up before the next one
    rm(DT, f, f_new, e, allEqtl)
}

message(paste(Sys.time(), 'lines for each file'))
system('wc -l BrainSeqPhaseII*')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
