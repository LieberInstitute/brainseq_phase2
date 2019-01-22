module load fusion_twas/github

##
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.assoc_test.R \
    --sumstats /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2_sumstats/PGC2.SCZ.sumstats_hg38_ourname \
    --weights /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene/HIPPO_gene.pos \
    --weights_dir /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr 22 \
    --out test_PGC2.SCZ.22.dat
    
## companion plotting step
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.post_process.R \
    --sumstats /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2_sumstats/PGC2.SCZ.sumstats_hg38_ourname \
    --input test_PGC2.SCZ.22.dat \
    --out test_PGC2.SCZ.22.analysis \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr 22 \
    --plot --locus_win 100000 --verbose 2
    

## Original goal
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.assoc_test.R \
    --sumstats /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2_sumstats/PGC2.SCZ.sumstats \
    --weights /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene/HIPPO_gene.pos \
    --weights_dir /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr 22 \
    --out test_PGC2.SCZ.22.dat
    
## Removed similar warnings for all the other chr22 genes...
# WARNING :  out_files/gene_9997.wgt.RDat ENSG00000272940.1 22 50597152 50597599 had 0 overlapping SNPs, but none with non-zero expression weights, try more SNPS or a different model
# WARNING :  out_files/gene_9997.wgt.RDat ENSG00000272940.1 22 50597152 50597599  had no overlapping SNPs
# Analysis completed.
# NOTE: 146 / 146 genes were skipped
# If a large number of genes were skipped, verify that your GWAS Z-scores, expression weights, and LDREF data use the same SNPs (or nearly)
# Or consider pre-imputing your summary statistics to the LDREF markers using summary-imputation software such as [http://bogdan.bioinformatics.ucla.edu/software/impg/]
    
## companion plotting step (doesn't work for another reason)
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.post_process.R \
    --sumstats /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2_sumstats/PGC2.SCZ.sumstats \
    --input test_PGC2.SCZ.22.dat \
    --out test_PGC2.SCZ.22.analysis \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr 22 \
    --plot --locus_win 100000

## Use filtered PGC2 sumstats to those snps present in the ported hg38 reference
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.assoc_test.R \
    --sumstats /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2_sumstats/PGC2.SCZ.sumstats_filt \
    --weights /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene/HIPPO_gene.pos \
    --weights_dir /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr 22 \
    --out test_PGC2.SCZ.22.dat

## Try with another chromosome (also doesn't work)
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.assoc_test.R \
    --sumstats /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2_sumstats/PGC2.SCZ.sumstats \
    --weights /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene/HIPPO_gene.pos \
    --weights_dir /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr 10 \
    --out test_PGC2.SCZ.10.dat

## Try with the original hg19 reference SNPs (although gene coords are in hg38)
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.assoc_test.R \
    --sumstats /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2_sumstats/PGC2.SCZ.sumstats \
    --weights /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene/HIPPO_gene.pos \
    --weights_dir /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LDREF/1000G.EUR. \
    --chr 22 \
    --out test_PGC2.SCZ.22.dat #--force_model "lasso"


        
        
load('out_files/gene_12892.wgt.RDat', verbose = TRUE)
library(data.table)
pgc2 <- fread('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2_sumstats/PGC2.SCZ.sumstats')

present <- snps$V2 %in% pgc2$SNP
table(present)
# present
# FALSE  TRUE
#   348   601
head(snps$V2[!present])

their_bims <- dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LDREF', '.*bim$', full.names = TRUE)
names(their_bims) <- dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LDREF', '.*bim$')

ldref_bim <- do.call(rbind, lapply(their_bims, function(input_bim) {
    message(paste(Sys.time(), 'reading file', input_bim))
    res <- fread(input_bim,
        col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2'),
        colClasses = c('character', 'character', 'numeric', 'integer', 'character', 'character')
    )
    #setkey(res, 'chr', 'basepair')
    return(res)
}))

present_ref <- snps$V2 %in% ldref_bim$snp
table(present_ref)
# present_ref
# FALSE  TRUE
#   289   660

their_bims_hg38 <- dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38', '.*bim$', full.names = TRUE)
names(their_bims_hg38) <- dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38', '.*bim$')

ldref_bim_hg38 <- do.call(rbind, lapply(their_bims_hg38, function(input_bim) {
    message(paste(Sys.time(), 'reading file', input_bim))
    res <- fread(input_bim,
        col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2'),
        colClasses = c('character', 'character', 'numeric', 'integer', 'character', 'character')
    )
    #setkey(res, 'chr', 'basepair')
    return(res)
}))

present_ref_hg38 <- snps$V2 %in% ldref_bim_hg38$snp
table(present_ref_hg38)
# FALSE  TRUE
#   338   611

present_z_ref <- pgc2$SNP %in% ldref_bim$snp
table(present_z_ref)
#    TRUE
# 1083014

present_z_ref_hg38 <- pgc2$SNP %in% ldref_bim_hg38$snp
table(present_z_ref_hg38)
#  FALSE   TRUE
# 114105 968909

## Ok, filter the PGC2 sumstats to the SNPs we have in the ported hg38 reference
pgc2_filt <- pgc2[present_z_ref_hg38, ]
fwrite(pgc2_filt, file = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2_sumstats/PGC2.SCZ.sumstats_filt', sep = '\t')
