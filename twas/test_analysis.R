module load fusion_twas/github


## Using summary stats from http://walters.psycm.cf.ac.uk/

## Choose chr
chr="22"
mkdir -p test_psycm

## gotta debug chr 3

for chr in {1..22}
do
    echo "*************************"
    echo ""
    echo "processing chromosome ${chr}"
    date
    echo ""

## Create summarized analysis
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.assoc_test.R \
    --sumstats /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/psycm/clozuk_pgc2.meta.reformatted.sumstats_hg38_ourname \
    --weights /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene/HIPPO_gene.pos \
    --weights_dir /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr ${chr} \
    --out test_psycm/test_psycm.${chr}.dat
    
    echo ""
    echo "making plots for chromosome ${chr}"
    date
    echo ""

## companion plotting step
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.post_process.R \
    --sumstats /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/psycm/clozuk_pgc2.meta.reformatted.sumstats_hg38_ourname \
    --input test_psycm/test_psycm.${chr}.dat \
    --out test_psycm/test_psycm.${chr}.analysis \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr ${chr} \
    --plot --locus_win 100000 --verbose 2 --plot_individual --plot_eqtl --plot_corr

done









## Choose chr
chr="21"

for chr in {1..22}
do
    echo "*************************"
    echo ""
    echo "processing chromosome ${chr}"
    date
    echo ""

## Create summarized analysis
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.assoc_test.R \
    --sumstats /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2_sumstats/PGC2.SCZ.sumstats_hg38_ourname \
    --weights /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene/HIPPO_gene.pos \
    --weights_dir /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr ${chr} \
    --out test_PGC2.SCZ.${chr}.dat
    
    echo ""
    echo "making plots for chromosome ${chr}"
    date
    echo ""

## companion plotting step
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.post_process.R \
    --sumstats /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2_sumstats/PGC2.SCZ.sumstats_hg38_ourname \
    --input test_PGC2.SCZ.${chr}.dat \
    --out test_PGC2.SCZ.${chr}.analysis \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr ${chr} \
    --plot --locus_win 100000 --verbose 2 --plot_individual --plot_eqtl --plot_corr

done

library('readr')

## Weird issue with chr 6
info_all <- lapply(c(1:5, 7:22), function(chr) {
    message(paste(Sys.time(), 'reading chromosome', chr))
    patt <- paste0('test_PGC2.SCZ.', chr, '\\..*dat')
    files <- dir(pattern = patt)
    if(length(files) == 0) {
        warning(paste('no files for chromosome', chr))
        return(NULL)
    }
    info <- lapply(files, read_tsv)
    names(info) <- gsub(paste0('test_PGC2.SCZ.', chr, '.|analysis.'), '', dir(pattern = patt))

    # table(is.na(info[['dat']]$EQTL.ID))
    # sapply(info, dim)
    
    # chr <- 5
    # ids <- unlist(sapply(info[grepl('joint', names(info))], function(x) x$ID))
    # non_na_id <- info[['dat']]$ID[!is.na(info[['dat']]$EQTL.ID)]
    # s(subset(info[['dat']], ID %in% non_na_id[which(!non_na_id %in% ids)]))
    # non_na_id[133]

    stopifnot(sum(sapply(info[grepl('joint', names(info))], nrow)) == sum(!is.na(info[['dat']]$EQTL.ID)) - sum(info[['dat']]$TWAS.P == 1, na.rm = TRUE))
    return(info)
})
names(info_all) <- c(1:5, 7:22)

info_all <- info_all[!sapply(info_all, is.null)]


library('purrr')

info_all <- purrr::map(purrr::transpose(info_all), function(x) do.call(rbind, x))
## Missing chrs: 3, 6 and 13
which(!1:22 %in% unique(info_all[['dat']]$CHR))

library('GenomicRanges')
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_GWAS_riskSNPs/eqtl_tables/mergedEqtl_output_hippo_raggr_4features.rda', verbose = TRUE)


dim(allEqtl)
summary(allEqtl$FDR)

eGene <- subset(allEqtl, Type == 'Gene' & FDR < 0.01)
length(unique(eGene$gene))

table(info_all[['joint_dropped.dat']]$ID %in% unique(eGene$gene))
table(info_all[['joint_included.dat']]$ID %in% unique(eGene$gene))
table(info_all[['dat']]$ID %in% unique(eGene$gene))


table(info_all[['joint_dropped.dat']]$ID %in% unique(allEqtl$gene))
table(info_all[['joint_included.dat']]$ID %in% unique(allEqtl$gene))
table(info_all[['dat']]$ID %in% unique(allEqtl$gene))


summary(info_all[['joint_included.dat']])
summary(subset(info_all[['joint_included.dat']], ID %in% unique(eGene$gene)))

x <- subset(info_all[['joint_included.dat']], ID %in% unique(eGene$gene))
x
y <- subset(info_all[['dat']], ID %in% x$ID)
as.data.frame(y)

as.data.frame(y)[which.min(y$TWAS.P), ]

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
