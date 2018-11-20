## Load plink before starting R:
# module load plink/1.90b6.6
## Also load twas fusion code
# module load fusion_twas/github

library('SummarizedExperiment')
library('jaffelab')
library('data.table')
library('sessioninfo')
library('getopt')

## Specify parameters
spec <- matrix(c(
    'region', 'r', 1, 'character', 'Either DLPFC or HIPPO',
	'feature', 'f', 1, 'character', 'One of: gene, exon, jxn, tx',
    'cores', 'c', 1, 'integer', 'Number of cores to use. Use a small number',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## For testing
if(FALSE) {
    opt <- list(region = 'HIPPO', feature = 'gene', cores = 1)
    feat = opt$feature; reg = opt$reg
}

stopifnot(opt$region %in% c('HIPPO', 'DLPFC'))
stopifnot(opt$feature %in% c('gene', 'exon', 'jxn', 'tx'))

dir.create(opt$region, showWarnings = FALSE)
dir.create(file.path(opt$region, opt$feature), showWarnings = FALSE)

load_rse <- function(feat, reg) {
    message(paste(Sys.time(), 'loading expression data'))
    if(feat == 'gene') {
        load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
        rse <- rse_gene
        assays(rse)$raw_expr <- assays(rse_gene)$rpkm
    } else if (feat == 'exon') {
        load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata", verbose = TRUE)
        rse <- rse_exon
        assays(rse)$raw_expr <- assays(rse_exon)$rpkm
    } else if (feat == 'jxn') {
        load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata", verbose = TRUE)
        rse <- rse_jxn
        assays(rse)$raw_expr <- assays(rse_jxn)$rp10m
    } else if (feat == 'tx') {
        load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata", verbose = TRUE)
        rse <- rse_tx
        assays(rse)$raw_expr <- assays(rse_tx)$tpm
    }
    
    message(paste(Sys.time(), 'subsetting to age and region data'))
    keepInd <- which(colData(rse)$Age > 13 & colData(rse)$Region == reg)
    rse <- rse[, keepInd]
    
    message(paste(Sys.time(), 'loading model pieces'))
    if(reg == 'HIPPO') {
        load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/rdas/pcs_hippo_4features_filtered_over13.rda', verbose = TRUE)
    } else if (reg == 'DLPFC') {
        load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/rdas/pcs_dlpfc_4features_filtered_over13.rda', verbose = TRUE)
    }
    
    message(paste(Sys.time(), 'building model'))
    if(feat == 'gene') {
        pcs <- genePCs
    } else if (feat == 'exon') {
        pcs <- exonPCs
    } else if (feat == 'jxn') {
        pcs <- jxnPCs
    } else if (feat == 'tx') {
        pcs <- txPCs
    }
    colData(rse) <- cbind(colData(rse), pcs)
    mod <- model.matrix(~Dx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + pcs,
    	data = colData(rse))
        
    message(paste(Sys.time(), 'cleaning expression'))
    assays(rse) <- List(
        'raw_expr' = assays(rse)$raw_expr,
        'clean_expr' = cleaningY(log2(assays(rse)$raw_expr + 1), mod, P = 2)
    )
    
    message(paste(Sys.time(), 'switch column names to BrNum'))
    stopifnot(!any(duplicated(colnames(rse))))
    colnames(rse) <- colData(rse)$BrNum
        
    return(rse)
}

rse_file <- file.path(opt$region, opt$feature, 'working_rse.Rdata')
if(!file.exists(rse_file)) {
    rse <- load_rse(opt$feature, opt$region)
    message(paste(Sys.time(), 'saving the rse file for later at', rse_file))
    save(rse, file = rse_file)
} else {
    message(paste(Sys.time(), 'loading previous rse file', rse_file))
    load(rse_file, verbose = TRUE)
}


bim_file <- paste0(
    '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered_',
    opt$region
)
bim <- fread(
    paste0(bim_file, '.bim'),
    col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2')
)
bim_gr <- GRanges(
    paste0('chr', bim$chr),
    IRanges(bim$basepair, width = 1)
)
mcols(bim_gr) <- bim[, - c('chr', 'basepair')]

## Based on http://gusevlab.org/projects/fusion/#computing-your-own-functional-weights
## they restrict to 500kb on each side
rse_window <- resize(rowRanges(rse), width(rowRanges(rse)) + 500000 * 2, fix = 'center')


if(FALSE) {
    load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2/scz2_anneal_gr_hg19.Rdata', verbose = TRUE)
    
    ## Which of the best snps did we actually observe?
    table(countOverlaps(scz2_anneal_gr_hg19, bim_gr))
    #  0  1
    # 77 31
    scz2_anneal_gr_hg19 <- scz2_anneal_gr_hg19[countOverlaps(scz2_anneal_gr_hg19, bim_gr) > 0]
    
    ## Keep just the feature windows that include the PCG SNPs that we observed
    keep_feat <- unique(subjectHits(findOverlaps(scz2_anneal_gr_hg19, rse_window)))
    rse <- rse[keep_feat, ]
    rse_window <- rse_window[keep_feat, ]
    stopifnot(nrow(rse) == length(rse_window))
}

## Number of SNPs per feature window
table(countOverlaps(rse_window, bim_gr))

## Keep only those feature windows with some SNPs nearby
keep_feat <- which(countOverlaps(rse_window, bim_gr) > 0)
rse <- rse[keep_feat, ]
rse_window <- rse_window[keep_feat, ]
stopifnot(nrow(rse) == length(rse_window))

print('Final feature dimensions:')
print(dim(rse))

## Subset to features
setwd(file.path(opt$region, opt$feature))
dir.create('snp_files', showWarnings = FALSE)
dir.create('bim_files', showWarnings = FALSE)
if(FALSE) i <- 1
lapply(seq_len(nrow(rse)), function(i) {
    
    message(paste(Sys.time(), 'processing i =', i, 'corresponding to feature', rownames(rse)[i]))
    j <- subjectHits(findOverlaps(rse_window[i], bim_gr))
    
    filt_snp <- paste0('LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered_', opt$region, '_', opt$feature, '_', i, '.txt')
    filt_bim <- gsub('.txt', '', filt_snp)
    filt_snp <- file.path('snp_files', filt_snp)
    filt_bim <- file.path('bim_files', filt_bim)
    
    fwrite(
        bim[j, 'snp'],
        file = filt_snp,
        sep = '\t', col.names = FALSE
    )
    message(paste(Sys.time(), 'running bfile extract'))
    system(paste("plink --bfile", bim_file, '--extract', filt_snp, 
    	"--make-bed --out", filt_bim, " --memory 225000"))
    
    ## Gotta edit the last column of the fam, I think.
    ## To be continued...
    ## system('head bim_files/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered_HIPPO_gene_1.fam')
})


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

system('plink --version')