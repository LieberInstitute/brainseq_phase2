## Load plink before starting R:
# module load plink/1.90b6.6
## Also load twas fusion code
# module load fusion_twas/github

library('SummarizedExperiment')
library('jaffelab')
library('data.table')
library('sessioninfo')
library('getopt')
library('BiocParallel')

## Without this, the memory use blows up
## getDTthreads() will detect 64 threads ins ome cases here
getDTthreads()
setDTthreads(threads = 1)
getDTthreads()

## Specify parameters
spec <- matrix(c(
    'region', 'r', 1, 'character', 'Either DLPFC or HIPPO',
	'feature', 'f', 1, 'character', 'One of: gene, exon, jxn, tx',
    'cores', 'c', 1, 'integer', 'Number of cores to use. Use a small number',
    'pgconly', 'p', 1, 'logical', 'Subset to only PGC loci?',
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
    opt <- list(region = 'HIPPO', feature = 'gene', cores = 1, 'pgconly' = TRUE)
    opt <- list(region = 'HIPPO', feature = 'gene', cores = 1, 'pgconly' = FALSE)
    opt <- list(region = 'DLPFC', feature = 'gene', cores = 1, 'pgconly' = FALSE)
    opt <- list(region = 'HIPPO', feature = 'gene', cores = 3, 'pgconly' = FALSE)
    # feat = opt$feature; reg = opt$reg
    
    opt <- list(region = 'HIPPO', feature = 'exon', cores = 3, 'pgconly' = FALSE)
}

stopifnot(opt$region %in% c('HIPPO', 'DLPFC'))
stopifnot(opt$feature %in% c('gene', 'exon', 'jxn', 'tx'))

dir.create(opt$region, showWarnings = FALSE)
dir.create(file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, '_pgconly', ''))), showWarnings = FALSE)

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

rse_file <- file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, '_pgconly', '')), 'working_rse.Rdata')
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


if(opt$pgconly) {
    load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2/scz2_anneal_gr_hg19.Rdata', verbose = TRUE)

    bim_hg19 <- fread(
        paste0(bim_file, '.bim.original'),
        col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2')
    )
    bim_gr_hg19 <- GRanges(
        paste0('chr', bim_hg19$chr),
        IRanges(bim_hg19$basepair, width = 1)
    )
    mcols(bim_gr_hg19) <- bim_hg19[, - c('chr', 'basepair')]
    
    ## Which of the best snps did we actually observe?
    table(countOverlaps(scz2_anneal_gr_hg19, bim_gr_hg19))
    #  0  1
    # 77 31
    scz2_anneal_gr_hg19 <- scz2_anneal_gr_hg19[countOverlaps(scz2_anneal_gr_hg19, bim_gr_hg19) > 0]
    
    ## Would need to fix this part for this to work (since the rse coordinates are in hg38...)
    
    ## Keep just the feature windows that include the PCG SNPs that we observed
    keep_feat <- unique(subjectHits(findOverlaps(scz2_anneal_gr_hg19, rse_window)))
    rse <- rse[keep_feat, ]
    rse_window <- rse_window[keep_feat, ]
    stopifnot(nrow(rse) == length(rse_window))
}

## Number of SNPs per feature window
# table(countOverlaps(rse_window, bim_gr))

## Keep only those feature windows with some SNPs nearby
keep_feat <- which(countOverlaps(rse_window, bim_gr) > 0)
message(paste(Sys.time(), 'number of features kept:', length(keep_feat)))
rse <- rse[keep_feat, ]
rse_window <- rse_window[keep_feat, ]
stopifnot(nrow(rse) == length(rse_window))

print('Final RSE feature dimensions:')
print(dim(rse))

rse_file <- file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, '_pgconly', '')), 'subsetted_rse.Rdata')
message(paste(Sys.time(), 'saving the subsetted rse file for later at', rse_file))
save(rse, file = rse_file)

## Subset to features
setwd(file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, '_pgconly', ''))))
dir.create('snp_files', showWarnings = FALSE)
dir.create('bim_files', showWarnings = FALSE)
dir.create('tmp_files', showWarnings = FALSE)
dir.create('out_files', showWarnings = FALSE)

## For testing
if(FALSE) rse <- rse[1:20, ]
    rse <- rse[1:1000, ]

## Work in parallel
output_status <- bpmapply(function(i, feat_id) {
    
    if(i == 1 || i %% 1000 == 0) {
        message('*******************************************************************************')
        message(paste(Sys.time(), 'pre-processing i =', i, 'corresponding to feature', feat_id))
    }
    
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
    
    system(paste("plink --bfile", bim_file, '--extract', filt_snp, 
    	"--make-bed --out", filt_bim, '--memory 2000 --threads 1 --silent'))
    
    ## Edit the "phenotype" column of the fam file
    filt_fam <- fread(paste0(filt_bim, '.fam'),
        col.names = c('famid', 'w_famid', 'w_famid_fa', 'w_famid_mo', 'sex_code', 'phenotype')
    )
    
    ## Note BrNums might be duplicated, hence the use of match()
    m <- match(filt_fam$famid, colData(rse)$BrNum)
    
    ## Use cleaned expression for now. Could be an argument for the code.
    filt_fam$phenotype <- assays(rse)$clean_expr[i, m]
    
    ## Ovewrite fam file (for the phenotype info)
    fwrite(filt_fam, file = paste0(filt_bim, '.fam'), sep = ' ', col.names = FALSE)    
    
    
    ## Specify input/output files
    tmp_file <- file.path(getwd(), 'tmp_files', paste0(opt$feature, '_', i))
    out_file <- file.path(getwd(), 'out_files', paste0(opt$feature, '_', i))
    
    system(paste(
        'Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.compute_weights.R --bfile',
        filt_bim,
        '--tmp', tmp_file,
        '--out', out_file,
        '--models top1,blump,lasso,enet'
        ## Using --noclean TRUE for testing right now
        # '--models top1,blump,lasso,enet --noclean TRUE --hsq_p 0.9'
    ))   
    
    ## Clean up
    unlink(filt_snp)
    system(paste0('rm ', filt_bim, '.*'))
    
    ## The first five genes (pgconly = TRUE) failed the hsq_p default of 0.01 https://github.com/gusevlab/fusion_twas/blob/master/FUSION.compute_weights.R#L26
    ## Just tested the 5th one with hsq_p 0.9 and it works
    # '--models top1,blump,lasso,enet --noclean TRUE --hsq_p 0.9'
    
    return(file.exists(paste0(out_file, '.wgt.RDat')))
    
}, seq_len(nrow(rse)), rownames(rse), BPPARAM = MulticoreParam(workers = opt$cores), SIMPLIFY = FALSE)
output_status <- unlist(output_status)

message('*******************************************************************************')
message(paste(Sys.time(), 'summary output status (TRUE means that there is a file)'))
table(output_status)

message(paste(Sys.time(), 'saving the output_status.Rdata file'))
names(output_status) <- paste0(opt$feature, '_', seq_len(length(rowRanges(rse))))
save(output_status, file = 'output_status.Rdata')

message(paste(Sys.time(), 'creating the wgt profile files'))
rdat_files <- dir('out_files', '.wgt.RDat', full.names = TRUE)
stopifnot(length(rdat_files) == sum(output_status))

wglist <- paste0(opt$region, '_', opt$feature, '.list')
write.table(rdat_files, file = wglist, row.names = FALSE, col.names = FALSE, quote = FALSE)
system(paste0('Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/utils/FUSION.profile_wgt.R ', wglist, ' > ', opt$region, '_', opt$feature, '.profile.err 2>&1'))
## Rename the wglist_summary.txt file to keep the naming convention consistent
system(paste0('mv wglist_summary.txt ', opt$region, '_', opt$feature, '.profile'))

message(paste(Sys.time(), 'creating the .pos file'))
if(opt$feature %in% c('gene', 'exon')) {
    var <- 'gencodeID'
    # vars <- 'Symbol'
} else if (opt$feature == 'jxn') {
    var <- 'gencodeGeneID'
    # vars <- 'Symbol'
} else {
    var <- 'gene_id'
    # vars <- 'gene_name'
}

pos_match <- match(gsub('out_files/|\\.wgt\\.RDat', '', rdat_files), names(output_status))
stopifnot(all(!is.na(pos_match)))

stopifnot(identical(gsub('out_files/|\\.wgt\\.RDat', '', rdat_files), names(output_status)[pos_match]))


pos_info <- data.frame(
    'WGT' = rdat_files,
    ## use the feature id
    'ID' = names(rowRanges(rse))[pos_match],
    'CHR' = gsub('chr', '', seqnames(rowRanges(rse)[pos_match])),
    'P0' = start(rowRanges(rse)[pos_match]), 
    'P1' = end(rowRanges(rse)[pos_match]),
    'geneID' = mcols(rowRanges(rse))[, var][pos_match],
    stringsAsFactors = FALSE
)
sapply(pos_info, function(x) sum(is.na(x)))
pos_info$ID[pos_info$ID == ''] <- NA
pos_info$geneID[pos_info$geneID == ''] <- NA
write.table(pos_info[, -which(colnames(pos_info) == 'geneID')],
    file = paste0(opt$region, '_', opt$feature, '.pos'),
    row.names = FALSE, col.names = TRUE, quote = FALSE
)
save(pos_info, file = 'pos_info.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

message('-- plink version information --')
system('plink --version')

