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

## Subset to features
setwd(file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, '_pgconly', ''))))

## Start from the row data
load('rse_row.Rdata')

output_status <- bpmapply(function(i, feat_id, feature, region, basepath) {
    
    if(i == 1 || i %% 1000 == 0) {
        message('*******************************************************************************')
        message(paste(Sys.time(), 'pre-processing i =', i, 'corresponding to feature', feat_id))
    }
    
    
    filt_snp <- paste0('LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered_', region, '_', feature, '_', i, '.txt')
    filt_bim <- gsub('.txt', '', filt_snp)
    filt_snp <- file.path('snp_files', filt_snp)
    filt_bim <- file.path('bim_files', filt_bim)
    
    ## Specify input/output files
    tmp_file <- file.path(basepath, 'tmp_files', paste0(feature, '_', i))
    out_file <- file.path(basepath, 'out_files', paste0(feature, '_', i))
    
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
    
}, seq_len(length(rse_row)), names(rse_row), opt$feature, opt$region, getwd(), BPPARAM = SnowParam(workers = opt$cores), SIMPLIFY = FALSE)
output_status <- unlist(output_status) 



message('*******************************************************************************')
message(paste(Sys.time(), 'summary output status (TRUE means that there is a file)'))
table(output_status)

message(paste(Sys.time(), 'saving the output_status.Rdata file'))
names(output_status) <- paste0(opt$feature, '_', seq_len(length(rse_row)))
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
    'ID' = names(rse_row)[pos_match],
    'CHR' = gsub('chr', '', seqnames(rse_row[pos_match])),
    'P0' = start(rse_row[pos_match]), 
    'P1' = end(rse_row[pos_match]),
    'geneID' = mcols(rse_row)[, var][pos_match],
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

