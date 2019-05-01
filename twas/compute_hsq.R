## Adapted from compute_weights.R

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
dir.create(file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, '_pgconly', '')), 'heritability'), showWarnings = FALSE)

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


rse_file <- file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, '_pgconly', '')), 'subsetted_rse.Rdata')
message(paste(Sys.time(), 'loading', rse_file))
load(rse_file, verbose = TRUE)
print('RSE dimensions:')
dim(rse)

## Based on http://gusevlab.org/projects/fusion/#computing-your-own-functional-weights
## they restrict to 500kb on each side
rse_window <- resize(rowRanges(rse), width(rowRanges(rse)) + 500000 * 2, fix = 'center')

## Subset to features
setwd(file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, '_pgconly', '')), 'heritability'))
dir.create('snp_files', showWarnings = FALSE)
dir.create('bim_files', showWarnings = FALSE)
dir.create('tmp_files', showWarnings = FALSE)
dir.create('out_files', showWarnings = FALSE)

## For testing
if(FALSE) rse <- rse[1:20, ]

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
        'Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/FUSION.compute_hsq.R --bfile',
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
    
    return(file.exists(paste0(out_file, '.hsq.Rdata')))
    
}, seq_len(nrow(rse)), rownames(rse), BPPARAM = MulticoreParam(workers = opt$cores), SIMPLIFY = FALSE)
output_status <- unlist(output_status)

message('*******************************************************************************')
message(paste(Sys.time(), 'summary output status (TRUE means that there is a file)'))
table(output_status)

message(paste(Sys.time(), 'saving the output_status.Rdata file'))
names(output_status) <- paste0(opt$feature, '_', seq_len(length(rowRanges(rse))))
save(output_status, file = 'output_status.Rdata')

message(paste(Sys.time(), 'creating the heritability table file'))
rdat_files <- dir('out_files', '.hsq.Rdata', full.names = TRUE)
stopifnot(length(rdat_files) == sum(output_status))

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

pos_match <- match(gsub('out_files/|\\.hsq\\.Rdata', '', rdat_files), names(output_status))
stopifnot(all(!is.na(pos_match)))

stopifnot(identical(gsub('out_files/|\\.hsq\\.Rdata', '', rdat_files), names(output_status)[pos_match]))


hsq_info <- data.frame(
    'HSQ' = rdat_files,
    ## use the feature id
    'ID' = names(rowRanges(rse))[pos_match],
    'CHR' = gsub('chr', '', seqnames(rowRanges(rse)[pos_match])),
    'P0' = start(rowRanges(rse)[pos_match]), 
    'P1' = end(rowRanges(rse)[pos_match]),
    'geneID' = mcols(rowRanges(rse))[, var][pos_match],
    'hsq' = NA,
    'hsq.se' = NA,
    'hsq.pv' = NA,
    'region' = opt$region,
    'feature' = opt$feature,
    stringsAsFactors = FALSE
)

hsq_info$ID[hsq_info$ID == ''] <- NA
hsq_info$geneID[hsq_info$geneID == ''] <- NA

## Adapted from https://github.com/LieberInstitute/fusion_twas/blob/jhpce/utils/FUSION.profile_wgt.R
for ( i in seq_len(nrow(hsq_info)) ) {
    load(hsq_info$HSQ[i])
    hsq_info$hsq[i] <- hsq[1]
    hsq_info$hsq.se[i] <- hsq[2]
    hsq_info$hsq.pv[i] <- hsq.pv
}

## Does it fail the filters at
## https://github.com/LieberInstitute/fusion_twas/blob/jhpce/FUSION.compute_weights.R#L238 ?

hsq_info$failFilter <- hsq_info$hsq < 0 | hsq_info$hsq.pv > 0.01

print('hsq_info summary information')
sapply(hsq_info, function(x) sum(is.na(x)))
summary(hsq_info)

save(hsq_info, file = 'hsq_info.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

message('-- plink version information --')
system('plink --version')

