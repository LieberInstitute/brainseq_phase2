## Required libraries
library('derfinder')
library('BiocParallel')
library('jaffelab')
library('getopt')

## Specify parameters
spec <- matrix(c(
	'organism', 'o', 1, 'character', 'Either rn6, mm10 or human',
	'maindir', 'm', 1, 'character', 'Main directory',
	'experiment', 'e', 1, 'character', 'Experiment',
	'prefix', 'p', 1, 'character', 'Prefix',
    'paired', 'l', 1, 'logical', 'Whether the reads are paired-end or not',
    'fullcov', 'f', 1, 'logical', 'Whether to create the full coverage object or not',
    'cores', 'c', 1, 'integer', 'Number of cores to use',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

EXPNAME <- paste0(opt$experiment, "_", opt$prefix)

if(opt$fullcov) {
    ## read in pheno	
    manifest <- read.table(file.path(opt$maindir, 'samples.manifest'),
        sep = '\t', header = FALSE, stringsAsFactors = FALSE)
    info <- data.frame('SAMPLE_ID' = manifest[, ncol(manifest)],
        stringsAsFactors = FALSE)
    N <- length(info$SAMPLE_ID)

    ### add bigwig and bam files
    info$bamFile <- file.path(opt$maindir, 'HISAT2_out', paste0(info$SAMPLE_ID,
        '_accepted_hits.sorted.bam'))
    info$bwFile <- file.path(opt$maindir, 'Coverage', paste0(info$SAMPLE_ID,
        '.bw'))

    ## Chrs to use, mitocondrial chromosome has to be the last one for the code
    ## to work later on
    if (opt$organism == "rn6") {
        CHR <- c(1:20,"X","Y","MT")
    } else if (opt$organism == "mm10") {
        CHR <- paste0("chr",c(1:19,"X","Y","M"))
    } else {
        CHR <- paste0("chr",c(1:22,"X","Y","M"))
    }
    stopifnot(grepl('M', CHR[length(CHR)]))

    ###################################################################
    
    ## Uses BAM files if the bigwigs are strand specific
    strandrule <- readLines(file.path(opt$maindir,
        'inferred_strandness_pattern.txt'))
    
    if(strandrule == 'none') {
        fullCov <- fullCoverage(files = info$bwFile, chrs = CHR,
            mc.cores = opt$cores)
    } else {
        warning('Using the BAM files instead of the strand-specific BigWigs. You might want to run fullCoverage on the strand-specific BigWigs for your analysis purposes')
        fullCov <- fullCoverage(files = info$bamFile, chrs = CHR,
            mc.cores = opt$cores)
    }
    
    save(fullCov, file = file.path(opt$maindir, paste0('fullCoverage_', EXPNAME,
        '_n', N, '.rda')))
}


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
gotDevtools <- requireNamespace('devtools', quietly = TRUE)
if(gotDevtools) {
    devtools::session_info()
} else {
    sessionInfo()
}
