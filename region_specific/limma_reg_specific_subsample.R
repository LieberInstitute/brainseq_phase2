library('SummarizedExperiment')
library('getopt')
library('limma')
library('edgeR')
library('devtools')

## Specify parameters
spec <- matrix(c(
    'type', 't', 1, 'character', 'Either gene, exon, tx or jxn',
    'age', 'a', 1, 'character', 'Either adult or fetal',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## For testing
if(FALSE){
    opt <- list('type' = 'gene', 'age' = 'fetal')
}

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

stopifnot(opt$age %in% c('adult', 'fetal'))
source('load_funs.R')
dir.create('rda', showWarnings = FALSE)
dir.create('pdf', showWarnings = FALSE)

## Load data
rse <- load_foo(opt$type, opt$age)

## To simplify later code
pd <- as.data.frame(colData(rse))
pd <- pd[, match(c('Age', 'Sex', 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5', 'Region', 'Race', 'mean_mitoRate', 'mean_totalAssignedGene'), colnames(pd))]

## Interaction model
mods <-  get_mods( colData(rse) )
sapply(mods, colnames)

## Get pieces needed for running duplication correlation
brnum <- colData(rse)$BrNum
design <- mods$mod
stopifnot(is.fullrank(design))

if(opt$type != 'tx') {
    dge <- DGEList(counts = assays(rse)$counts)
    dge <- calcNormFactors(dge)
    pdf(paste0('pdf/limma_region_specific_', opt$age, '_', opt$type, '.pdf'))
    v <- voom(dge, design, plot = TRUE)
    dev.off()

    system.time( corfit <- duplicateCorrelation(v$E, design[, c('(Intercept)',
        'RegionHIPPO')], block=brnum) )

    ## Main fit steps
    system.time( fit <- lmFit(v, design, block=brnum,
        correlation = corfit$consensus.correlation) )

    exprsNorm <- v$E
} else {
    system.time( corfit <- duplicateCorrelation(log2(assays(rse)$tpm + 0.5),
        design[, c('(Intercept)', 'RegionHIPPO')], block=brnum) )

    ## Main fit steps
    system.time( fit <- lmFit(log2(assays(rse)$tpm + 0.5), design, block=brnum,
        correlation = corfit$consensus.correlation) )
    exprsNorm <- log2(assays(rse)$tpm + 0.5)
}
system.time( fit <- eBayes(fit) )

print('Consensus correlation and summary (also after tanh transform)')
corfit$consensus.correlation
summary(corfit$atanh.correlations)
summary(tanh(corfit$atanh.correlations))


## Extract top results
colnames(design)[grep('Region', colnames(design))]
top <- topTable(fit, coef = grep('Region', colnames(design)), n = nrow(rse),
    sort.by = 'none')

save(corfit, fit, top, exprsNorm,
    file = paste0('rda/limma_region_specific_', opt$age, '_', opt$type, '.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
