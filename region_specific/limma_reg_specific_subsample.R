library('SummarizedExperiment')
library('getopt')
library('limma')
library('edgeR')
library('sessioninfo')

## Specify parameters
spec <- matrix(c(
    'type', 't', 1, 'character', 'Either gene, exon, tx or jxn',
    'age', 'a', 1, 'character', 'Either adult or fetal',
    'iteration', 'i', 1, 'numeric', 'An interation number',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## For testing
if(FALSE){
    opt <- list('type' = 'gene', 'age' = 'adult', 'iteration' = 1)
}

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

stopifnot(opt$age %in% c('adult', 'fetal'))
source('load_funs.R')
dir.create('subsample', showWarnings = FALSE)
dir.create('subsample/rda', showWarnings = FALSE)
dir.create('subsample/pdf', showWarnings = FALSE)

## Load data
rse <- load_foo(opt$type, opt$age)

## With:
# opt <- list('type' = 'gene', 'age' = 'fetal', 'iteration' = 1)
##
# > table(rse$Region)
#
# DLPFC HIPPO
#    28    28

## Subsample
index <- split(seq_len(ncol(rse)), rse$Region)
message(paste(Sys.time(), 'selected samples'))
set.seed(opt$iteration + 20190312)
new_index <- unlist(purrr::map(index, sample, size = 28))
new_index

message(paste(Sys.time(), 'final subsetted dimensions'))
rse <- rse[, new_index]
stopifnot(all(table(rse$Region) == 28))
dim(rse)


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
    pdf(paste0('subsample/pdf/limma_region_specific_', opt$age, '_', opt$type, '_', opt$iteration, '.pdf'))
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
    file = paste0('subsample/rda/limma_region_specific_', opt$age, '_', opt$type, '_', opt$iteration, '.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
