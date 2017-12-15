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
    opt <- list('type' = 'gene', 'age' = 'adult')
    opt <- list('type' = 'jxn', 'age' = 'adult')
}

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

stopifnot(opt$age %in% c('adult', 'fetal'))

## Load data
load_foo <- function(type, age) {
    load_file <- file.path(
        '/dcl01/lieber/ajaffe/lab/brainseq_phase2/brainspan',
        paste0('rse_span_', type, '.Rdata'))
    stopifnot(file.exists(load_file))
    load(load_file)
    
    ## Get the appropriate object
    if(type == 'gene') {
        rse <- rse_span_gene
    } else if (type == 'exon') {
        rse <- rse_span_exon
    } else if (type == 'jxn') {
        rse <- rse_span_jxn
    } else if (type == 'tx') {
        rse <- rse_span_tx
    }
    
    ## Keep the corresponding age group
    if(age == 'adult') {
        rse <- rse[, colData(rse)$Age >= 18]
    } else if (age == 'fetal') {
        rse <- rse[, colData(rse)$Age <= 0]
    }
    
    print('Dimensions of the data used')
    print(dim(rse))
    
    return(rse)
}

rse <- load_foo(opt$type, opt$age)

## To simplify later code
pd <- as.data.frame(colData(rse))
pd <- pd[, match(c('Age', 'Sex', 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5', 'Region', 'Race', 'mean_mitoRate', 'mean_totalAssignedGene', 'mean_rRNA_rate'), colnames(pd))]

## Define models
fm_mod <- ~Region + Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_rRNA_rate + mean_RIN
fm_mod0 <- ~Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_rRNA_rate + mean_RIN


get_mods <- function(pd) {    
    mod = model.matrix(fm_mod, data=pd)
    mod0 = model.matrix(fm_mod0, data=pd)
    
    return(list(mod = mod, mod0 = mod0))
}

## Interaction model
mods <-  get_mods( colData(rse) )
sapply(mods, colnames)

## Get pieces needed for running duplication correlation
brnum <- colData(rse)$Braincode
design <- mods$mod



if(opt$type != 'tx') {
    if(opt$type == 'jxn') {
        ## as.matrix(assays(rse)$counts) fails due to some rownames length issue
        ## even after dropping the rownames, so here I build the matrix manually
        cts <- matrix(0, nrow = nrow(rse), ncol = ncol(rse))
        rownames(cts) <- rownames(rse)
        colnames(cts) <- colnames(rse)
        for(i in seq_len(ncol(rse))) cts[, i] <- as.vector(assays(rse)$counts[, i])
        rm(i)
        dge <- DGEList(counts = cts)
    } else {
        dge <- DGEList(counts = assays(rse)$counts)
    }
    
    dge <- calcNormFactors(dge)
    pdf(paste0('pdf/span_limma_region_specific_', opt$age, '_', opt$type, '.pdf'))
    v <- voom(dge, design, plot = TRUE)
    dev.off()
        
    system.time( corfit <- duplicateCorrelation(v$E, design[, c('(Intercept)',
        'RegionHIPPO')], block=brnum) )
    
    ## Main fit steps
    system.time( fit <- lmFit(v, design, block=brnum,
        correlation = corfit$consensus.correlation) )
        
    exprsNorm <- v$E
} else {
    system.time( corfit <- duplicateCorrelation(assays(rse)$tpm,
        design[, c('(Intercept)', 'RegionHIPPO')], block=brnum) )

    ## Main fit steps
    system.time( fit <- lmFit(assays(rse)$tpm, design, block=brnum,
        correlation = corfit$consensus.correlation) )
    exprsNorm <- assays(rse)$tpm
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
    file = paste0('rda/span_limma_region_specific_', opt$age, '_', opt$type, '.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
