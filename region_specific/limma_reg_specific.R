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

## Load data
load_foo <- function(type, age) {
    load_file <- file.path(
        '/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff',
        paste0('rse_', type, '.Rdata'))
    stopifnot(file.exists(load_file))
    load(load_file)
    
    ## Get the appropriate object
    if(type == 'gene') {
        rse <- rse_gene
    } else if (type == 'exon') {
        rse <- rse_exon
    } else if (type == 'jxn') {
        rse <- rse_jxn
    } else if (type == 'tx') {
        rse <- rse_tx
    }
    ## Keep controls only
    rse <- rse[, colData(rse)$Dx == 'Control']
    
    ## Keep the corresponding age group
    if(age == 'adult') {
        rse <- rse[, colData(rse)$Age >= 18]
    } else if (age == 'fetal') {
        rse <- rse[, colData(rse)$Age <= 0]
    }
    
    ## Add mds info
     load(file.path('/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data', 
         'mds_extracted_from_BrainSeq_Phase2_RiboZero_Genotypes_n551.Rdata'))
    m <- match(colData(rse)$BrNum, rownames(mds))
    print('Number of missing brains in the MDS data')
    print(table(is.na(m)))

    ## Drop those that don't match
    if(any(is.na(m))) {
        print(colData(rse)$BrNum[which(is.na(m))])
    }
    rse <- rse[, !is.na(m)]
    colData(rse) <- cbind(colData(rse), mds[m[!is.na(m)], ])
    
    ## Set as factor
    colData(rse)$Region <- relevel(factor(colData(rse)$Region), 'DLPFC')
    colData(rse)$Race <- relevel(factor(colData(rse)$Race), ref = 'CAUC')
    colData(rse)$Sex <- relevel(factor(colData(rse)$Sex), ref = 'F')
    
    print('Dimensions of the data used')
    print(dim(rse))
    
    return(rse)
}

rse <- load_foo(opt$type, opt$age)

## To simplify later code
pd <- as.data.frame(colData(rse))
pd <- pd[, match(c('Age', 'Sex', 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5', 'Region', 'Race'), colnames(pd))]

## Define models
fm_mod <- ~Region + Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5
fm_mod0 <- ~Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5


get_mods <- function(pd) {    
    mod = model.matrix(fm_mod, data=pd)
    mod0 = model.matrix(fm_mod0, data=pd)
    
    return(list(mod = mod, mod0 = mod0))
}

## Interaction model
mods <-  get_mods( colData(rse), int = TRUE)
sapply(mods, colnames)


design <- mods$mod

if(opt$type != 'tx') {
    dge <- DGEList(counts = assays(rse)$counts)
    dge <- calcNormFactors(dge)
    pdf(paste0('limma_region_specific_', opt$age, '_', opt$type, '.pdf'))
    v <- voom(dge, design, plot = TRUE)
    dev.off()
        
    ## Main fit steps
    system.time( fit <- lmFit(v, design) )
} else {
    ## Main fit steps
    system.time( fit <- lmFit(assays(rse)$tpm, design) )
}
system.time( fit <- eBayes(fit) )


## Extract top results
colnames(design)[grep('Region', colnames(design))]
top <- topTable(fit, coef = grep('Region', colnames(design)), n = nrow(rse),
    sort.by = 'none')

save(fit, top,
    file = paste0('limma_region_specific_', opt$age, '_', opt$type, '.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
