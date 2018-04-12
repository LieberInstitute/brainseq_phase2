library('SummarizedExperiment')
library('getopt')
library('limma')
library('edgeR')
library('devtools')

## Specify parameters
spec <- matrix(c(
    'type', 't', 1, 'character', 'Either gene, exon, tx or jxn',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## For testing
if(FALSE){
    opt <- list('type' = 'gene')
}

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Load data
load_foo <- function(type) {
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

    print('Dimensions of the data used')
    print(dim(rse))

    ## Add age linear splines
    fetal <- ifelse(colData(rse)$Age < 0, 1,0)
    birth <- colData(rse)$Age
    birth[birth < 0] <- 0 # linear spline
    infant <- colData(rse)$Age - 1
    infant[infant < 0] <- 0 # linear spline
    child <- colData(rse)$Age - 10
    child[child < 0] <- 0 # linear spline
    teen <- colData(rse)$Age - 20
    teen[teen < 0] <- 0 # linear spline
    adult <- colData(rse)$Age - 50
    adult[adult < 0] <- 0 # linear spline

    colData(rse)$fetal <- fetal
    colData(rse)$birth <- birth
    colData(rse)$infant <- infant
    colData(rse)$child <- child
    colData(rse)$teen <- teen

    ## None are adult in the BrainSpan data set
    colData(rse)$adult <- adult

    return(rse)
}

rse <- load_foo(opt$type)

## To simplify later code
pd <- as.data.frame(colData(rse))
pd <- pd[, match(c('Age', 'fetal', 'birth', 'infant', 'child', 'teen', 'adult', 'Sex', 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5', 'Region', 'Race', 'mean_mitoRate', 'mean_totalAssignedGene'), colnames(pd))]

## Define models
fm_mod <-  ~Age + fetal + birth + infant + child + teen + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN #+ adult
fm_mod0 <- ~ Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN
fm_mod_all <- ~Age *Region + fetal * Region + birth *Region + infant *Region + child * Region + teen * Region + Sex + Region + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN # + adult * Region
fm_mod0_all <- ~ Age + fetal + birth + infant + child + teen + Sex + Region + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN # + adult


get_mods <- function(pd, int = FALSE) {
    ### adjust for race, sex
    if(int) {
        mod = model.matrix(fm_mod_all, data=pd)
        mod0 = model.matrix(fm_mod0_all, data=pd)
    } else {
        mod = model.matrix(fm_mod, data=pd)
        mod0 = model.matrix(fm_mod0, data=pd)
    }
    return(list(mod = mod, mod0 = mod0))
}

## Interaction model
mods <-  get_mods( colData(rse), int = TRUE)
sapply(mods, colnames)

## Get pieces needed for running duplication correlation
brnum <- colData(rse)$Braincode
design <- mods$mod
stopifnot(is.fullrank(design))

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
    pdf(paste0('pdf/span_limma_dev_interaction_', opt$type, '.pdf'))
    v <- voom(dge, design, plot = TRUE)
    dev.off()

    system.time( corfit <- duplicateCorrelation(v$E, design[, c('(Intercept)',
        'RegionHIPPO')], block=brnum) )

    ## Main fit steps
    system.time( fit <- lmFit(v, design, block=brnum,
        correlation = corfit$consensus.correlation) )

    exprsNorm <- v$E
} else {
    system.time( corfit <- duplicateCorrelation(assays(rse)$tpm, design[,
        c('(Intercept)', 'RegionHIPPO')], block=brnum) )

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
colnames(design)[grep(':', colnames(design))]
top <- topTable(fit, coef = grep(':', colnames(design)), n = nrow(rse),
    sort.by = 'none')

save(corfit, fit, top, exprsNorm,
    file = paste0('rda/span_limma_dev_interaction_', opt$type, '.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
