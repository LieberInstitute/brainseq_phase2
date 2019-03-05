library('SummarizedExperiment')
library('getopt')
library('limma')
library('edgeR')
library('devtools')
library('jaffelab')

## Specify parameters
spec <- matrix(c(
    'type', 't', 1, 'character', 'Either gene, exon, tx or jxn',
    'region', 'r', 1, 'character', 'Either DLPFC or HIPPO',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## For testing
if(FALSE){
    opt <- list('type' = 'gene', 'region' = 'DLPFC')
}

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Load data
load_foo <- function(type, region) {
    load_file <- file.path(
        '/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff',
        paste0('rse_', type, '.Rdata'))
    stopifnot(file.exists(load_file))
    load(load_file)
    
    ## Get the appropriate object
    if(type == 'gene') {
        rse <- rse_gene
    } else if (type == 'exon') {
        ## Drop those 4 exons not present in BrainSpan
        rse <- rse_exon[-c(175584, 175585, 175586, 175604), ]
    } else if (type == 'jxn') {
        rse <- rse_jxn
    } else if (type == 'tx') {
        rse <- rse_tx
    }
    ## Keep controls only
    rse <- rse[, colData(rse)$Dx == 'Control' & colData(rse)$Region == region]
    
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
    colData(rse)$Race <- relevel(factor(colData(rse)$Race), ref = 'CAUC')
    colData(rse)$Sex <- relevel(factor(colData(rse)$Sex), ref = 'F')

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
    colData(rse)$adult <- adult
    
    ## Add means
    colData(rse)$mean_mitoRate <- mean(colData(rse)$mitoRate)
    colData(rse)$mean_totalAssignedGene <- mean(colData(rse)$totalAssignedGene)
    ## Makes the design matrix not full rank in one of the models
#    colData(rse)$mean_rRNA_rate <- mean(colData(rse)$rRNA_rate)
    colData(rse)$mean_RIN <- mean(colData(rse)$RIN)
    
    return(rse)
}

rse <- load_foo(opt$type, opt$region)

## To simplify later code
pd <- as.data.frame(colData(rse))
pd <- pd[, match(c('Age', 'fetal', 'birth', 'infant', 'child', 'teen', 'adult', 'Sex', 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5', 'Region', 'Race', 'mean_mitoRate', 'mean_totalAssignedGene'), colnames(pd))]

## Define models
fm_mod <-  ~Age + fetal + birth + infant + child + teen + adult + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN
fm_mod0 <- ~ Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN

get_mods <- function(pd, int=FALSE) {    
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
mods <-  get_mods( colData(rse), int = FALSE)
sapply(mods, colnames)

## Get pieces needed for running duplication correlation
design <- mods$mod
stopifnot(is.fullrank(design))

if(opt$type != 'tx') {
    dge <- DGEList(counts = assays(rse)$counts)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot = FALSE)
    system.time( fit <- lmFit(v, design) )
    exprsNorm <- v$E
} else {
    system.time( fit <- lmFit(log2(assays(rse)$tpm + 0.5), design ) )
    exprsNorm <- log2(assays(rse)$tpm + 0.5)
}
system.time( fit <- eBayes(fit) )

## Extract top results
top <- topTable(fit, coef = 2:8, n = nrow(rse), sort.by = 'none')

## corr to age
cc = t(cor(matrix(pd$Age, nc=1), t(exprsNorm)))
top$ageCorr = cc

## fetal up/down

if(opt$type != 'tx') {
    dge <- DGEList(counts = assays(rse)$counts)
    dge <- calcNormFactors(dge)
    v <- voom(dge, model.matrix(~rse$fetal), plot = FALSE)
    system.time( fitFetal <- lmFit(v, design) )
    exprsNorm <- v$E
} else {
    system.time( fitFetal <- lmFit(log2(assays(rse)$tpm + 0.5), design ) )
    exprsNorm <- log2(assays(rse)$tpm + 0.5)
}
system.time( fitFetal <- eBayes(fitFetal) )
topFetal = topTable(fitFetal, coef = 2, n = nrow(rse), sort.by = 'none')

top$fetal_logFC = topFetal$logFC
top$fetal_t = topFetal$t
top$fetal_pval = topFetal$P.Value
top$fetal_qval = topFetal$adj.P.Val

save(fit, top, exprsNorm,
    file = paste0('rda/limma_dev_', opt$type, '_only', opt$region, '.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
