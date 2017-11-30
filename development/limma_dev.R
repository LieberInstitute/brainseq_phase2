library('SummarizedExperiment')
library('recount')
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

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Load data
load_foo <- function(type) {
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
        stop('Type tx is not supported yet')
        rse <- rse_tx
    }
    ## Keep controls only
    rse <- rse[, colData(rse)$Dx == 'Control']
    
    ## Add mds info
     load(file.path('/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data', 
         'mds_extracted_from_BrainSeq_Phase2_RiboZero_Genotypes_n551.Rdata'))
    m <- match(colData(rse)$BrNum, rownames(mds))

    ## Drop those that don't match
    colData(rse)$BrNum[which(is.na(m))]

    rse <- rse[, !is.na(m)]
    colData(rse) <- cbind(colData(rse), mds[m[!is.na(m)], ])
    
    ## Set as factor
    colData(rse)$Region <- relevel(factor(colData(rse)$Region), 'DLPFC')
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
}

rse <- load_foo(opt$type)

## To simplify later code
pd <- as.data.frame(colData(rse))
pd <- pd[, match(c('Age', 'fetal', 'birth', 'infant', 'child', 'teen', 'adult', 'Sex', 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5', 'Region', 'Race'), colnames(pd))]

## Define models
fm_mod <-  ~Age + fetal + birth + infant + child + teen + adult + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5
fm_mod0 <- ~ Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5
fm_mod_all <- ~Age *Region + fetal * Region + birth *Region + infant *Region + child * Region + teen * Region + adult * Region + Sex + Region + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5
fm_mod0_all <- ~ Age + fetal + birth + infant + child + teen + adult + Sex + Region + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5


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
ind <- data.frame(
    brnum = colData(rse)$BrNum,
    rnum = colData(rse)$RNum,
    stringsAsFactors = FALSE
)

brnum <- ind$brnum[match(rownames(pd), ind$rnum)]
design <- mods$mod

dge <- DGEList(counts = assays(rse)$counts)
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot = TRUE)
system.time( corfit <- duplicateCorrelation(v$E, design, block=brnum) )

print('Consensus correlation and summary (also after tanh transform)')
corfit$consensus.correlation
summary(corfit$atanh.correlations)
summary(tanh(corfit$atanh.correlations))

## Main fit steps
system.time( fit <- lmFit(v, design, block=brnum,
    correlation = corfit$consensus.correlation) )
system.time( fit <- eBayes(fit) )

## Extract top results
colnames(design)[grep(':', colnames(design))]
top <- topTable(fit, coef = grep(':', colnames(design)), n = nrow(rse),
    sort.by = 'none')

save(corfit, fit, top,
    file = paste0('limma_dev_interaction_', opt$type, '.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
