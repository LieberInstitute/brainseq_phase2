## Based on https://github.com/LieberInstitute/brainseq_phase2/blob/master/development/limma_dev.R
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

## Load DNAm pd prop
load('rda/methprop_pd.Rdata', verbose = TRUE)
meth_pd <- pd
meth_pd$Brain.Region <- toupper(meth_pd$Brain.Region)

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
        ## Drop those 4 exons not present in BrainSpan
        rse <- rse_exon[-c(175584, 175585, 175586, 175604), ]
    } else if (type == 'jxn') {
        rse <- rse_jxn
    } else if (type == 'tx') {
        rse <- rse_tx
    }
    ## Keep controls only
    rse <- rse[, colData(rse)$Dx == 'Control']
    
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
    
    colData(rse)$agegroup <- factor(with(colData(rse), dplyr::case_when(
            Age < 0 ~ 'Prenatal',
            Age >= 0 & Age < 1 ~ 'Infant',
            Age >= 1 & Age < 10 ~ 'Child',
            Age >= 10 & Age < 20 ~ 'Teen',
            Age >= 20 & Age < 50 ~ 'Adult',
            Age >= 50 ~ 'OlderAdult'
    )), levels = c('Prenatal', 'Infant', 'Child', 'Teen', 'Adult', 'OlderAdult'))
    
    ## Add cell type info
    m_meth <- match(
        paste0(colData(rse)$Region, '-', colData(rse)$BrNum),
        paste0(meth_pd$Brain.Region, '-', meth_pd$BrNum)
    )
    colData(rse)$NeuN <- meth_pd$NeuN_pos[m_meth]
    
    print(with(colData(rse), addmargins(table('Age group' = agegroup, 'Missing NeuN' = is.na(m_meth)))))
#             Missing NeuN
# Age group    FALSE TRUE Sum
#   Prenatal      23   33  56
#   Infant        19    9  28
#   Child         16    6  22
#   Teen          60   17  77
#   Adult        270  109 379
#   OlderAdult   218  120 338
#   Sum          606  294 900
    
    ## Drop those that have no cell type estimation info
    rse <- rse[, !is.na(m_meth)]
    
    return(rse)
}

rse <- load_foo(opt$type)

## Number of samples when loading genes
# rse <- load_foo('gene')
#             Missing NeuN
# Age group    FALSE TRUE Sum
#   Prenatal      23   33  56
#   Infant        19    9  28
#   Child         16    6  22
#   Teen          56   14  70
#   Adult        186   68 254
#   OlderAdult   112   72 184
#   Sum          412  202 614

with(colData(rse), addmargins(table('Age group' = agegroup, 'Region' = Region)))
#             Region
# Age group    DLPFC HIPPO Sum
#   Prenatal       8    15  23
#   Infant         9    10  19
#   Child          7     9  16
#   Teen          25    31  56
#   Adult         85   101 186
#   OlderAdult    55    57 112
#   Sum          189   223 412

tapply(colData(rse)$NeuN, colData(rse)$Region, summary)
# $DLPFC
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  0.2302  0.2918  0.2760  0.3300  0.4264
#
# $HIPPO
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.01978 0.13316 0.16724 0.17122 0.19662 0.37130

## To simplify later code
pd <- as.data.frame(colData(rse))
pd <- pd[, match(c('Age', 'fetal', 'birth', 'infant', 'child', 'teen', 'adult', 'Sex', 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5', 'Region', 'Race', 'mean_mitoRate', 'mean_totalAssignedGene', 'NeuN'), colnames(pd))]

## Define models
fm_mod <-  ~Age + fetal + birth + infant + child + teen + adult + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN + NeuN
fm_mod0 <- ~ Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN + NeuN
fm_mod_all <- ~Age *Region + fetal * Region + birth *Region + infant *Region + child * Region + teen * Region + adult * Region + Sex + Region + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN + NeuN
fm_mod0_all <- ~ Age + fetal + birth + infant + child + teen + adult + Sex + Region + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN + NeuN


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
brnum <- colData(rse)$BrNum
design <- mods$mod
stopifnot(is.fullrank(design))

if(opt$type != 'tx') {
    dge <- DGEList(counts = assays(rse)$counts)
    dge <- calcNormFactors(dge)
    pdf(paste0('pdf/limma_dev_interaction_adjNeunProp_', opt$type, '.pdf'))
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
colnames(design)[grep(':', colnames(design))]
top <- topTable(fit, coef = grep(':', colnames(design)), n = nrow(rse),
    sort.by = 'none')

save(corfit, fit, top, exprsNorm,
    file = paste0('rda/limma_dev_interaction_adjNeunProp_', opt$type, '.Rdata'))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
