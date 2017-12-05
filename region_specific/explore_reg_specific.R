library('limma')
library('SummarizedExperiment')
library('jaffelab')
library('devtools')

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
    
    ## Add means
    colData(rse)$mean_mitoRate <- mean(colData(rse)$mitoRate)
    colData(rse)$mean_totalAssignedGene <- mean(colData(rse)$totalAssignedGene)
    colData(rse)$mean_rRNA_rate <- mean(colData(rse)$rRNA_rate)
    colData(rse)$mean_RIN <- mean(colData(rse)$RIN)
    
    print('Dimensions of the data used')
    print(dim(rse))
    
    return(rse)
}

## Define models
fm_mod <- ~Region + Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_rRNA_rate + mean_RIN
fm_mod0 <- ~Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_rRNA_rate + mean_RIN


get_mods <- function(pd) {    
    mod = model.matrix(fm_mod, data=pd)
    mod0 = model.matrix(fm_mod0, data=pd)
    
    return(list(mod = mod, mod0 = mod0))
}





top <- mapply(function(age, type) {
    load(paste0('limma_region_specific_', age, '_', type, '.Rdata'))
    top$age <- age
    top$type <- type
    return(top)
}, rep(c('adult', 'fetal'), each = 4), rep(c('gene', 'exon', 'jxn', 'tx'), 2), SIMPLIFY = FALSE)
names(top) <- paste0(rep(c('adult', 'fetal'), each = 4), "_", rep(c('gene', 'exon', 'jxn', 'tx'), 2))

fit <- mapply(function(age, type) {
    load(paste0('limma_region_specific_', age, '_', type, '.Rdata'))
    return(fit)
}, rep(c('adult', 'fetal'), each = 4), rep(c('gene', 'exon', 'jxn', 'tx'), 2), SIMPLIFY = FALSE)

exprsNorm <- mapply(function(age, type) {
    load(paste0('limma_region_specific_', age, '_', type, '.Rdata'))
    return(exprsNorm)
}, rep(c('adult', 'fetal'), each = 4), rep(c('gene', 'exon', 'jxn', 'tx'), 2), SIMPLIFY = FALSE)

names(exprsNorm) <- names(fit) <- names(top)

sapply(top, function(x) {
    table(x$adj.P.Val < 0.05)
})

sapply(top, function(x) {
    round(table(x$adj.P.Val < 0.05) / nrow(x) * 100, 2)
})

sapply(top, function(x) {
    table(x$adj.P.Val < 0.01)
})

sapply(top, function(x) {
    round(table(x$adj.P.Val < 0.01) / nrow(x) * 100, 2)
})



sapply(top, function(x) {
    table(p.adjust(x$P.Value, 'bonferroni') < 0.05)
})

sapply(top, function(x) {
    round(table(p.adjust(x$P.Value, 'bonferroni') < 0.05) / nrow(x) * 100, 2)
})



rse <- load_foo('gene', 'adult')
rse_tx <- load_foo('tx', 'adult')
rse_exon <- load_foo('exon', 'adult')
rse_jxn <- load_foo('jxn', 'adult')


## Reg specific model
design <- get_mods( colData(rse) )$mod

min(top$adult_gene$adj.P.Val)
which.min(top$adult_gene$adj.P.Val)
which(rank(top$adult_gene$adj.P.Val) == 1)



cleanedVoom <- cleaningY(exprsNorm$adult_gene, design, 2)

boxplot(v$E[which(rank(top$adult_gene$adj.P.Val) == 1), ] ~ colData(rse)$Region)
t.test(v$E[which(rank(top$adult_gene$adj.P.Val) == 1), ] ~ colData(rse)$Region)

boxplot(cleanedVoom[which(rank(top$adult_gene$adj.P.Val) == 1), ] ~ colData(rse)$Region)
t.test(cleanedVoom[which(rank(top$adult_gene$adj.P.Val) == 1), ] ~ colData(rse)$Region)

top$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 2), ]
boxplot(cleaned[which(rank(top$adult_gene$adj.P.Val) == 2), ] ~ colData(rse)$Region)




## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()