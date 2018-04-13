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
        ## Drop those 4 exons not present in BrainSpan
        rse <- rse_exon[-c(175584, 175585, 175586, 175604), ]
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

load_span <- function(type, age) {
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

## Define models
fm_mod <- ~Region + Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN
fm_mod0 <- ~Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN


get_mods <- function(pd) {
    mod = model.matrix(fm_mod, data=pd)
    mod0 = model.matrix(fm_mod0, data=pd)

    return(list(mod = mod, mod0 = mod0))
}
