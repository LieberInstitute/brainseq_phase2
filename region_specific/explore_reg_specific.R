library('limma')
library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('ggplot2')

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

## Load BrainSeq model results
raw <- mapply(function(age, type) {
    load(paste0('rda/limma_region_specific_', age, '_', type, '.Rdata'))
    top$age <- age
    top$type <- type
    return(list(top = top, fit = fit, exprsNorm = exprsNorm))
}, rep(c('adult', 'fetal'), each = 4), rep(c('gene', 'exon', 'jxn', 'tx'), 2), SIMPLIFY = FALSE)
names(raw) <- paste0(rep(c('adult', 'fetal'), each = 4), "_", rep(c('gene', 'exon', 'jxn', 'tx'), 2))
top <- lapply(raw, '[[', 'top')
fit <- lapply(raw, '[[', 'fit')
exprsNorm <- lapply(raw, '[[', 'exprsNorm')

## Load BrainSpan model results
raw_span <- mapply(function(age, type) {
    load(paste0('rda/span_limma_region_specific_', age, '_', type, '.Rdata'))
    top$age <- age
    top$type <- type
    return(list(top = top, fit = fit, exprsNorm = exprsNorm))
}, rep(c('adult', 'fetal'), each = 4), rep(c('gene', 'exon', 'jxn', 'tx'), 2), SIMPLIFY = FALSE)
names(raw_span) <- paste0(rep(c('adult', 'fetal'), each = 4), "_", rep(c('gene', 'exon', 'jxn', 'tx'), 2))
top_span <- lapply(raw_span, '[[', 'top')
fit_span <- lapply(raw_span, '[[', 'fit')
exprsNorm_span <- lapply(raw_span, '[[', 'exprsNorm')

get_pcheck <- function(top_table) {
    pcheck <- do.call(rbind, lapply(top_table, function(x) {
        x$P.Bonf <- p.adjust(x$P.Value, 'bonferroni')
        return(x)
    }))
    pcheck$global_fdr <- p.adjust(pcheck$P.Value, 'fdr')
    pcheck$global_bonf <- p.adjust(pcheck$P.Value, 'bonferroni')
    return(pcheck)
}

pcheck <- get_pcheck(top)
pcheck_span <- get_pcheck(top_span)
pcheck_span_tmp <- pcheck_span[, -which(colnames(pcheck_span) %in% c('global_fdr', 'global_bonf'))]
colnames(pcheck_span_tmp) <- paste0('span_', colnames(pcheck_span_tmp))
pcheck_both <- cbind(pcheck, pcheck_span_tmp)
rm(pcheck_span_tmp)





p_summary <- function(pvar = 'FDR', cut = 0.05, pchk) {
    top_table <- split(pchk, paste0(pchk$age, '_', pchk$type))

    num <- sapply(top_table, function(x) {
        if (pvar == 'FDR') {
            table(factor(x$adj.P.Val < cut, levels = c('FALSE', 'TRUE')))
        } else if (pvar == 'bonf') {
            table(factor(x$P.Bonf < cut, levels = c('FALSE', 'TRUE')))
        } else if (pvar == 'global_bonf') {
            table(factor(x$global_bonf < cut, levels = c('FALSE', 'TRUE')))
        } else if (pvar == 'global_fdr') {
            table(factor(x$global_fdr < cut, levels = c('FALSE', 'TRUE')))
        }

    })
    perc <- sweep(num, 2, colSums(num), function(x, y) { round(x / y * 100, 2)} )
    res <- rbind(num, perc)
    status <- as.logical(rownames(res))
    rownames(res) <- NULL
    res <- as.data.frame(res)
    res$DEstatus <- status
    res$method <- pvar
    res$cutoff <- cut
    res$unit <- rep(c('number', 'percent'), each = 2)
    return(res)
}

p_summ_run <- function(pchk) {
    do.call(rbind, mapply(p_summary, rep(c('FDR', 'bonf', 'global_bonf'),
        each = 2), rep(c(0.05, 0.01), 3), MoreArgs = list(pchk = pchk),
        SIMPLIFY = FALSE, USE.NAMES = FALSE))
}


p_sum <- p_summ_run(pcheck)
options(width = 140)
p_sum

p_sum_span <- p_summ_run(pcheck_span)
p_sum_span


summary(pcheck_both$logFC)
summary(pcheck_both$span_logFC)
table('brainseq' = abs(pcheck_both$logFC) > 10, 'span' = abs(pcheck_both$span_logFC) > 10)
table('brainseq' = abs(pcheck_both$logFC) > 20, 'span' = abs(pcheck_both$span_logFC) > 20)
table('brainseq' = abs(pcheck_both$logFC) > 100, 'span' = abs(pcheck_both$span_logFC) > 100)
table('brainseq' = abs(pcheck_both$logFC) > 1000, 'span' = abs(pcheck_both$span_logFC) > 1000)

weird <- which(abs(pcheck_both$logFC) > 10 | abs(pcheck_both$span_logFC) > 10)
table(pcheck_both$type[weird])

summary(c(pcheck_both$logFC[pcheck_both$type == 'jxn'], pcheck_both$span_logFC[pcheck_both$type == 'jxn']))
max(abs(c(pcheck_both$logFC[pcheck_both$type == 'jxn'], pcheck_both$span_logFC[pcheck_both$type == 'jxn'])))

table(pcheck_both$type[which(abs(pcheck_both$logFC) > 16 | abs(pcheck_both$span_logFC) > 16)])



# pdf('pdf/compare_with_span_logFC.pdf')
# ggplot(pcheck_both, aes(x = logFC, y = span_logFC, alpha = 1/10)) +
#     facet_grid(age ~ type) + ylab('BrainSpan log FC') +
#     xlab('BrainSeq log FC') + geom_point() + xlim(-16, 16) + ylim(-16, 16) +
#     geom_smooth(method=lm, se=FALSE)
# dev.off()

png('pdf/compare_with_span_logFC.png', type = 'cairo'
)
ggplot(pcheck_both, aes(x = logFC, y = span_logFC, alpha = 1/20)) +
    facet_grid(age ~ type) + ylab('BrainSpan log FC') +
    xlab('BrainSeq log FC') + geom_point() + xlim(-16, 16) + ylim(-16, 16) +
    geom_smooth(method=lm, se=FALSE)
dev.off()

png('pdf/compare_with_span_logFC_noTx.png', type = 'cairo')
ggplot(subset(pcheck_both, type != 'tx'), aes(x = logFC, y = span_logFC,
    alpha = 1/20)) +
    facet_grid(age ~ type) + ylab('BrainSpan log FC') +
    xlab('BrainSeq log FC') + geom_point() + xlim(-7.5, 7.5) + ylim(-7.5, 7.5) +
    geom_smooth(method=lm, se=FALSE)
dev.off()

# pdf('pdf/compare_with_span_logFC_noTx.pdf')
# ggplot(subset(pcheck_both, type != 'tx'), aes(x = logFC, y = span_logFC,
#     alpha = 1/10)) +
#     facet_grid(age ~ type) + ylab('BrainSpan log FC') +
#     xlab('BrainSeq log FC') + geom_point() + xlim(-7.5, 7.5) + ylim(-7.5, 7.5) +
#     geom_smooth(method=lm, se=FALSE)
# dev.off()


pdf('pdf/compare_with_span_logFC_density.pdf')
ggplot(pcheck_both, aes(x = logFC, y = span_logFC)) +
    facet_grid(age ~ type) + ylab('BrainSpan log FC') +
    xlab('BrainSeq log FC') + stat_density2d() + xlim(-16, 16) + ylim(-16, 16) +
    geom_smooth(method=lm, se=FALSE)
dev.off()

pdf('pdf/compare_with_span_logFC_density_noTx.pdf')
ggplot(subset(pcheck_both, type != 'tx'), aes(x = logFC, y = span_logFC)) +
    facet_grid(age ~ type) + ylab('BrainSpan log FC') +
    xlab('BrainSeq log FC') + stat_density2d() + xlim(-7.5, 7.5) +
    ylim(-7.5, 7.5) +
    geom_smooth(method=lm, se=FALSE)
dev.off()


## Replication (p < 0.05) & same direction

rep_span <- do.call(rbind, mapply(function(pvar, cut, type_sub, age_sub) {
    pinfo <- subset(pcheck_both, type == type_sub & age == age_sub)
    n <- sum(sign(pinfo$logFC) == sign(pinfo$span_logFC) & pinfo$span_P.Value < 0.05 & pinfo[, pvar] < cut)
    n_sign <- sum(sign(pinfo$logFC) == sign(pinfo$span_logFC) & pinfo[, pvar] < cut)
    data.frame(pvar = pvar, cutoff = cut, type = type_sub, age = age_sub, replicated = n, number_de = sum(pinfo[, pvar] < cut), total = nrow(pinfo), replicated_sign = n_sign, stringsAsFactors = FALSE)
}, pvar = rep(rep(c('adj.P.Val', 'P.Bonf'), each = 6), 4 * 2), cut = rep(c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001), 2 * 4 * 2), type_sub = rep(rep(unique(pcheck_both$type), each = 6 * 2), 2), age_sub = rep(unique(pcheck_both$age), each = 4 * 6 * 2), SIMPLIFY = FALSE, USE.NAMES = FALSE))

pdf('pdf/replication_exploration.pdf', width = 14)
ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated / number_de, color = pvar)) + facet_grid(age ~ type) + ylab('Replication rate') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(color='P-value method') + ylim(c(0, 1))

ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = number_de, color = pvar)) + facet_grid(age ~ type) + ylab('Number of DE features') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10() + labs(color='P-value method')

ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = number_de / total * 100, color = pvar)) + facet_grid(age ~ type) + ylab('Percent of DE features') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(color='P-value method') + ylim(c(0, 100))

ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated_sign / number_de, color = pvar)) + facet_grid(age ~ type) + ylab('Replication rate (sign only)') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(color='P-value method') + ylim(c(0, 1))
dev.off()



## Explore some of the adjustment variables
head(fit$fetal_gene$design)
head(fit_span$fetal_gene$design)
nrow(fit$fetal_gene$design)
nrow(fit_span$fetal_gene$design)
summary(fit$fetal_gene$design)
summary(fit_span$fetal_gene$design)
summary(fit$fetal_gene$design[, 'mean_RIN'])
summary(fit_span$fetal_gene$design[, 'mean_RIN'])
sort(fit$fetal_gene$design[, 'mean_RIN'])
sort(fit_span$fetal_gene$design[, 'mean_RIN'])

summary(fit$fetal_gene$design[, 'Age'])
summary(fit_span$fetal_gene$design[, 'Age'])
t.test(fit$fetal_gene$design[, 'Age'], fit_span$fetal_gene$design[, 'Age'])
t.test(fit$fetal_gene$design[, 'mean_mitoRate'], fit_span$fetal_gene$design[, 'mean_mitoRate'])
t.test(fit$fetal_gene$design[, 'mean_totalAssignedGene'], fit_span$fetal_gene$design[, 'mean_totalAssignedGene']) ## Diff
t.test(fit$fetal_gene$design[, 'mean_RIN'], fit_span$fetal_gene$design[, 'mean_RIN'])

plot(fit$fetal_gene$design[, 'snpPC1'], fit$fetal_gene$design[, 'snpPC2'])
plot(fit_span$fetal_gene$design[, 'snpPC1'], fit_span$fetal_gene$design[, 'snpPC2'])

t.test(fit$fetal_gene$design[fit$fetal_gene$design[, 'RegionHIPPO'] == 1, 'mean_totalAssignedGene'], fit_span$fetal_gene$design[fit_span$fetal_gene$design[, 'RegionHIPPO'] == 1, 'mean_totalAssignedGene'])
t.test(fit$fetal_gene$design[fit$fetal_gene$design[, 'RegionHIPPO'] == 0, 'mean_totalAssignedGene'], fit_span$fetal_gene$design[fit_span$fetal_gene$design[, 'RegionHIPPO'] == 0, 'mean_totalAssignedGene'])
t.test(fit$fetal_gene$design[, 'mean_totalAssignedGene'] ~ fit$fetal_gene$design[, 'RegionHIPPO'])
t.test(fit_span$fetal_gene$design[, 'mean_totalAssignedGene'] ~ fit_span$fetal_gene$design[, 'RegionHIPPO'])



plot(-log10(pcheck$global_fdr), -log10(pcheck$adj.P.Val), col = c('gene' = 'blue', 'exon' = 'orange', 'jxn' = 'grey20', 'tx' = 'light blue')[pcheck$type], pch = c('adult' = 21, 'fetal' = 22)[pcheck$age])
abline(a = 0, b = 1, col = 'red')

table('Global FDR' = pcheck$global_fdr < 0.05, 'FDR' = pcheck$adj.P.Val < 0.05, 'Age group' = pcheck$age, 'Feature type' = pcheck$type)
table('Global FDR' = pcheck$global_fdr < 0.01, 'FDR' = pcheck$adj.P.Val < 0.01, 'Age group' = pcheck$age, 'Feature type' = pcheck$type)

table('Global Bonf' = pcheck$global_bonf < 0.05, 'Bonf' = pcheck$P.Bonf < 0.05, 'Age group' = pcheck$age, 'Feature type' = pcheck$type)
table('Global Bonf' = pcheck$global_bonf < 0.01, 'Bonf' = pcheck$P.Bonf < 0.01, 'Age group' = pcheck$age, 'Feature type' = pcheck$type)


## Numbers for BrainSeq phase 2 update
pinfo <- lapply(unique(pcheck_both$type), function(feat)  {
    res <- lapply(unique(pcheck_both$age), function(agegrp) {
        pinfo <- subset(pcheck_both, type == feat & age == agegrp)
        pinfo[sign(pinfo$t) == sign(pinfo$span_t) & pinfo$span_P.Value < 0.05 & pinfo$P.Bonf < 0.01, ]
    })
    names(res) <- unique(pcheck_both$age)
    return(res)
})
names(pinfo) <- unique(pcheck_both$type)
sapply(pinfo, function(x) sapply(x, nrow))
#       gene  exon  jxn   tx
# adult 1612 15442 5561 1750
# fetal   29    65   16    0

rses <- lapply(unique(pcheck_both$type), function(type) {
    load_foo(type, 'fetal')
})
names(rses) <- unique(pcheck_both$type)

## If all of them are needed
# rses <- lapply(unique(pcheck_both$type), function(type) {
#     res <- lapply(unique(pcheck_both$age), function(agegrp) load_foo(type, agegrp))
#     names(res) <- unique(pcheck_both$age)
#     return(res)
# })


de_genes <- lapply(names(pinfo[[1]]), function(agegrp) {
     res2 <- lapply(names(pinfo), function(feat) {
        m <- match(gsub('.*gene.|.*exon.|.*jxn.|.*tx.', '', rownames(pinfo[[feat]][[agegrp]])), names(rses[[feat]]))
        print(table(!is.na(m)))
        if(feat %in% c('gene', 'exon')) {
            res <- rowRanges(rses[[feat]])$gencodeID[m]
        } else if(feat == 'jxn') {
            res <- rowRanges(rses[[feat]])$gencodeGeneID[m]
            res <- res[!is.na(res)]
        } else {
            res <- rowRanges(rses[[feat]])$gene_id[m]
        }
        return(unique(res))
    })
    names(res2) <- names(pinfo)
    return(res2)
})
names(de_genes) <- names(pinfo$gene)
sapply(de_genes, function(x) sapply(x, length))
#      adult fetal
# gene  1612    29
# exon  2686    11
# jxn   1897     7
# tx    1438     0

library(gplots)

pdf('pdf/venn_de_features.pdf')
venn(de_genes$adult) + title('DE features grouped by gene id (adult)')
venn(de_genes$adult[c('gene', 'exon', 'jxn')]) + title('DE features grouped by gene id (adult)')
venn(de_genes$fetal) + title('DE features grouped by gene id (fetal)')
venn(de_genes$fetal[c('gene', 'exon', 'jxn')]) + title('DE features grouped by gene id (fetal)')
dev.off()



##
rse_gene <- load_foo('gene', 'fetal')


corfit <- duplicateCorrelation(exprsNorm$fetal_gene, fit$fetal_gene$design[, c('(Intercept)',
       'RegionHIPPO')], block=colData(rse_gene)$BrNum)
print('Consensus correlation and summary (also after tanh transform)')
corfit$consensus.correlation
summary(corfit$atanh.correlations)
summary(tanh(corfit$atanh.correlations))

rse_tx <- load_foo('tx', 'adult')
rse_jxn <- load_foo('jxn', 'adult')

rse_gene_span <- load_span('gene', 'adult')
rse_tx_span <- load_span('tx', 'adult')
rse_jxn_span <- load_span('jxn', 'adult')




## Reg specific model
design <- get_mods( colData(rse) )$mod

min(top$adult_gene$adj.P.Val)
which.min(top$adult_gene$adj.P.Val)
which(rank(top$adult_gene$adj.P.Val) == 1)




cleanedVoom <- cleaningY(exprsNorm$adult_gene, design, 2)


top$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 1), ]
boxplot(exprsNorm$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 1), ] ~ colData(rse)$Region)
t.test(exprsNorm$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 1), ] ~ colData(rse)$Region)

boxplot(cleanedVoom[which(rank(top$adult_gene$adj.P.Val) == 1), ] ~ colData(rse)$Region)
t.test(cleanedVoom[which(rank(top$adult_gene$adj.P.Val) == 1), ] ~ colData(rse)$Region)

top$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 2), ]
boxplot(exprsNorm$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 2), ] ~ colData(rse)$Region)
t.test(exprsNorm$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 2), ] ~ colData(rse)$Region)

boxplot(cleanedVoom[which(rank(top$adult_gene$adj.P.Val) == 2), ] ~ colData(rse)$Region)
t.test(cleanedVoom[which(rank(top$adult_gene$adj.P.Val) == 2), ] ~ colData(rse)$Region)

boxplot(assays(rse)$rpkm[which(rank(top$adult_gene$adj.P.Val) == 2), ] ~ colData(rse)$Region)



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
