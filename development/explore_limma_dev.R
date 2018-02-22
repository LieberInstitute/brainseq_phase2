library('limma')
library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('ggplot2')

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
fm_mod <-  ~Age + fetal + birth + infant + child + teen + adult + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN
fm_mod0 <- ~ Age + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN
fm_mod_all <- ~Age *Region + fetal * Region + birth *Region + infant *Region + child * Region + teen * Region + adult * Region + Sex + Region + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN
fm_mod0_all <- ~ Age + fetal + birth + infant + child + teen + adult + Sex + Region + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 + mean_mitoRate + mean_totalAssignedGene + mean_RIN


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


## Load BrainSeq model results
raw <- lapply(c('gene', 'exon', 'jxn', 'tx'), function(type) {
    load(paste0('rda/limma_dev_interaction_', type, '.Rdata'))
    top$type <- type
    return(list(top = top, fit = fit, exprsNorm = exprsNorm))
})
names(raw) <- c('gene', 'exon', 'jxn', 'tx')
top <- lapply(raw, '[[', 'top')
fit <- lapply(raw, '[[', 'fit')
exprsNorm <- lapply(raw, '[[', 'exprsNorm')

## Load BrainSpan model results
raw_span <- lapply(c('gene', 'exon', 'jxn', 'tx'), function(type) {
    load(paste0('rda/span_limma_dev_interaction_', type, '.Rdata'))
    top$type <- type
    return(list(top = top, fit = fit, exprsNorm = exprsNorm))
})
names(raw_span) <- c('gene', 'exon', 'jxn', 'tx')
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
    top_table <- split(pchk, pchk$type)

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
    do.call(rbind, mapply(p_summary, rep(c('FDR', 'bonf', 'global_bonf'), each = 2),
        rep(c(0.05, 0.01), 3), MoreArgs = list(pchk = pchk),
        SIMPLIFY = FALSE, USE.NAMES = FALSE))
}


p_sum <- p_summ_run(pcheck)
options(width = 140)
p_sum

p_sum_span <- p_summ_run(pcheck_span)
p_sum_span


summary(pcheck_both$Age.RegionHIPPO)
summary(pcheck_both$span_Age.RegionHIPPO)
table('brainseq' = abs(pcheck_both$Age.RegionHIPPO) > 10, 'span' = abs(pcheck_both$span_Age.RegionHIPPO) > 10)
table('brainseq' = abs(pcheck_both$Age.RegionHIPPO) > 20, 'span' = abs(pcheck_both$span_Age.RegionHIPPO) > 20)
table('brainseq' = abs(pcheck_both$Age.RegionHIPPO) > 100, 'span' = abs(pcheck_both$span_Age.RegionHIPPO) > 100)
table('brainseq' = abs(pcheck_both$Age.RegionHIPPO) > 1000, 'span' = abs(pcheck_both$span_Age.RegionHIPPO) > 1000)

weird <- which(abs(pcheck_both$Age.RegionHIPPO) > 10 | abs(pcheck_both$span_Age.RegionHIPPO) > 10)
table(pcheck_both$type[weird])

summary(c(pcheck_both$Age.RegionHIPPO[pcheck_both$type == 'jxn'], pcheck_both$span_Age.RegionHIPPO[pcheck_both$type == 'jxn']))
max(abs(c(pcheck_both$Age.RegionHIPPO[pcheck_both$type == 'jxn'], pcheck_both$span_Age.RegionHIPPO[pcheck_both$type == 'jxn'])))

table(pcheck_both$type[which(abs(pcheck_both$Age.RegionHIPPO) > 16 | abs(pcheck_both$span_Age.RegionHIPPO) > 16)])

png('pdf/compare_with_span_logFC.png', type = 'cairo'
)
ggplot(pcheck_both, aes(x = F, y = span_F, alpha = 1/20)) +
    facet_grid(. ~ type) + ylab('BrainSpan F-stat') +
    xlab('BrainSeq F-stat') + geom_point()
    geom_smooth(method=lm, se=FALSE)
dev.off()

png('pdf/compare_with_span_logFC_noTx.png', type = 'cairo')
ggplot(subset(pcheck_both, type != 'tx'), aes(x = logFC, y = span_logFC,
                                              alpha = 1/20)) +
    facet_grid(age ~ type) + ylab('BrainSpan log FC') +
    xlab('BrainSeq log FC') + geom_point() + xlim(-7.5, 7.5) + ylim(-7.5, 7.5) +
    geom_smooth(method=lm, se=FALSE)
dev.off()


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
rep_span <- do.call(rbind, mapply(function(pvar, cut, type_sub) {
    pinfo <- subset(pcheck_both, type == type_sub)
    vars <- gsub('span_', '', colnames(pcheck_both)[grep('Region', colnames(pcheck_both))])
    vars <- vars[which(duplicated(vars))]
    do.call(rbind, lapply(c('F', vars), function(var) {
        n <- sum(sign(pinfo[, var]) == sign(pinfo[, paste0('span_', var)]) & pinfo$span_P.Value < 0.05 & pinfo[, pvar] < cut)
        n_sign <- sum(sign(pinfo[, var]) == sign(pinfo[, paste0('span_', var)]) & pinfo[, pvar] < cut)
        data.frame(pvar = pvar, cutoff = cut, type = type_sub, replicated = n, number_de = sum(pinfo[, pvar] < cut), total = nrow(pinfo), replicated_sign = n_sign, term = var, stringsAsFactors = FALSE)
    }))
}, pvar = rep(rep(c('adj.P.Val', 'P.Bonf'), each = 6), 4), cut = rep(c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001), 2 * 4), type_sub = rep(unique(pcheck_both$type), each = 6 * 2), SIMPLIFY = FALSE, USE.NAMES = FALSE))
rep_span$term_clean <- tolower(gsub('\\.|RegionHIPPO', '', rep_span$term))

pdf('pdf/replication_exploration.pdf', width = 14, height =  14)
ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated / number_de, color = pvar)) + facet_grid(term_clean ~ type) + ylab('Replication rate') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(color='P-value method') + ylim(c(0, 1))

ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = number_de, color = pvar)) + facet_grid(term_clean ~ type) + ylab('Number of DE features') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10() + labs(color='P-value method')

ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = number_de / total * 100, color = pvar)) + facet_grid(term_clean ~ type) + ylab('Percent of DE features') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(color='P-value method') + ylim(c(0, 100))


ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated_sign / number_de, color = pvar)) + facet_grid(term_clean ~ type) + ylab('Replication rate (sign only)') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(color='P-value method') + ylim(c(0, 1))
dev.off()


## Explore some of the adjustment variables
head(fit$gene$design)
head(fit_span$gene$design)
nrow(fit$gene$design)
nrow(fit_span$gene$design)
summary(fit$gene$design)
summary(fit_span$gene$design)
summary(fit$gene$design[, 'mean_RIN'])
summary(fit_span$gene$design[, 'mean_RIN'])
sort(fit$gene$design[, 'mean_RIN'])
sort(fit_span$gene$design[, 'mean_RIN'])

summary(fit$gene$design[, 'Age'])
summary(fit_span$gene$design[, 'Age'])
t.test(fit$gene$design[, 'Age'], fit_span$gene$design[, 'Age'])
t.test(fit$gene$design[, 'mean_mitoRate'], fit_span$gene$design[, 'mean_mitoRate'])
t.test(fit$gene$design[, 'mean_totalAssignedGene'], fit_span$gene$design[, 'mean_totalAssignedGene']) ## Diff
t.test(fit$gene$design[, 'mean_RIN'], fit_span$gene$design[, 'mean_RIN'])

plot(fit$gene$design[, 'snpPC1'], fit$gene$design[, 'snpPC2'])
plot(fit_span$gene$design[, 'snpPC1'], fit_span$gene$design[, 'snpPC2'])

t.test(fit$gene$design[fit$gene$design[, 'RegionHIPPO'] == 1, 'mean_totalAssignedGene'], fit_span$gene$design[fit_span$gene$design[, 'RegionHIPPO'] == 1, 'mean_totalAssignedGene'])
t.test(fit$gene$design[fit$gene$design[, 'RegionHIPPO'] == 0, 'mean_totalAssignedGene'], fit_span$gene$design[fit_span$gene$design[, 'RegionHIPPO'] == 0, 'mean_totalAssignedGene'])
t.test(fit$gene$design[, 'mean_totalAssignedGene'] ~ fit$gene$design[, 'RegionHIPPO'])
t.test(fit_span$gene$design[, 'mean_totalAssignedGene'] ~ fit_span$gene$design[, 'RegionHIPPO'])



plot(-log10(pcheck$global_fdr), -log10(pcheck$adj.P.Val), col = c('gene' = 'blue', 'exon' = 'orange', 'jxn' = 'grey20', 'tx' = 'light blue')[pcheck$type], pch = c('adult' = 21, 'fetal' = 22)[pcheck$age])
abline(a = 0, b = 1, col = 'red')

table('Global FDR' = pcheck$global_fdr < 0.05, 'FDR' = pcheck$adj.P.Val < 0.05, 'Feature type' = pcheck$type)
table('Global FDR' = pcheck$global_fdr < 0.01, 'FDR' = pcheck$adj.P.Val < 0.01, 'Feature type' = pcheck$type)

table('Global Bonf' = pcheck$global_bonf < 0.05, 'Bonf' = pcheck$P.Bonf < 0.05, 'Feature type' = pcheck$type)
table('Global Bonf' = pcheck$global_bonf < 0.01, 'Bonf' = pcheck$P.Bonf < 0.01,  'Feature type' = pcheck$type)








rse <- load_foo('gene')


# corfit <- duplicateCorrelation(exprsNorm$fetal_gene, fit$fetal_gene$design[, c('(Intercept)',
#                                                                                'RegionHIPPO')], block=colData(rse_gene)$BrNum)
# print('Consensus correlation and summary (also after tanh transform)')
# corfit$consensus.correlation
# summary(corfit$atanh.correlations)
# summary(tanh(corfit$atanh.correlations))
#
# rse_tx <- load_foo('tx')
# rse_jxn <- load_foo('jxn')
#
# rse_gene_span <- load_span('gene')
# rse_tx_span <- load_span('tx')
# rse_jxn_span <- load_span('jxn')



## Reg specific model
design <- get_mods( colData(rse), int = TRUE)$mod

min(top$gene$adj.P.Val)
which.min(top$gene$adj.P.Val)
which(rank(top$gene$adj.P.Val) == 1)


pinfo <- subset(pcheck_both, type == 'gene')
pinfo <- pinfo[order(pinfo$P.Bonf), ]
pinfo <- pinfo[sign(pinfo$F) == sign(pinfo$span_F) & pinfo$span_P.Value < 0.05 & pinfo$P.Bonf < 0.000001, ]


tocheck <- match(gsub('gene.', '', head(rownames(pinfo))), rownames(rse))
tocheck


## Re-organize design before cleaningY


cleanedVoom <- cleaningY(exprsNorm$gene, design, 2)

i <- 1

top$gene[tocheck[i], ]


plot(exprsNorm$gene[tocheck[i], ] ~ colData(rse)$Age, col = ifelse(colData(rse)$Region == 'HIPPO', 'blue', 'orange'))
t.test(exprsNorm$gene[tocheck[i], ] ~ colData(rse)$Region)

boxplot(cleanedVoom[tocheck[i], ] ~ colData(rse)$Region)
t.test(cleanedVoom[tocheck[i], ] ~ colData(rse)$Region)

boxplot(assays(rse)$rpkm[which(rank(top$gene$adj.P.Val) == 2), ] ~ colData(rse)$Region)


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


