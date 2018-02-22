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

load_span <- function(type) {
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
pinfo <- pinfo[sign(pinfo$F) == sign(pinfo$span_F) & pinfo$span_P.Value < 0.05 & pinfo$P.Bonf < 0.01, ]


tocheck <- match(gsub('gene.', '', head(rownames(pinfo), 100)), rownames(rse))
tocheck


## Re-organize design before cleaningY
protect <- grepl(':RegionHIPPO|RegionHIPPO:|Intercept', colnames(design))
design <- cbind(design[, protect], design[, !protect])

cleanedVoom <- cleaningY(exprsNorm$gene, design, sum(protect))

i <- 1
top$gene[tocheck[i], ]

get_main <- function(i) {
    j <- tocheck[i]
    paste(with(rowRanges(rse)[j, ], paste(gencodeID, Symbol)), 'p-bonf', signif(pinfo$P.Bonf[i], 3))
}

plot_age_mod <-
    design[, c(
        '(Intercept)',
        'Age',
        'RegionHIPPO',
        'fetal',
        'birth',
        'infant',
        'child',
        'teen',
        'adult'
    )]
p_cols <- ifelse(colData(rse)$Region == 'HIPPO', 'blue', 'orange')
l_cols <- c('lightgoldenrod', 'light blue')
age_brks <- c(-1, 0, 1, 10, 20, 50, 100)

pdf('pdf/top_gene_replicated_exprNorm.pdf', width = 14)
for(i in 1:100) {
    agePlotter(
        exprsNorm$gene[tocheck[i],],
        colData(rse)$Age,
        pointColor = p_cols,
        ageBreaks = age_brks,
        mainText = get_main(i),
        lineColor = l_cols,
        mod = plot_age_mod,
        ylab = 'Voom-normalized expression'
    )
    legend('bottom', c('DLPFC', 'HIPPO'), col = l_cols, lwd = 3, bty = 'n', ncol = 2)
}
dev.off()

pdf('pdf/top_gene_replicated_cleanedVoom.pdf', width = 14)
for(i in 1:100) {
    agePlotter(
        cleanedVoom[tocheck[i],],
        colData(rse)$Age,
        pointColor = p_cols,
        ageBreaks = age_brks,
        mainText = get_main(i),
        lineColor = l_cols,
        mod = plot_age_mod,
        ylab = 'Voom-normalized expression with covariate effects removed'
    )
    legend('bottom', c('DLPFC', 'HIPPO'), col = l_cols, lwd = 3, bty = 'n', ncol = 2)
}
dev.off()


## Numbers for BOG2018 abstract
pinfo <- lapply(unique(pcheck_both$type), function(feat)  {
    pinfo <- subset(pcheck_both, type == feat)
    pinfo[sign(pinfo$F) == sign(pinfo$span_F) & pinfo$span_P.Value < 0.05 & pinfo$P.Bonf < 0.01, ]
})
names(pinfo) <- unique(pcheck_both$type)


rses <- lapply(unique(pcheck_both$type), load_foo)
names(rses) <- unique(pcheck_both$type)

sapply(pinfo, nrow)

de_genes <- lapply(names(pinfo), function(feat) {
    m <- match(gsub('gene.|exon.|jxn.|tx.', '', rownames(pinfo[[feat]])), names(rses[[feat]]))
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
names(de_genes) <- names(pinfo)
sapply(de_genes, length)


library(gplots)

pdf('pdf/venn_de_features.pdf')
venn(de_genes) + title('DE features grouped by gene id')
venn(de_genes[c('gene', 'exon', 'jxn')]) + title('DE features grouped by gene id')
dev.off()

load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/caseControl/dxStats_hippo_filtered_qSVA.rda')

outInfo <- list('gene' = outGene, 'exon' = outExon, 'jxn' = outJxn, 'tx' = outTx)
case_genes <- lapply(names(pinfo), function(feat) {
    m <- match(rownames(outInfo[[feat]][outInfo[[feat]]$adj.P.Val < 0.05, ]), names(rses[[feat]]))
    print(table(!is.na(m)))
    m <- m[!is.na(m)]
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
names(case_genes) <- names(pinfo)
sapply(case_genes, length)

pdf('pdf/venn_de_features_caseControl.pdf')
venn(case_genes) + title('DE features grouped by gene id')
venn(case_genes[c('gene', 'exon', 'jxn')]) + title('DE features grouped by gene id')
dev.off()


load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eqtl_tables/matrixEqtl_output_interaction_4features.rda', verbose = TRUE)

me <- list('gene'= meGene$cis$eqtls, 'exon' = meExon$cis$eqtls, 'jxn' = meJxn$cis$eqtls, 'tx' = meTx$cis$eqtls)

me_genes <- lapply(names(pinfo), function(feat) {
    m <- match(me[[feat]]$gene[me[[feat]]$FDR < 0.01], names(rses[[feat]]))
    print(table(!is.na(m)))
    m <- m[!is.na(m)]
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
names(me_genes) <- names(pinfo)
sapply(me_genes, length)

pdf('pdf/venn_eQTL_interaction.pdf')
venn(me_genes) + title('eQTLs grouped by gene id')
venn(me_genes[c('gene', 'exon', 'jxn')]) + title('eQTLs grouped by gene id')
dev.off()

length(unique(unlist(me_genes[c('gene', 'exon', 'jxn')])))


riskLoci = read.csv("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/PGC_risk_loci.csv", stringsAsFactors=FALSE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda", verbose = TRUE)
keepIndex = which(!is.na(snpMap$chr_hg38))
snpMap = snpMap[keepIndex,]
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

in_risk <- lapply(me, function(eqtl) {
    m <- match(eqtl$snps, rownames(snpMap))
    print(table(!is.na(m)))
    m <- m[!is.na(m)]
    snpMap$pos_hg19[m] %in% riskLoci$hg19POS
    res <- snpMap$pos_hg19[m] %in% riskLoci$hg19POS
    names(res) <- eqtl$gene
    return(res)
})
sapply(in_risk, sum)
sapply(in_risk, mean) * 100

in_riskFDR <- lapply(me, function(eqtl) {
    if(length(eqtl$snps[eqtl$FDR < 0.01]) == 0) return(NULL)
    m <- match(eqtl$snps[eqtl$FDR < 0.01], rownames(snpMap))
    print(table(!is.na(m)))
    m <- m[!is.na(m)]
    res <- snpMap$pos_hg19[m] %in% riskLoci$hg19POS
    names(res) <- eqtl$gene[eqtl$FDR < 0.01]
    return(res)
})
sapply(in_riskFDR, sum)
sapply(in_riskFDR, mean) * 100

qtls <- lapply(me, function(eqtl) {
    eqtl$snps[eqtl$FDR < 0.01]
})
venn(qtls)
venn(qtls[c('gene', 'exon', 'jxn')])
qtl_v <- venn(qtls[c('gene', 'exon', 'jxn')], show.plot = FALSE)
sum(sapply(attr(qtl_v, 'intersections'), length))




sapply(in_riskFDR, sum) / sapply(in_risk, sum) / ( sapply(in_riskFDR, length) /  sapply(in_risk, length))


v_me <- venn(me_genes[c('gene', 'exon', 'jxn')], show.plot = FALSE)
risk3 <- lapply(list(all = in_risk$gene, FDR = in_riskFDR$gene), function(rk) {
    names(rk[rk]) %in% c(attr(v_me, 'intersections')[['gene:exon:jxn']], attr(v_me, 'intersections')[['jxn']])
} )
sapply(risk3, sum)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


