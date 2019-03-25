library('limma')
library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('ggplot2')
library('gplots')
library('VennDiagram')
library('RColorBrewer')
library('clusterProfiler')

source('load_funs.R')
dir.create('rda', showWarnings = FALSE)
dir.create('pdf', showWarnings = FALSE)

## Load BrainSeq model results
if(!file.exists('rda/raw.Rdata')) {
    raw <- lapply(c('gene', 'exon', 'jxn', 'tx'), function(type) {
        f <- paste0('rda/limma_dev_interaction_', type, '.Rdata')
        message(paste(Sys.time(), 'loading', f))
        load(f, verbose = TRUE)
        top$type <- type
        return(list(top = top, fit = fit, exprsNorm = exprsNorm))
    })
    names(raw) <- c('gene', 'exon', 'jxn', 'tx')
    message(paste(Sys.time(), 'saving rda/raw.Rdata'))
    save(raw, file = 'rda/raw.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/raw.Rdata'))
    load('rda/raw.Rdata', verbose = TRUE)
}

top <- lapply(raw, '[[', 'top')
fit <- lapply(raw, '[[', 'fit')
exprsNorm <- lapply(raw, '[[', 'exprsNorm')

## Load BrainSpan model results
if(!file.exists('rda/raw_span.Rdata')) {
    raw_span <- lapply(c('gene', 'exon', 'jxn', 'tx'), function(type) {
        f <- paste0('rda/span_limma_dev_interaction_', type, '.Rdata')
        message(paste(Sys.time(), 'loading', f))
        load(f, verbose = TRUE)
        top$type <- type
        return(list(top = top, fit = fit, exprsNorm = exprsNorm))
    })
    names(raw_span) <- c('gene', 'exon', 'jxn', 'tx')
    message(paste(Sys.time(), 'saving rda/raw_span.Rdata'))
    save(raw_span, file = 'rda/raw_span.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/raw_span.Rdata'))
    load('rda/raw_span.Rdata', verbose = TRUE)
}

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

if(!file.exists('rda/pcheck_both.Rdata')) {
    pcheck <- get_pcheck(top)
    pcheck_span <- get_pcheck(top_span)
    pcheck_span_tmp <- pcheck_span[, -which(colnames(pcheck_span) %in% c('global_fdr', 'global_bonf'))]
    colnames(pcheck_span_tmp) <- paste0('span_', colnames(pcheck_span_tmp))
    pcheck_both <- cbind(pcheck, pcheck_span_tmp)
    rm(pcheck_span_tmp)
    save(pcheck, pcheck_span, pcheck_both, file = 'rda/pcheck_both.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/pcheck_both.Rdata'))
    load('rda/pcheck_both.Rdata', verbose = TRUE)
}

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

if(!file.exists('rda/p_sum.Rdata')) {
    p_sum <- p_summ_run(pcheck)
    p_sum_span <- p_summ_run(pcheck_span)
    save(p_sum, p_sum_span, file = 'rda/p_sum.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/p_sum.Rdata'))
    load('rda/p_sum.Rdata', verbose = TRUE)
}

options(width = 140)
p_sum
p_sum_span

summary(pcheck_both$Age.RegionHIPPO)
summary(pcheck_both$span_Age.RegionHIPPO)

# pdf('pdf/compare_with_span_F.pdf', useDingbats = FALSE)
png('pdf/compare_with_span_F.png', type = 'cairo')
ggplot(pcheck_both, aes(x = F, y = span_F, alpha = 1/20)) +
    facet_grid(. ~ type, scales = 'free') + ylab('BrainSpan F-stat') +
    xlab('BrainSeq F-stat') + geom_point() +
    geom_smooth(method=lm, se=FALSE) + scale_x_log10() + scale_y_log10()
dev.off()

# pdf('pdf/compare_with_span_F_noTx.pdf', useDingbats = FALSE)
png('pdf/compare_with_span_F_noTx.png', type = 'cairo')
ggplot(subset(pcheck_both, type != 'tx'), aes(x = F, y = span_F,
                                              alpha = 1/20)) +
    facet_grid(. ~ type, scales = 'free') + ylab('BrainSpan F-stat') +
    xlab('BrainSeq F-stat') + geom_point() +
    geom_smooth(method=lm, se=FALSE) + scale_x_log10() + scale_y_log10()
dev.off()


pdf('pdf/compare_with_span_F_density.pdf', useDingbats = FALSE)
ggplot(pcheck_both, aes(x = F, y = span_F)) +
    facet_grid(. ~ type, scales = 'free') + ylab('BrainSpan F-stat') +
    xlab('BrainSeq F-stat') + stat_density2d() +
    geom_smooth(method=lm, se=FALSE) + scale_x_log10() + scale_y_log10()
dev.off()

pdf('pdf/compare_with_span_F_density_noTx.pdf', useDingbats = FALSE)
ggplot(subset(pcheck_both, type != 'tx'), aes(x = F, y = span_F)) +
    facet_grid(. ~ type, scales = 'free') + ylab('BrainSpan F-stat') +
    xlab('BrainSeq F-stat') + stat_density2d() +
    geom_smooth(method=lm, se=FALSE) + scale_x_log10() + scale_y_log10()
dev.off()



## Replication (p < 0.05) & same direction
if(!file.exists('rda/rep_span.Rdata')) {
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
    save(rep_span, file = 'rda/rep_span.Rdata')
} else {
    message(paste(Sys.time(), 'rda/rep_span.Rdata'))
    load('rda/rep_span.Rdata', verbose = TRUE)
}



pdf('pdf/replication_exploration.pdf', width = 14, height =  14, useDingbats = FALSE)
ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated / number_de, color = pvar)) + facet_grid(term_clean ~ type) + ylab('Replication rate') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(color='P-value method') + ylim(c(0, 1))

ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = number_de, color = pvar)) + facet_grid(term_clean ~ type) + ylab('Number of DE features') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10() + labs(color='P-value method')

ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = number_de / total * 100, color = pvar)) + facet_grid(term_clean ~ type) + ylab('Percent of DE features') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(color='P-value method') + ylim(c(0, 100))


ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated_sign / number_de, color = pvar)) + facet_grid(term_clean ~ type) + ylab('Replication rate (sign only)') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(color='P-value method') + ylim(c(0, 1))
dev.off()


pdf('pdf/replication_exploration_subset.pdf', width = 14, height =  10, useDingbats = FALSE)
ggplot(subset(rep_span, pvar == 'P.Bonf' & type != 'tx'), aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated / number_de)) + facet_grid(type ~ toupper(term_clean)) + ylab('Replication rate') + xlab('p-value threshold') + geom_point() + theme_bw(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(c(0, 1))
dev.off()

pdf('pdf/replication_exploration_subset_F.pdf', width = 8, height =  10, useDingbats = FALSE)
ggplot(subset(rep_span, pvar == 'P.Bonf' & type != 'tx' & term_clean == 'f'), aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated / number_de)) + facet_grid(type ~ toupper(term_clean)) + ylab('Replication rate') + xlab('p-value threshold') + geom_point() + theme_bw(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(c(0, 1)) + geom_line(aes(y = replicated / number_de, x = rep(1:6, 3)))

ggplot(subset(rep_span, pvar == 'P.Bonf' & type != 'tx' & term_clean == 'f'), aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated / number_de)) + facet_grid(type ~ toupper(term_clean)) + ylab('Replication rate') + xlab('p-value threshold') + geom_point() + theme_bw(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_line(aes(y = replicated / number_de, x = rep(1:6, 3)))
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

# plot(fit$gene$design[, 'snpPC1'], fit$gene$design[, 'snpPC2'])
# plot(fit_span$gene$design[, 'snpPC1'], fit_span$gene$design[, 'snpPC2'])

t.test(fit$gene$design[fit$gene$design[, 'RegionHIPPO'] == 1, 'mean_totalAssignedGene'], fit_span$gene$design[fit_span$gene$design[, 'RegionHIPPO'] == 1, 'mean_totalAssignedGene'])
t.test(fit$gene$design[fit$gene$design[, 'RegionHIPPO'] == 0, 'mean_totalAssignedGene'], fit_span$gene$design[fit_span$gene$design[, 'RegionHIPPO'] == 0, 'mean_totalAssignedGene'])
t.test(fit$gene$design[, 'mean_totalAssignedGene'] ~ fit$gene$design[, 'RegionHIPPO'])
t.test(fit_span$gene$design[, 'mean_totalAssignedGene'] ~ fit_span$gene$design[, 'RegionHIPPO'])



# plot(-log10(pcheck$global_fdr), -log10(pcheck$adj.P.Val), col = c('gene' = 'blue', 'exon' = 'orange', 'jxn' = 'grey20', 'tx' = 'light blue')[pcheck$type], pch = c('adult' = 21, 'fetal' = 22)[pcheck$age])
# abline(a = 0, b = 1, col = 'red')

table('Global FDR' = pcheck$global_fdr < 0.05, 'FDR' = pcheck$adj.P.Val < 0.05, 'Feature type' = pcheck$type)
table('Global FDR' = pcheck$global_fdr < 0.01, 'FDR' = pcheck$adj.P.Val < 0.01, 'Feature type' = pcheck$type)

table('Global Bonf' = pcheck$global_bonf < 0.05, 'Bonf' = pcheck$P.Bonf < 0.05, 'Feature type' = pcheck$type)
table('Global Bonf' = pcheck$global_bonf < 0.01, 'Bonf' = pcheck$P.Bonf < 0.01,  'Feature type' = pcheck$type)



## Numbers for BOG2018 abstract
pinfo <- lapply(unique(pcheck_both$type), function(feat)  {
    pinfo <- subset(pcheck_both, type == feat)
    pinfo[sign(pinfo$F) == sign(pinfo$span_F) & pinfo$span_P.Value < 0.05 & pinfo$P.Bonf < 0.01, ]
})
names(pinfo) <- unique(pcheck_both$type)


rses <- lapply(unique(pcheck_both$type), load_foo)
names(rses) <- unique(pcheck_both$type)

sapply(pinfo, nrow)

if(!file.exists('rda/de_genes.Rdata')) {
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
    save(pinfo, de_genes, file = 'rda/de_genes.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/de_genes.Rdata'))
    load('rda/de_genes.Rdata', verbose = TRUE)
}

## Pretty venn code
venn_cols <- brewer.pal('Set1', n = 4)
names(venn_cols) <- names(de_genes)
make_venn <- function(genes, title = 'DE features grouped by gene id') {
    v <- venn.diagram(genes, filename = NULL,
        main = title,
        col = 'transparent', fill = venn_cols[names(genes)],
        alpha = 0.5, margin = 0,
        main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
        cat.col = venn_cols[names(genes)])
    grid.newpage()
    grid.draw(v)
}

pdf('pdf/venn_de_features.pdf', useDingbats = FALSE)
make_venn(de_genes)
make_venn(de_genes[c('gene', 'exon', 'jxn')])
dev.off()


if(!file.exists('rda/case_genes.Rdata')) {
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

    case_genes_fdr1 <- lapply(names(pinfo), function(feat) {
        m <- match(rownames(outInfo[[feat]][outInfo[[feat]]$adj.P.Val < 0.1, ]), names(rses[[feat]]))
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
    names(case_genes_fdr1) <- names(pinfo)

    save(case_genes, case_genes_fdr1, file = 'rda/case_genes.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/case_genes.Rdata'))
    load('rda/case_genes.Rdata', verbose = TRUE)
}
sapply(case_genes, length)
sapply(case_genes_fdr1, length)


pdf('pdf/venn_de_features_caseControl.pdf', useDingbats = FALSE)
make_venn(case_genes)
make_venn(case_genes[c('gene', 'exon', 'jxn')])
dev.off()

pdf('pdf/venn_de_features_caseControl_fdr1.pdf', useDingbats = FALSE)
make_venn(case_genes_fdr1)
make_venn(case_genes_fdr1[c('gene', 'exon', 'jxn')])
dev.off()




if(!file.exists('rda/me_genes.Rdata')) {
    load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/matrixEqtl_output_interaction_4features.rda', verbose = TRUE)
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

    me_genes_fdr05 <- lapply(names(pinfo), function(feat) {
        m <- match(me[[feat]]$gene[me[[feat]]$FDR < 0.05], names(rses[[feat]]))
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
    names(me_genes_fdr05) <- names(pinfo)

    save(me_genes, me_genes_fdr05, file = 'rda/me_genes.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/me_genes.Rdata'))
    load('rda/me_genes.Rdata', verbose = TRUE)
}
sapply(me_genes, length)
sapply(me_genes_fdr05, length)
length(unique(unlist(me_genes[c('gene', 'exon', 'jxn')])))

pdf('pdf/venn_eQTL_interaction.pdf', useDingbats = FALSE)
make_venn(me_genes, title = 'eQTLs grouped by gene id')
make_venn(me_genes[c('gene', 'exon', 'jxn')], title = 'eQTLs grouped by gene id')
dev.off()

pdf('pdf/venn_eQTL_interaction_fdr05.pdf', useDingbats = FALSE)
make_venn(me_genes_fdr05, title = 'eQTLs grouped by gene id')
make_venn(me_genes_fdr05[c('gene', 'exon', 'jxn')], title = 'eQTLs grouped by gene id')
dev.off()


## HIPPO only
if(!file.exists('rda/me_genes_HIPPO.Rdata')) {
    load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_hippo_4features.rda', verbose = TRUE)
    me_HIPPO <- split(allEqtl, allEqtl$Type)
    names(me_HIPPO) <- tolower(names(me_HIPPO))
    me_genes_HIPPO <- lapply(names(pinfo), function(feat) {
        m <- match(me_HIPPO[[feat]]$gene[me_HIPPO[[feat]]$FDR < 0.01], names(rses[[feat]]))
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
    names(me_genes_HIPPO) <- names(pinfo)

    me_genes_fdr05_HIPPO <- lapply(names(pinfo), function(feat) {
        m <- match(me_HIPPO[[feat]]$gene[me_HIPPO[[feat]]$FDR < 0.05], names(rses[[feat]]))
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
    names(me_genes_fdr05_HIPPO) <- names(pinfo)
    
    with(allEqtl, addmargins(table('FDR 1%' = FDR < 0.01, 'FDR 5%' = FDR < 0.05)))
    
#        FDR 5%
# FDR 1%     FALSE     TRUE      Sum
#   FALSE  3417626  6245573  9663199
#   TRUE         0 12939932 12939932
#   Sum    3417626 19185505 22603131

    with(subset(allEqtl, Type != 'Tx'), addmargins(table('FDR 1%' = FDR < 0.01, 'FDR 5%' = FDR < 0.05)))
    
#        FDR 5%
# FDR 1%     FALSE     TRUE      Sum
#   FALSE  3044949  5445766  8490715
#   TRUE         0 11237357 11237357
#   Sum    3044949 16683123 19728072

    lapply(me_HIPPO, function(x) {
        with(x, addmargins(table('FDR 1%' = FDR < 0.01, 'FDR 5%' = FDR < 0.05)))
    })
    
    # $Exon
  #          FDR 5%
  #   FDR 1%     FALSE     TRUE      Sum
  #     FALSE  1743337  3220128  4963465
  #     TRUE         0  6206375  6206375
  #     Sum    1743337  9426503 11169840
  #
  #   $Gene
  #          FDR 5%
  #   FDR 1%     TRUE     Sum
  #     FALSE  404706  404706
  #     TRUE  1075733 1075733
  #     Sum   1480439 1480439
  #
  #   $Jxn
  #          FDR 5%
  #   FDR 1%    FALSE    TRUE     Sum
  #     FALSE 1301612 1820932 3122544
  #     TRUE        0 3955249 3955249
  #     Sum   1301612 5776181 7077793
  #
  #   $Tx
  #          FDR 5%
  #   FDR 1%    FALSE    TRUE     Sum
  #     FALSE  372677  799807 1172484
  #     TRUE        0 1702575 1702575
  #     Sum    372677 2502382 2875059

    with(subset(allEqtl, Type != 'Tx'), addmargins(table('FDR 1%' = FDR < 0.01, 'FDR 5%' = FDR < 0.05)))

    save(me_genes_HIPPO, me_genes_fdr05_HIPPO, file = 'rda/me_genes_HIPPO.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/me_genes_HIPPO.Rdata'))
    load('rda/me_genes_HIPPO.Rdata', verbose = TRUE)
}
sapply(me_genes_HIPPO, length)
#  gene  exon   jxn    tx
# 11463 14206 10312 12294
sapply(me_genes_fdr05_HIPPO, length)
#  gene  exon   jxn    tx
# 17044 19090 15030 19810
length(unique(unlist(me_genes_HIPPO[c('gene', 'exon', 'jxn')])))
# [1] 17719

pdf('pdf/venn_eQTL_HIPPO.pdf', useDingbats = FALSE)
make_venn(me_genes_HIPPO, title = 'eQTLs grouped by gene id')
make_venn(me_genes_HIPPO[c('gene', 'exon', 'jxn')], title = 'eQTLs grouped by gene id')
dev.off()

pdf('pdf/venn_eQTL_HIPPO_fdr05.pdf', useDingbats = FALSE)
make_venn(me_genes_fdr05_HIPPO, title = 'eQTLs grouped by gene id')
make_venn(me_genes_fdr05_HIPPO[c('gene', 'exon', 'jxn')], title = 'eQTLs grouped by gene id')
dev.off()



## DLPFC only
if(!file.exists('rda/me_genes_DLPFC.Rdata')) {
    load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/eqtl_tables/mergedEqtl_output_dlpfc_4features.rda', verbose = TRUE)
    me_DLPFC <- split(allEqtl, allEqtl$Type)
    names(me_DLPFC) <- tolower(names(me_DLPFC))
    me_genes_DLPFC <- lapply(names(pinfo), function(feat) {
        m <- match(me_DLPFC[[feat]]$gene[me_DLPFC[[feat]]$FDR < 0.01], names(rses[[feat]]))
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
    names(me_genes_DLPFC) <- names(pinfo)

    me_genes_fdr05_DLPFC <- lapply(names(pinfo), function(feat) {
        m <- match(me_DLPFC[[feat]]$gene[me_DLPFC[[feat]]$FDR < 0.05], names(rses[[feat]]))
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
    names(me_genes_fdr05_DLPFC) <- names(pinfo)
    
    with(allEqtl, addmargins(table('FDR 1%' = FDR < 0.01, 'FDR 5%' = FDR < 0.05)))
    #
    # FDR 1%     FALSE     TRUE      Sum
    #   FALSE  2215761  8610116 10825877
    #   TRUE         0 18037662 18037662
    #   Sum    2215761 26647778 28863539
    
    with(subset(allEqtl, Type != 'Tx'), addmargins(table('FDR 1%' = FDR < 0.01, 'FDR 5%' = FDR < 0.05)))
    
#        FDR 5%
# FDR 1%     FALSE     TRUE      Sum
#   FALSE  1996632  7547106  9543738
#   TRUE         0 15766398 15766398
#   Sum    1996632 23313504 25310136
    
    lapply(me_DLPFC, function(x) {
        with(x, addmargins(table('FDR 1%' = FDR < 0.01, 'FDR 5%' = FDR < 0.05)))
    })
    
    # $Exon
   #         FDR 5%
   #  FDR 1%     FALSE     TRUE      Sum
   #    FALSE   963327  4645078  5608405
   #    TRUE         0  8928309  8928309
   #    Sum     963327 13573387 14536714
   #
   #  $Gene
   #         FDR 5%
   #  FDR 1%     TRUE     Sum
   #    FALSE  418559  418559
   #    TRUE  1577963 1577963
   #    Sum   1996522 1996522
   #
   #  $Jxn
   #         FDR 5%
   #  FDR 1%    FALSE    TRUE     Sum
   #    FALSE 1033305 2483469 3516774
   #    TRUE        0 5260126 5260126
   #    Sum   1033305 7743595 8776900
   #
   #  $Tx
   #         FDR 5%
   #  FDR 1%    FALSE    TRUE     Sum
   #    FALSE  219129 1063010 1282139
   #    TRUE        0 2271264 2271264
   #    Sum    219129 3334274 3553403

    save(me_genes_DLPFC, me_genes_fdr05_DLPFC, file = 'rda/me_genes_DLPFC.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/me_genes_DLPFC.Rdata'))
    load('rda/me_genes_DLPFC.Rdata', verbose = TRUE)
}
sapply(me_genes_DLPFC, length)
#  gene  exon   jxn    tx
# 14261 16047 11732 14379
sapply(me_genes_fdr05_DLPFC, length)
#  gene  exon   jxn    tx
# 18281 20171 15769 21578
length(unique(unlist(me_genes_DLPFC[c('gene', 'exon', 'jxn')])))
# [1] 19482

pdf('pdf/venn_eQTL_DLPFC.pdf', useDingbats = FALSE)
make_venn(me_genes_DLPFC, title = 'eQTLs grouped by gene id')
make_venn(me_genes_DLPFC[c('gene', 'exon', 'jxn')], title = 'eQTLs grouped by gene id')
dev.off()

pdf('pdf/venn_eQTL_DLPFC_fdr05.pdf', useDingbats = FALSE)
make_venn(me_genes_fdr05_DLPFC, title = 'eQTLs grouped by gene id')
make_venn(me_genes_fdr05_DLPFC[c('gene', 'exon', 'jxn')], title = 'eQTLs grouped by gene id')
dev.off()



if(!file.exists('rda/in_risk.Rdata')) {
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

    in_riskFDR <- lapply(me, function(eqtl) {
        if(length(eqtl$snps[eqtl$FDR < 0.01]) == 0) return(NULL)
        m <- match(eqtl$snps[eqtl$FDR < 0.01], rownames(snpMap))
        print(table(!is.na(m)))
        m <- m[!is.na(m)]
        res <- snpMap$pos_hg19[m] %in% riskLoci$hg19POS
        names(res) <- eqtl$gene[eqtl$FDR < 0.01]
        return(res)
    })
    qtls <- lapply(me, function(eqtl) {
        eqtl$snps[eqtl$FDR < 0.01]
    })

    save(in_risk, in_riskFDR, qtls, file = 'rda/in_risk.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/in_risk.Rdata'))
    load('rda/in_risk.Rdata', verbose = TRUE)
}


sapply(in_risk, sum)
sapply(in_risk, mean) * 100


sapply(in_riskFDR, sum)
sapply(in_riskFDR, mean) * 100

pdf('pdf/venn_eQTL_interaction_bySNP.pdf', useDingbats = FALSE)
make_venn(qtls, title = 'eQTLs grouped by SNP id')
make_venn(qtls[c('gene', 'exon', 'jxn')], title = 'eQTLs grouped by SNP id')
dev.off()

qtl_v <- venn(qtls[c('gene', 'exon', 'jxn')], show.plot = FALSE)
sum(sapply(attr(qtl_v, 'intersections'), length))


sapply(in_riskFDR, sum) / sapply(in_risk, sum) / ( sapply(in_riskFDR, length) /  sapply(in_risk, length))


v_me <- venn(me_genes[c('gene', 'exon', 'jxn')], show.plot = FALSE)
risk3 <- lapply(list(all = in_risk$gene, FDR = in_riskFDR$gene), function(rk) {
    names(rk[rk]) %in% c(attr(v_me, 'intersections')[['gene:exon:jxn']], attr(v_me, 'intersections')[['jxn']])
} )
sapply(risk3, sum)


## Go analysis
gene_ens <- lapply(names(rses), function(feat) {
    if(feat %in% c('gene', 'exon', 'jxn')) {
        res <- rowRanges(rses[[feat]])$ensemblID
    } else {
        res <- gsub('\\..*', '', rowRanges(rses[[feat]])$gene_id)
    }
    res <- unique(res[!is.na(res)])
    return(res)
})
names(gene_ens) <- names(rses)

sapply(gene_ens, length)

pdf('pdf/venn_expressed_genes_ensembl_ids.pdf', useDingbats = FALSE)
make_venn(gene_ens)
make_venn(gene_ens[c('gene', 'exon', 'jxn')])
dev.off()

## Use genes expressed in all 3 features (aka, exclude tx)
uni <- unique(unlist(gene_ens[c('gene', 'exon', 'jxn')]))
length(uni)


run_go <- function(genes, ont = c('BP', 'MF', 'CC')) {
    ## Change to ENSEMBL ids
    genes_ens <- sapply(genes, function(x) { gsub('\\..*', '', x) })

    ## Group by jxn, exon, gene
    genes_venn <- venn(genes_ens, show.plot = FALSE)

    ## Run GO analysis
    go_cluster <- lapply(ont, function(bp) {
        message(paste(Sys.time(), 'running GO analysis for', bp))
        tryCatch(compareCluster(attr(genes_venn, 'intersections'), fun = "enrichGO",
            universe = uni, OrgDb = 'org.Hs.eg.db',
            ont = bp, pAdjustMethod = "BH",
            pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
            readable = TRUE, keyType = 'ENSEMBL'),
            error = function(e) { return(NULL) })
    })
    names(go_cluster) <- ont
    
    genes_ncbi <- lapply(lapply(attr(genes_venn, 'intersections'), bitr, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db'), function(x) x$ENTREZID)
    
    uni_ncbi <- bitr(uni, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID
    
    go_cluster$KEGG <- tryCatch(compareCluster(genes_ncbi, fun = 'enrichKEGG',
        universe = uni_ncbi, organism = 'hsa', pAdjustMethod = 'BH',
        pvalueCutoff = 0.1, qvalueCutoff = 0.05, keyType = 'ncbi-geneid'),
        error = function(e) { return(NULL) })
        
    return(go_cluster)
}

run_go_novenn <- function(genes, ont = c('BP', 'MF', 'CC')) {
    ## Change to ENSEMBL ids
    genes_ens <- sapply(genes, function(x) { gsub('\\..*', '', x) })

    ## Run GO analysis
    go_cluster <- lapply(ont, function(bp) {
        message(paste(Sys.time(), 'running GO analysis for', bp))
        tryCatch(compareCluster(genes_ens, fun = "enrichGO",
            universe = uni, OrgDb = 'org.Hs.eg.db',
            ont = bp, pAdjustMethod = "BH",
            pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
            readable = TRUE, keyType = 'ENSEMBL'),
            error = function(e) { return(NULL) })
    })
    names(go_cluster) <- ont
    
    genes_ncbi <- lapply(lapply(genes_ens, bitr, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db'), function(x) x$ENTREZID)
    
    uni_ncbi <- bitr(uni, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID
    
    go_cluster$KEGG <- tryCatch(compareCluster(genes_ncbi, fun = 'enrichKEGG',
        universe = uni_ncbi, organism = 'hsa', pAdjustMethod = 'BH',
        pvalueCutoff = 0.1, qvalueCutoff = 0.05, keyType = 'ncbi-geneid'),
        error = function(e) { return(NULL) })
        
    return(go_cluster)
}



## Development DE genes
if(!file.exists('rda/go_de_genes.Rdata')) {
    system.time( go_de_genes <- run_go(de_genes[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_de_genes.Rdata'))
    save(go_de_genes, file = 'rda/go_de_genes.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_de_genes.Rdata'))
    load('rda/go_de_genes.Rdata', verbose = TRUE)
}
sapply(go_de_genes, class)

if(!file.exists('rda/go_de_genes_novenn.Rdata')) {
    system.time( go_de_genes_novenn <- run_go_novenn(de_genes[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_de_genes_novenn.Rdata'))
    save(go_de_genes_novenn, file = 'rda/go_de_genes_novenn.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_de_genes_novenn.Rdata'))
    load('rda/go_de_genes_novenn.Rdata', verbose = TRUE)
}
sapply(go_de_genes_novenn, class)

## Case-control DE genes
if(!file.exists('rda/go_case_genes.Rdata')) {
    system.time( go_case_genes <- run_go(case_genes[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_case_genes.Rdata'))
    save(go_case_genes, file = 'rda/go_case_genes.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_case_genes.Rdata'))
    load('rda/go_case_genes.Rdata', verbose = TRUE)
}
sapply(go_case_genes, class)

if(!file.exists('rda/go_case_genes_novenn.Rdata')) {
    system.time( go_case_genes_novenn <- run_go_novenn(case_genes[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_case_genes_novenn.Rdata'))
    save(go_case_genes_novenn, file = 'rda/go_case_genes_novenn.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_case_genes_novenn.Rdata'))
    load('rda/go_case_genes_novenn.Rdata', verbose = TRUE)
}
sapply(go_case_genes_novenn, class)

if(!file.exists('rda/go_case_genes_fdr1.Rdata')) {
    system.time( go_case_genes_fdr1 <- run_go(case_genes_fdr1[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_case_genes_fdr1.Rdata'))
    save(go_case_genes_fdr1, file = 'rda/go_case_genes_fdr1.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_case_genes_fdr1.Rdata'))
    load('rda/go_case_genes_fdr1.Rdata', verbose = TRUE)
}
sapply(go_case_genes_fdr1, class)

if(!file.exists('rda/go_case_genes_fdr1_novenn.Rdata')) {
    system.time( go_case_genes_fdr1_novenn <- run_go_novenn(case_genes_fdr1[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_case_genes_fdr1_novenn.Rdata'))
    save(go_case_genes_fdr1_novenn, file = 'rda/go_case_genes_fdr1_novenn.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_case_genes_fdr1_novenn.Rdata'))
    load('rda/go_case_genes_fdr1_novenn.Rdata', verbose = TRUE)
}
sapply(go_case_genes_fdr1_novenn, class)

## eQTL interaction genes
if(!file.exists('rda/go_me_genes.Rdata')) {
    system.time( go_me_genes <- run_go(me_genes[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_me_genes.Rdata'))
    save(go_me_genes, file = 'rda/go_me_genes.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_me_genes.Rdata'))
    load('rda/go_me_genes.Rdata', verbose = TRUE)
}
sapply(go_me_genes, class)

if(!file.exists('rda/go_me_genes_novenn.Rdata')) {
    system.time( go_me_genes_novenn <- run_go_novenn(me_genes[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_me_genes_novenn.Rdata'))
    save(go_me_genes_novenn, file = 'rda/go_me_genes_novenn.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_me_genes_novenn.Rdata'))
    load('rda/go_me_genes_novenn.Rdata', verbose = TRUE)
}
sapply(go_me_genes_novenn, class)


if(!file.exists('rda/go_me_genes_fdr05.Rdata')) {
    system.time( go_me_genes_fdr05 <- run_go(me_genes_fdr05[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_me_genes_fdr05.Rdata'))
    save(go_me_genes_fdr05, file = 'rda/go_me_genes_fdr05.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_me_genes_fdr05.Rdata'))
    load('rda/go_me_genes_fdr05.Rdata', verbose = TRUE)
}
sapply(go_me_genes_fdr05, class)

if(!file.exists('rda/go_me_genes_fdr05_novenn.Rdata')) {
    system.time( go_me_genes_fdr05_novenn <- run_go_novenn(me_genes_fdr05[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_me_genes_fdr05_novenn.Rdata'))
    save(go_me_genes_fdr05_novenn, file = 'rda/go_me_genes_fdr05_novenn.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_me_genes_fdr05_novenn.Rdata'))
    load('rda/go_me_genes_fdr05_novenn.Rdata', verbose = TRUE)
}
sapply(go_me_genes_fdr05_novenn, class)


simplify_go <- function(x) {
    gsub('jxn', 'J', gsub('exon', 'E', gsub('gene', 'G', x)))
}

plot_go <- function(go_cluster, cat = 10) {
    lapply(names(go_cluster), function(bp) {
        go <- go_cluster[[bp]]
        if(is.null(go)) {
            message(paste(Sys.time(), 'found no results for', bp))
            return(NULL)
        }

        ## Simplify names
        go@compareClusterResult$Cluster <- simplify_go(go@compareClusterResult$Cluster)
        names(go@geneClusters) <- simplify_go(names(go@geneClusters))

        print(plot(go, title = paste('ontology:', bp), font.size = 18, showCategory = cat, includeAll = TRUE))
        return(NULL)
    })
}


pdf('pdf/go_de_genes.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_go(go_de_genes)
dev.off()

pdf('pdf/go_case_genes.pdf', width = 14, height = 8, useDingbats = FALSE)
plot_go(go_case_genes)
dev.off()

pdf('pdf/go_case_genes_fdr1.pdf', width = 14, height = 8, useDingbats = FALSE)
plot_go(go_case_genes_fdr1)
dev.off()

pdf('pdf/go_me_genes.pdf', width = 14, height = 8, useDingbats = FALSE)
plot_go(go_me_genes)
dev.off()

pdf('pdf/go_me_genes_fdr05.pdf', width = 14, height = 8, useDingbats = FALSE)
plot_go(go_me_genes_fdr05)
dev.off()


## all
pdf('pdf/go_all_de_genes.pdf', width = 14, height = 50, useDingbats = FALSE)
plot_go(go_de_genes, cat = NULL)
dev.off()

# pdf('pdf/go_all_case_genes.pdf', width = 14, height = 8, useDingbats = FALSE)
# plot_go(go_case_genes, cat = NULL)
# dev.off()

pdf('pdf/go_all_case_genes_fdr1.pdf', width = 18, height = 35, useDingbats = FALSE)
plot_go(go_case_genes_fdr1, cat = NULL)
dev.off()

pdf('pdf/go_all_me_genes.pdf', width = 14, height = 8, useDingbats = FALSE)
plot_go(go_me_genes, cat = NULL)
dev.off()

pdf('pdf/go_all_me_genes_fdr05.pdf', width = 14, height = 8, useDingbats = FALSE)
plot_go(go_me_genes_fdr05, cat = NULL)
dev.off()



## novenn

shorten_go <- function(go, query, short = query) {
    i <- grep(query, go@compareClusterResult$Description)
    go@compareClusterResult$Description[i] <- paste0(go@compareClusterResult$ID[i], ': ', short)
    return(go)
}

go_de_genes_novenn$BP <- shorten_go(go_de_genes_novenn$BP, 'G-protein coupled receptor signaling pathway')

pdf('pdf/go_de_genes_novenn.pdf', width = 11, height = 9, useDingbats = FALSE)
plot_go(go_de_genes_novenn)
dev.off()

pdf('pdf/go_case_genes_novenn.pdf', width = 14, height = 8, useDingbats = FALSE)
plot_go(go_case_genes_novenn)
dev.off()

pdf('pdf/go_case_genes_fdr1_novenn.pdf', width = 14, height = 8, useDingbats = FALSE)
plot_go(go_case_genes_fdr1_novenn)
dev.off()

pdf('pdf/go_me_genes_novenn.pdf', width = 14, height = 8, useDingbats = FALSE)
plot_go(go_me_genes_novenn)
dev.off()

pdf('pdf/go_me_genes_fdr05_novenn.pdf', width = 14, height = 8, useDingbats = FALSE)
plot_go(go_me_genes_fdr05_novenn)
dev.off()

## all
pdf('pdf/go_all_de_genes_novenn.pdf', width = 14, height = 50, useDingbats = FALSE)
plot_go(go_de_genes_novenn, cat = NULL)
dev.off()

pdf('pdf/go_all_case_genes_fdr1_novenn.pdf', width = 18, height = 15, useDingbats = FALSE)
plot_go(go_case_genes_fdr1_novenn, cat = NULL)
dev.off()

pdf('pdf/go_all_me_genes_novenn.pdf', width = 14, height = 11, useDingbats = FALSE)
plot_go(go_me_genes_novenn, cat = NULL)
dev.off()

pdf('pdf/go_all_me_genes_fdr05_novenn.pdf', width = 14, height = 13, useDingbats = FALSE)
plot_go(go_me_genes_fdr05_novenn, cat = NULL)
dev.off()


## Volcano plots
rep_brain <- function(x) {
    ifelse(x, 'BrainSpan rep.', 'BrainSpan not rep.')
}

pdf('pdf/volcano_plots_F.pdf', width = 14, height = 8, useDingbats = FALSE)

ggplot(subset(pcheck_both, type != 'tx'), aes(x = F, y = -log10(P.Value), color = P.Bonf < 0.01)) + geom_point() + facet_grid(rep_brain(span_P.Value < 0.05) ~ type, scale = 'free') + theme_bw(base_size = 18) + scale_colour_manual(values = c('grey20', 'red'), name = 'F Bonf. <1%') + ylab('-log10 p value') + scale_x_log10()

## Subset to smaller F values
ggplot(subset(pcheck_both, type != 'tx' & F < 1100), aes(x = F, y = -log10(P.Value), color = P.Bonf < 0.01)) + geom_point() + facet_grid(rep_brain(span_P.Value < 0.05) ~ type, scale = 'free') + theme_bw(base_size = 18) + scale_colour_manual(values = c('grey20', 'red'), name = 'F Bonf. <1%') + ylab('-log10 p value') + scale_x_log10()

## from limma::topTableF() the rowSums of the abs coefs is the lfc for the F
## although it doesn't look like the F-volcano in Figure 5 (page 16) of
## https://www.bioconductor.org/packages/release/bioc/vignettes/maanova/inst/doc/maanova.pdf
ggplot(subset(pcheck_both, type != 'tx'), aes(x = abs(Age.RegionHIPPO) + abs(RegionHIPPO.fetal) + abs(RegionHIPPO.birth) + abs(RegionHIPPO.infant) + abs(RegionHIPPO.child) + abs(RegionHIPPO.teen) + abs(RegionHIPPO.adult), y = -log10(P.Value), color = P.Bonf < 0.01)) + geom_point() + facet_grid(rep_brain(span_P.Value < 0.05) ~ type, scale = 'free') + theme_bw(base_size = 18) + scale_colour_manual(values = c('grey20', 'red'), name = 'F Bonf. <1%') + ylab('-log10 p value') + xlab('LFC: sum of absolute log2 FC for interaction terms')

dev.off()


## This one separates the BrainSeq and the BrainSpan effects
make_volcano2 <- function(df, v) {
    ggplot(df, aes(x = logfc, y = -log10(pval), color = sigBrainSeq)) + geom_point() + facet_grid(sigBrainSpan ~ type, scales = 'free') + theme_bw(base_size = 18) + scale_colour_manual(values = c('grey20', 'red'), name = 'F Bonf. <1%') + ylab('-log10 p value') + xlab(paste('log2 FC', v))
}

volcano_custom <- function(var, foo = make_volcano) {
    df <- do.call(rbind, lapply(names(fit)[-4], function(f) {
        data.frame(
            pval = fit[[f]]$p.value[, var],
            logfc = fit[[f]]$coefficients[, var],
            type = f,
            stringsAsFactors = FALSE
        )
    }))
    df$name <- rownames(df)
    m <- match(paste0(df$type, '.', df$name), rownames(pcheck_both))
    df$sig <- with(pcheck_both[m, ], P.Bonf < 0.01 & span_P.Value < 0.05)
    df$sigBrainSeq <- with(pcheck_both[m, ], P.Bonf < 0.01)
    df$sigBrainSpan <- rep_brain(with(pcheck_both[m, ], span_P.Value < 0.05))

    g <- foo(df, v = var)
    print(g)
    return(df)
}

pdf('pdf/volcano_plots.pdf', width = 14, height = 8, useDingbats = FALSE)
volcano_info <- lapply(colnames(fit$gene$coef)[grep(':', colnames(fit$gene$coef))], volcano_custom, foo = make_volcano2)
dev.off()


## Some basic exploration (prior to DE)
source('pca_funs.R')

if(!file.exists('rda/pcas.Rdata')) {
    pcas <- lapply(names(exprsNorm), function(type) {
        message(paste(Sys.time(), 'processing', type))
        pc_function(exprsNorm[[type]])
    })
    names(pcas) <- names(exprsNorm)

    pcas_span <- lapply(names(exprsNorm_span), function(type) {
        message(paste(Sys.time(), 'processing', type))
        pc_function(exprsNorm_span[[type]])
    })
    names(pcas_span) <- names(exprsNorm_span)

    save(pcas, pcas_span, file = 'rda/pcas.Rdata')
} else {
    load('rda/pcas.Rdata', verbose = TRUE)
}

if(FALSE) {
    ## Check PCs and sex chrs
    load('rda/pcas.Rdata', verbose = TRUE)
    source('load_funs.R')
    library('SummarizedExperiment')
    rse_gene <- load_foo('gene')
    rse_exon <- load_foo('exon')
    
    tapply(pcas$gene$u[, 7], seqnames(rse_gene), summary)
    
    # $chr1
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0137080 -0.0022928 -0.0003403 -0.0005273  0.0012507  0.0208097
    #
    # $chr2
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -1.209e-02 -1.996e-03 -3.299e-05 -2.229e-04  1.538e-03  2.820e-02
    #
    # $chr3
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0163721 -0.0020139 -0.0002685 -0.0004326  0.0014694  0.0159849
    #
    # $chr4
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -1.136e-02 -1.845e-03  2.400e-04  1.145e-05  1.880e-03  1.597e-02
    #
    # $chr5
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -1.499e-02 -2.043e-03  3.624e-05 -2.097e-04  1.665e-03  1.612e-02
    #
    # $chr6
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0124551 -0.0019655 -0.0001863 -0.0002823  0.0014483  0.0183057
    #
    # $chr7
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0162665 -0.0023595 -0.0003755 -0.0004568  0.0012753  0.0282686
    #
    # $chr8
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -1.079e-02 -1.751e-03  1.908e-04 -3.222e-05  1.730e-03  1.241e-02
    #
    # $chr9
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0134702 -0.0024986 -0.0005146 -0.0007157  0.0011238  0.0128030
    #
    # $chr10
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -1.418e-02 -2.132e-03  2.646e-06 -3.565e-04  1.470e-03  1.682e-02
    #
    # $chr11
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0155661 -0.0026353 -0.0006692 -0.0008716  0.0009299  0.0137021
    #
    # $chr12
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0121895 -0.0022321 -0.0001950 -0.0004622  0.0013628  0.0109990
    #
    # $chr13
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0148956 -0.0015537  0.0004171  0.0003582  0.0021570  0.0237661
    #
    # $chr14
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -1.050e-02 -1.817e-03  1.712e-04 -4.305e-05  1.801e-03  1.291e-02
    #
    # $chr15
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0180643 -0.0024705 -0.0002243 -0.0004981  0.0017077  0.0108856
    #
    # $chr16
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0102932 -0.0030572 -0.0013213 -0.0013450  0.0004023  0.0167246
    #
    # $chr17
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0176315 -0.0028449 -0.0010630 -0.0011587  0.0006508  0.0185407
    #
    # $chr18
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0173907 -0.0020566 -0.0003483 -0.0003282  0.0014831  0.0182248
    #
    # $chr19
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0149466 -0.0027946 -0.0010850 -0.0012304  0.0002332  0.0213610
    #
    # $chr20
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0325424 -0.0028331 -0.0006072 -0.0011224  0.0009100  0.0115885
    #
    # $chr21
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0151642 -0.0025810 -0.0005833 -0.0007184  0.0015340  0.0166101
    #
    # $chr22
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0197491 -0.0030777 -0.0011081 -0.0012165  0.0005894  0.0161412
    #
    # $chrX
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0324965 -0.0012971  0.0003904  0.0006222  0.0021764  0.2072613
    #
    # $chrY
    #      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
    # -0.210573 -0.176114 -0.113417 -0.111009 -0.039595  0.004763
    #
    # $chrM
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0107657 -0.0034202 -0.0013154 -0.0013057  0.0006736  0.0069107
    
    tapply(pcas$exon$u[, 7], seqnames(rse_exon), summary)
    
    # $chr1
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -5.889e-03 -2.716e-04  9.943e-05  1.685e-04  5.710e-04  4.681e-03
    #
    # $chr2
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0057058 -0.0002929  0.0001193  0.0001816  0.0005987  0.0041749
    #
    # $chr3
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0041222 -0.0002659  0.0001338  0.0001971  0.0006194  0.0047684
    #
    # $chr4
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -3.307e-03 -3.367e-04  5.573e-05  1.272e-04  5.242e-04  4.049e-03
    #
    # $chr5
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -4.280e-03 -3.333e-04  4.305e-05  1.196e-04  5.035e-04  3.989e-03
    #
    # $chr6
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0037068 -0.0002799  0.0001164  0.0001589  0.0005683  0.0037027
    #
    # $chr7
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -7.324e-03 -3.209e-04  8.862e-05  1.609e-04  5.767e-04  4.301e-03
    #
    # $chr8
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -3.278e-03 -3.622e-04 -4.507e-06  6.340e-05  4.467e-04  3.906e-03
    #
    # $chr9
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0036043 -0.0002298  0.0001718  0.0002335  0.0006404  0.0049567
    #
    # $chr10
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -3.911e-03 -3.293e-04  6.528e-05  1.405e-04  5.506e-04  3.577e-03
    #
    # $chr11
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0040893 -0.0002630  0.0001036  0.0001687  0.0005671  0.0044779
    #
    # $chr12
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -2.880e-03 -2.882e-04  9.691e-05  1.601e-04  5.591e-04  3.542e-03
    #
    # $chr13
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -5.766e-03 -3.183e-04  4.286e-05  7.250e-05  4.493e-04  3.347e-03
    #
    # $chr14
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -3.739e-03 -2.990e-04  7.764e-05  1.275e-04  5.177e-04  3.353e-03
    #
    # $chr15
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0027215 -0.0002474  0.0001430  0.0002163  0.0006188  0.0056051
    #
    # $chr16
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0038058 -0.0002275  0.0001547  0.0002158  0.0006137  0.0061025
    #
    # $chr17
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0049125 -0.0002455  0.0001423  0.0002115  0.0006143  0.0048135
    #
    # $chr18
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0025112 -0.0003131  0.0000906  0.0001373  0.0005266  0.0037207
    #
    # $chr19
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0062974 -0.0002367  0.0001496  0.0001981  0.0006164  0.0043090
    #
    # $chr20
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -3.737e-03 -2.774e-04  9.633e-05  1.803e-04  5.643e-04  8.490e-03
    #
    # $chr21
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -4.745e-03 -3.348e-04  3.195e-05  1.076e-04  4.958e-04  4.030e-03
    #
    # $chr22
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -0.0036130 -0.0002998  0.0001396  0.0001733  0.0006049  0.0038809
    #
    # $chrX
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -6.349e-02 -4.177e-04 -4.118e-05 -1.800e-04  4.215e-04  6.926e-03
    #
    # $chrY
    #      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
    # -0.001691  0.024516  0.040186  0.033494  0.046594  0.057698
    #
    # $chrM
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -1.558e-03 -3.314e-04  2.297e-04  9.497e-05  5.931e-04  1.813e-03
    
    pdf('pdf/pca_check_sex.pdf', useDingbats = FALSE)
    boxplot(pcas$gene$u[, 7] ~ as.factor(seqnames(rse_gene)), las = 2, main = 'BrainSeq Phase 2, PC7, gene')
    boxplot(pcas$exon$u[, 7] ~ as.factor(seqnames(rse_exon)), las = 2, main = 'BrainSeq Phase 2, PC7, exon')
    dev.off()
}



main_lab <- function(type) {
    res <- paste0('PC (', type, '): ')
    if(type != 'tx') {
        res <- paste0(res, 'log2(CPM + 0.5)')
    } else {
        res <- paste0(res, 'log2(TPM + 0.5)')
    }
    return(res)
}

rses_span <- lapply(unique(pcheck_both$type), load_span)
names(rses_span) <- unique(pcheck_both$type)


pal <- brewer.pal('Paired', n = 8)
plot_pca <- function(type, pc, loc = 'bottomleft', pan = 1, rs = rses, nc = 2) {
    leg <- with(colData(rs[[type]]), paste0(Region, '-', ifelse(Age < 0, 'pre', 'post')))
    leg2 <- with(colData(rs[[type]]), paste0(Region, '-', Sex))
    race <- as.character(colData(rs[[type]])$Race)
    race[race == 'African'] <- 'AA'
    race[!race %in% c('CAUC', 'AA')] <- 'Other'
    #leg4 <- with(colData(rs$gene), paste0(Region, '-', Dx))

    palette(pal[1:4])
    pc_plot(pca = pc[[type]], legend = leg, color = leg,
            main = main_lab(type),
            position = loc, type='variable', ptsize = 0.75,
            ncol = nc, legend.pan = pan, lwd = 3, cex = 0.75)
    palette(pal[5:8])
    pc_plot(pca = pc[[type]], legend = leg2, color = leg2,
            main = main_lab(type),
            position = loc, type='variable', ptsize = 0.75,
            ncol = nc, legend.pan = pan, lwd = 3, cex = 0.75)
    palette(brewer.pal('Set1', n = 5)[3:5])
    pc_plot(pca = pc[[type]], legend = race,
            color = race,
            main = main_lab(type),
            position = loc, type='variable', ptsize = 0.75,
            ncol = nc, legend.pan = pan, lwd = 3, cex = 0.75)
    var_plot(pca = pc[[type]])
    return(NULL)
}

## plot by dataset
pdf('pdf/pcas.pdf')
plot_pca('gene', pcas)
plot_pca('exon', pcas, 'bottomright')
plot_pca('jxn', pcas, 'topright')
plot_pca('tx', pcas, 'topleft', pan = 4, nc = 1)
dev.off()


pdf('pdf/pcas_span.pdf')
plot_pca('gene', pcas_span, 'topright', rs = rses_span)
plot_pca('exon', pcas_span, 'topright', rs = rses_span)
plot_pca('jxn', pcas_span, pan = 4, rs = rses_span)
plot_pca('tx', pcas_span, 'bottomright', pan = 4, nc = 1, rs = rses_span)
dev.off()


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


## Re-loading if necessary
if(FALSE) {
    f <- dir('rda', full.names = TRUE)
    f <- f[!grepl('limma', f)]
    for(ff in f) load(ff, verbose = TRUE)
    rm(ff, f)
}
