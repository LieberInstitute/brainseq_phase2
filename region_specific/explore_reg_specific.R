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
    raw <- mapply(function(age, type) {
        load(paste0('rda/limma_region_specific_', age, '_', type, '.Rdata'))
        top$age <- age
        top$type <- type
        return(list(top = top, fit = fit, exprsNorm = exprsNorm))
    }, rep(c('adult', 'fetal'), each = 4), rep(c('gene', 'exon', 'jxn', 'tx'), 2), SIMPLIFY = FALSE)
    names(raw) <- paste0(rep(c('adult', 'fetal'), each = 4), "_", rep(c('gene', 'exon', 'jxn', 'tx'), 2))
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
    raw_span <- mapply(function(age, type) {
        load(paste0('rda/span_limma_region_specific_', age, '_', type, '.Rdata'))
        top$age <- age
        top$type <- type
        return(list(top = top, fit = fit, exprsNorm = exprsNorm))
        }, rep(c('adult', 'fetal'), each = 4), rep(c('gene', 'exon', 'jxn', 'tx'), 2), SIMPLIFY = FALSE)
    names(raw_span) <- paste0(rep(c('adult', 'fetal'), each = 4), "_", rep(c('gene', 'exon', 'jxn', 'tx'), 2))
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


summary(pcheck_both$logFC)
summary(pcheck_both$span_logFC)

# pdf('pdf/compare_with_span_logFC.pdf', useDingbats = FALSE)
png('pdf/compare_with_span_logFC.png', type = 'cairo')
ggplot(pcheck_both, aes(x = logFC, y = span_logFC, alpha = 1/20)) +
    facet_grid(age ~ type, scales = 'free') + ylab('BrainSpan log FC') +
    xlab('BrainSeq log FC') + geom_point() +
    geom_smooth(method=lm, se=FALSE)
dev.off()

# pdf('pdf/compare_with_span_logFC_noTx.pdf', useDingbats = FALSE)
png('pdf/compare_with_span_logFC_noTx.png', type = 'cairo')
ggplot(subset(pcheck_both, type != 'tx'), aes(x = logFC, y = span_logFC,
    alpha = 1/20)) +
    facet_grid(age ~ type, scales = 'free') + ylab('BrainSpan log FC') +
    xlab('BrainSeq log FC') + geom_point() +
    geom_smooth(method=lm, se=FALSE)
dev.off()


pdf('pdf/compare_with_span_logFC_density.pdf', useDingbats = FALSE)
ggplot(pcheck_both, aes(x = logFC, y = span_logFC)) +
    facet_grid(age ~ type, scales = 'free') + ylab('BrainSpan log FC') +
    xlab('BrainSeq log FC') + stat_density2d() +
    geom_smooth(method=lm, se=FALSE)
dev.off()

pdf('pdf/compare_with_span_logFC_density_noTx.pdf', useDingbats = FALSE)
ggplot(subset(pcheck_both, type != 'tx'), aes(x = logFC, y = span_logFC)) +
    facet_grid(age ~ type, scales = 'free') + ylab('BrainSpan log FC') +
    xlab('BrainSeq log FC') + stat_density2d() +
    geom_smooth(method=lm, se=FALSE)
dev.off()


## Replication (p < 0.05) & same direction
if(!file.exists('rda/rep_span.Rdata')) {
    rep_span <- do.call(rbind, mapply(function(pvar, cut, type_sub, age_sub) {
        pinfo <- subset(pcheck_both, type == type_sub & age == age_sub)
        n <- sum(sign(pinfo$logFC) == sign(pinfo$span_logFC) & pinfo$span_P.Value < 0.05 & pinfo[, pvar] < cut)
        n_sign <- sum(sign(pinfo$logFC) == sign(pinfo$span_logFC) & pinfo[, pvar] < cut)
        data.frame(pvar = pvar, cutoff = cut, type = type_sub, age = age_sub, replicated = n, number_de = sum(pinfo[, pvar] < cut), total = nrow(pinfo), replicated_sign = n_sign, stringsAsFactors = FALSE)
    }, pvar = rep(rep(c('adj.P.Val', 'P.Bonf'), each = 6), 4 * 2), cut = rep(c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001), 2 * 4 * 2), type_sub = rep(rep(unique(pcheck_both$type), each = 6 * 2), 2), age_sub = rep(unique(pcheck_both$age), each = 4 * 6 * 2), SIMPLIFY = FALSE, USE.NAMES = FALSE))
    
    ## Fix order of features and ages
    rep_span$age <- factor(ifelse(rep_span$age == 'fetal', 'prenatal', 'adult'), levels = c('prenatal', 'adult'))
    rep_span$type <- factor(rep_span$type, levels = c('gene', 'exon', 'jxn', 'tx'))
    save(rep_span, file = 'rda/rep_span.Rdata')
} else {
    message(paste(Sys.time(), 'rda/rep_span.Rdata'))
    load('rda/rep_span.Rdata', verbose = TRUE)
}




pdf('pdf/replication_exploration.pdf', width = 14, useDingbats = FALSE)
ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated / number_de, color = pvar)) + facet_grid(age ~ type) + ylab('Replication rate') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(color='P-value method') + ylim(c(0, 1))

ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated_sign / number_de, color = pvar)) + facet_grid(age ~ type) + ylab('Replication rate (sign only)') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(color='P-value method') + ylim(c(0, 1))

ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = number_de, color = pvar)) + facet_grid(age ~ type) + ylab('Number of DE features') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10() + labs(color='P-value method')

ggplot(rep_span, aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = number_de / total * 100, color = pvar)) + facet_grid(age ~ type) + ylab('Percent of DE features') + xlab('p-threshold') + geom_point() + theme_grey(base_size = 18) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(color='P-value method') + ylim(c(0, 100))

dev.off()


pdf('pdf/replication_exploration_subset.pdf', width = 8, height =  10, useDingbats = FALSE)
ggplot(subset(rep_span, pvar == 'P.Bonf' & type != 'tx'), aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated / number_de)) + facet_grid(type ~ age) + ylab('Replication rate') + xlab('p-value threshold') + geom_point() + theme_bw(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(c(0, 1)) + geom_line(aes(y = replicated / number_de, x = rep(1:6, 6)))

ggplot(subset(rep_span, pvar == 'P.Bonf' & type != 'tx'), aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated_sign / number_de)) + facet_grid(type ~ age) + ylab('Replication rate (sign only)') + xlab('p-value threshold') + geom_point() + theme_bw(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(c(0, 1)) + geom_line(aes(y = replicated_sign / number_de, x = rep(1:6, 6)))

ggplot(subset(rep_span, pvar == 'P.Bonf' & type != 'tx'), aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = number_de)) + facet_grid(type ~ age) + ylab('Number of DE features') + xlab('p-value threshold') + geom_point() + theme_bw(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  geom_line(aes(y = number_de, x = rep(1:6, 6)))

ggplot(subset(rep_span, pvar == 'P.Bonf' & type != 'tx'), aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = number_de / total * 100)) + facet_grid(type ~ age) + ylab('Percent of DE features') + xlab('p-value threshold') + geom_point() + theme_bw(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_line(aes(y = number_de / total * 100, x = rep(1:6, 6))) + ylim(c(0, 100))

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

# plot(fit$fetal_gene$design[, 'snpPC1'], fit$fetal_gene$design[, 'snpPC2'])
# plot(fit_span$fetal_gene$design[, 'snpPC1'], fit_span$fetal_gene$design[, 'snpPC2'])

t.test(fit$fetal_gene$design[fit$fetal_gene$design[, 'RegionHIPPO'] == 1, 'mean_totalAssignedGene'], fit_span$fetal_gene$design[fit_span$fetal_gene$design[, 'RegionHIPPO'] == 1, 'mean_totalAssignedGene'])
t.test(fit$fetal_gene$design[fit$fetal_gene$design[, 'RegionHIPPO'] == 0, 'mean_totalAssignedGene'], fit_span$fetal_gene$design[fit_span$fetal_gene$design[, 'RegionHIPPO'] == 0, 'mean_totalAssignedGene'])
t.test(fit$fetal_gene$design[, 'mean_totalAssignedGene'] ~ fit$fetal_gene$design[, 'RegionHIPPO'])
t.test(fit_span$fetal_gene$design[, 'mean_totalAssignedGene'] ~ fit_span$fetal_gene$design[, 'RegionHIPPO'])



# plot(-log10(pcheck$global_fdr), -log10(pcheck$adj.P.Val), col = c('gene' = 'blue', 'exon' = 'orange', 'jxn' = 'grey20', 'tx' = 'light blue')[pcheck$type], pch = c('adult' = 21, 'fetal' = 22)[pcheck$age])
# abline(a = 0, b = 1, col = 'red')

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
# adult 1612 15442 5561 1739
# fetal   32    71   18    3

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

if(!file.exists('rda/de_genes.Rdata')) {
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
    save(pinfo, de_genes, file = 'rda/de_genes.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/de_genes.Rdata'))
    load('rda/de_genes.Rdata', verbose = TRUE)
}

sapply(de_genes, function(x) sapply(x, length))
#      adult fetal
# gene  1612    32
# exon  2686    11
# jxn   1897     6
# tx    1414     1


if(!file.exists('rda/de_genes_sign.Rdata')) {
    de_genes_sign <- lapply(names(pinfo[[1]]), function(agegrp) {
        res2 <- lapply(names(pinfo), function(feat) {
            res3 <- lapply(c('DLPFC', 'HIPPO'), function(reg) {
                curr <- pinfo[[feat]][[agegrp]]
                de <- rownames(curr[ sign(curr$logFC) == ifelse(reg == 'DLPFC', -1, 1), ])
                m <- match(gsub('.*gene.|.*exon.|.*jxn.|.*tx.', '', de),
                    names(rses[[feat]]))
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
            names(res3) <- c('DLPFC', 'HIPPO')
            return(res3)
        })
        names(res2) <- names(pinfo)
        return(do.call(c, res2))
    })
    names(de_genes_sign) <- names(pinfo$gene)
    save(de_genes_sign, file = 'rda/de_genes_sign.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/de_genes_sign.Rdata'))
    load('rda/de_genes_sign.Rdata', verbose = TRUE)
}

sapply(de_genes_sign, function(x) sapply(x, length))
#            adult fetal
# gene.DLPFC   823    30
# gene.HIPPO   789     2
# exon.DLPFC  1335     8
# exon.HIPPO  1368     3
# jxn.DLPFC    949     4
# jxn.HIPPO    959     2
# tx.DLPFC     827     0
# tx.HIPPO     596     1



## Pretty venn code
venn_cols <- brewer.pal('Set1', n = 4)
names(venn_cols) <- names(de_genes[[1]])
make_venn <- function(genes, title = 'DE features grouped by gene id') {
    v <- venn.diagram(genes, filename = NULL,
        main = title,
        col = 'transparent', fill = venn_cols[seq_len(length(genes))],
        alpha = 0.5, margin = 0,
        main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
        cat.col = venn_cols[seq_len(length(genes))])
    grid.newpage()
    grid.draw(v)
}

pdf('pdf/venn_de_features.pdf', useDingbats = FALSE)
make_venn(de_genes$adult, title = 'DE features grouped by gene id (adult)')
make_venn(de_genes$adult[c('gene', 'exon', 'jxn')], title = 'DE features grouped by gene id (adult)')
make_venn(de_genes$fetal, title = 'DE features grouped by gene id (prenatal)')
make_venn(de_genes$fetal[c('gene', 'exon', 'jxn')], title = 'DE features grouped by gene id (prenatal)')


make_venn(de_genes_sign$adult[c(1, 3, 5, 7)], title = 'DE features grouped by gene id (adult)')
make_venn(de_genes_sign$adult[c(1, 3, 5)], title = 'DE features grouped by gene id (adult)')
make_venn(de_genes_sign$fetal[c(1, 3, 5, 7)], title = 'DE features grouped by gene id (prenatal)')
make_venn(de_genes_sign$fetal[c(1, 3, 5)], title = 'DE features grouped by gene id (prenatal)')

make_venn(de_genes_sign$adult[c(2, 4, 6, 8)], title = 'DE features grouped by gene id (adult)')
make_venn(de_genes_sign$adult[c(2, 4, 6)], title = 'DE features grouped by gene id (adult)')
make_venn(de_genes_sign$fetal[c(2, 4, 6, 8)], title = 'DE features grouped by gene id (prenatal)')
make_venn(de_genes_sign$fetal[c(2, 4, 6)], title = 'DE features grouped by gene id (prenatal)')


make_venn(de_genes_sign$adult[c(1:2, 3:4)], title = 'DE features grouped by gene id (adult)')
make_venn(de_genes_sign$fetal[c(1:2, 3:4)], title = 'DE features grouped by gene id (prenatal)')

make_venn(de_genes_sign$adult[c(1:2, 5:6)], title = 'DE features grouped by gene id (adult)')
make_venn(de_genes_sign$fetal[c(1:2, 5:6)], title = 'DE features grouped by gene id (prenatal)')

make_venn(de_genes_sign$adult[c(3:4, 5:6)], title = 'DE features grouped by gene id (adult)')
make_venn(de_genes_sign$fetal[c(3:4, 5:6)], title = 'DE features grouped by gene id (prenatal)')

make_venn(de_genes_sign$adult[c(1:2, 7:8)], title = 'DE features grouped by gene id (adult)')
make_venn(de_genes_sign$fetal[c(1:2, 7:8)], title = 'DE features grouped by gene id (prenatal)')
dev.off()
system('rm VennDiagram*.log')

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
#  gene  exon   jxn    tx
# 24652 23943 18329 31688

pdf('pdf/venn_expressed_genes_ensembl_ids.pdf', useDingbats = FALSE)
make_venn(gene_ens)
make_venn(gene_ens[c('gene', 'exon', 'jxn')])
dev.off()

## Use genes expressed in all 3 features (aka, exclude tx)
uni <- unique(unlist(gene_ens[c('gene', 'exon', 'jxn')]))
length(uni)
# [1] 27514

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
    
    message(paste(Sys.time(), 'running GO analysis for KEGG'))
    genes_ncbi <- lapply(lapply(attr(genes_venn, 'intersections'), bitr, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db'), function(x) x$ENTREZID)
    
    uni_ncbi <- bitr(uni, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID
    
    go_cluster$KEGG <- tryCatch(compareCluster(genes_ncbi, fun = 'enrichKEGG',
        universe = uni_ncbi, organism = 'hsa', pAdjustMethod = 'BH',
        pvalueCutoff = 0.1, qvalueCutoff = 0.05, keyType = 'ncbi-geneid'),
        error = function(e) { return(NULL) })

    return(go_cluster)
}

run_go_novenn <- function(genes, ont = c('BP', 'MF', 'CC'), set = 'adult') {
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
    
    message(paste(Sys.time(), 'running GO analysis for KEGG'))
    genes_ncbi <- lapply(lapply(genes_ens, bitr, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db'), function(x) x$ENTREZID)
    
    uni_ncbi <- bitr(uni, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID
    
    go_cluster$KEGG <- tryCatch(compareCluster(genes_ncbi, fun = 'enrichKEGG',
        universe = uni_ncbi, organism = 'hsa', pAdjustMethod = 'BH',
        pvalueCutoff = 0.1, qvalueCutoff = 0.05, keyType = 'ncbi-geneid'),
        error = function(e) { return(NULL) })

    return(go_cluster)
}

## Reg specific genes
if(!file.exists('rda/go_de_genes.Rdata')) {
    system.time( go_de_genes_adult <- run_go(de_genes$adult[c('gene', 'exon', 'jxn')]) )
    system.time( go_de_genes_fetal <- run_go(de_genes$fetal[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_de_genes.Rdata'))
    save(go_de_genes_adult, go_de_genes_fetal, file = 'rda/go_de_genes.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_de_genes.Rdata'))
    load('rda/go_de_genes.Rdata', verbose = TRUE)
}
sapply(go_de_genes_adult, class)
sapply(go_de_genes_fetal, class)


## Reg specific genes - by brain region
if(!file.exists('rda/go_de_genes_brain.Rdata')) {
    system.time( go_de_genes_brain_adult <- run_go_novenn(de_genes_sign$adult[1:6]) )
    system.time( go_de_genes_brain_fetal <- run_go_novenn(de_genes_sign$fetal[1:6], set = 'fetal') )
    message(paste(Sys.time(), 'saving rda/go_de_genes.Rdata'))
    save(go_de_genes_brain_adult, go_de_genes_brain_fetal, file = 'rda/go_de_genes_brain.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_de_genes_brain.Rdata'))
    load('rda/go_de_genes_brain.Rdata', verbose = TRUE)
}
sapply(go_de_genes_brain_adult, class)
sapply(go_de_genes_brain_fetal, class)

simplify_go <- function(x) {
    gsub('LPFC|IPPO', '', gsub('jxn', 'J', gsub('exon', 'E', gsub('gene', 'G', x))))
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


pdf('pdf/go_de_genes_adult.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_go(go_de_genes_adult)
dev.off()


shorten_go <- function(go, query, short = query) {
    i <- grep(query, go@compareClusterResult$Description)
    go@compareClusterResult$Description[i] <- paste0(go@compareClusterResult$ID[i], ': ', short)
    return(go)
}

go_de_genes_fetal$MF <- shorten_go(go_de_genes_fetal$MF, 'oxidoreductase activity')

pdf('pdf/go_de_genes_fetal.pdf', width = 12, height = 9, useDingbats = FALSE)
plot_go(go_de_genes_fetal)
dev.off()

go_de_genes_brain_adult$MF <- shorten_go(go_de_genes_brain_adult$MF, 'transcriptional activator activity')

pdf('pdf/go_de_genes_brain_adult.pdf', width = 14, height = 11, useDingbats = FALSE)
plot_go(go_de_genes_brain_adult)
dev.off()

go_de_genes_brain_fetal$MF <- shorten_go(go_de_genes_brain_fetal$MF, 'oxidoreductase activity')
go_de_genes_brain_fetal$MF <- shorten_go(go_de_genes_brain_fetal$MF, 'transcriptional activator activity')
go_de_genes_brain_fetal$MF <- shorten_go(go_de_genes_brain_fetal$MF, 'transcription factor activity')

pdf('pdf/go_de_genes_brain_fetal.pdf', width = 14, height = 11, useDingbats = FALSE)
plot_go(go_de_genes_brain_fetal)
dev.off()

## all
pdf('pdf/go_all_de_genes_adult.pdf', width = 12, height = 50, useDingbats = FALSE)
plot_go(go_de_genes_adult, cat = NULL)
dev.off()

pdf('pdf/go_all_de_genes_fetal.pdf', width = 12, height = 50, useDingbats = FALSE)
plot_go(go_de_genes_fetal, cat = NULL)
dev.off()

pdf('pdf/go_all_de_genes_brain_adult.pdf', width = 14, height = 100, useDingbats = FALSE)
plot_go(go_de_genes_brain_adult, cat = NULL)
dev.off()

pdf('pdf/go_all_de_genes_brain_fetal.pdf', width = 14, height = 50, useDingbats = FALSE)
plot_go(go_de_genes_brain_fetal, cat = NULL)
dev.off()


## Volcano plots
rep_brain <- function(x) {
    ifelse(x, 'BrainSpan rep.', 'BrainSpan not rep.')
}

## This one separates the BrainSeq and the BrainSpan effects
make_volcano2 <- function(df, v) {
    ggplot(df, aes(x = logfc, y = -log10(pval), color = sigBrainSeq)) + geom_point() + facet_grid(sigBrainSpan ~ type + age, scales = 'free') + theme_bw(base_size = 18) + scale_colour_manual(values = c('grey20', 'red'), name = 'Bonf. <1%') + ylab('-log10 p value') + xlab(paste('log2 FC', v))
}

volcano_custom <- function(var, foo = make_volcano) {
    df <- do.call(rbind, lapply(names(fit)[-grep('_tx', names(fit))], function(f) {
        data.frame(
            pval = fit[[f]]$p.value[, var],
            logfc = fit[[f]]$coefficients[, var],
            type = ss(f, '_', 2),
            age = ss(f, '_'),
            stringsAsFactors = FALSE
        )
    }))
    df$name <- rownames(df)
    df$name[df$age == 'fetal'] <- gsub('.$', '', df$name[df$age == 'fetal'])
    m <- match(paste0(df$age, '_', df$type, '.', df$name), rownames(pcheck_both))

    df$sig <- with(pcheck_both[m, ], P.Bonf < 0.01 & span_P.Value < 0.05)
    df$sigBrainSeq <- with(pcheck_both[m, ], P.Bonf < 0.01)
    df$sigBrainSpan <- rep_brain(with(pcheck_both[m, ], span_P.Value < 0.05 & sign(t) == sign(span_t)))

    df <- df[complete.cases(df), ]
    df$age[df$age == 'fetal'] <- 'prenatal'
    
    ## Fix order
    df$age <- factor(df$age, levels = c('prenatal', 'adult'))
    df$type <- factor(df$type, levels = c('gene', 'exon', 'jxn'))

    g <- foo(df, v = var)
    print(g)
    return(df)
}

pdf('pdf/volcano_plots.pdf', width = 16, height = 8, useDingbats = FALSE)
volcano_info <- volcano_custom('RegionHIPPO', foo = make_volcano2)
dev.off()


## Some basic exploration (prior to DE)
source('../development/pca_funs.R')

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


main_lab <- function(type) {
    res <- paste0('PC (', gsub('fetal', 'prenatal', type), '): ')
    if(!grepl('tx', type)) {
        res <- paste0(res, 'log2(CPM + 0.5)')
    } else {
        res <- paste0(res, 'log2(TPM + 0.5)')
    }
    return(res)
}

rses_span <- lapply(unique(pcheck_both$type), function(type) {
    load_span(type, 'fetal')
})
names(rses_span) <- unique(pcheck_both$type)


pal <- brewer.pal('Paired', n = 8)
plot_pca <- function(type, pc, loc = 'bottomleft', pan = 1, rs = rses, nc = 2) {
    leg <- with(colData(rs[[type]]), Region)
    leg2 <- with(colData(rs[[type]]), paste0(Region, '-', Sex))
    race <- as.character(colData(rs[[type]])$Race)
    race[race == 'African'] <- 'AA'
    race[!race %in% c('CAUC', 'AA')] <- 'Other'
    #leg4 <- with(colData(rs$gene), paste0(Region, '-', Dx))

    palette(c('dark orange', 'skyblue3', 'black'))
    pc_plot(pca = pc[[paste0('adult_', type)]], legend = leg, color = leg,
            main = main_lab(paste(type, 'adult')),
            position = loc, type='variable', ptsize = 0.75,
            ncol = nc, legend.pan = pan, lwd = 3, cex = 0.75)
    palette(pal[5:8])
    pc_plot(pca = pc[[paste0('adult_', type)]], legend = leg2, color = leg2,
            main = main_lab(paste(type, 'adult')),
            position = loc, type='variable', ptsize = 0.75,
            ncol = nc, legend.pan = pan, lwd = 3, cex = 0.75)
    palette(brewer.pal('Set1', n = 5)[3:5])
    pc_plot(pca = pc[[paste0('adult_', type)]], legend = race,
            color = race,
            main = main_lab(paste(type, 'adult')),
            position = loc, type='variable', ptsize = 0.75,
            ncol = nc, legend.pan = pan, lwd = 3, cex = 0.75)
    var_plot(pca = pc[[paste0('adult_', type)]])

    palette(c('dark orange', 'skyblue3', 'black'))
    pc_plot(pca = pc[[paste0('fetal_', type)]], legend = leg, color = leg,
            main = main_lab(paste(type, 'prenatal')),
            position = loc, type='variable', ptsize = 0.75,
            ncol = nc, legend.pan = pan, lwd = 3, cex = 0.75)
    palette(pal[5:8])
    pc_plot(pca = pc[[paste0('fetal_', type)]], legend = leg2, color = leg2,
            main = main_lab(paste(type, 'prenatal')),
            position = loc, type='variable', ptsize = 0.75,
            ncol = nc, legend.pan = pan, lwd = 3, cex = 0.75)
    palette(brewer.pal('Set1', n = 5)[3:5])
    pc_plot(pca = pc[[paste0('fetal_', type)]], legend = race,
            color = race,
            main = main_lab(paste(type, 'prenatal')),
            position = loc, type='variable', ptsize = 0.75,
            ncol = nc, legend.pan = pan, lwd = 3, cex = 0.75)
    var_plot(pca = pc[[paste0('fetal_', type)]])
    return(NULL)
}

## plot by dataset
pdf('pdf/pcas.pdf')
plot_pca('gene', pcas)
plot_pca('exon', pcas, 'bottomleft', pan = 4)
plot_pca('jxn', pcas, 'bottomright', pan = 4)
plot_pca('tx', pcas, 'bottomright', pan = 4, nc = 1)
dev.off()


pdf('pdf/pcas_span.pdf')
plot_pca('gene', pcas_span, 'bottomright', pan = 2, rs = rses_span)
plot_pca('exon', pcas_span, 'bottomright', pan = 2, rs = rses_span)
plot_pca('jxn', pcas_span, 'right', pan = 1, rs = rses_span)
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

