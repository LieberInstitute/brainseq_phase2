# qrsh -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
# mkdir -p logs
# Rscript explore_limma_dev.R > logs/explore_limma_dev.txt 2>&1

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
    save(pcheck_both, file = 'rda/pcheck_both.Rdata')
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
table('brainseq' = abs(pcheck_both$Age.RegionHIPPO) > 10, 'span' = abs(pcheck_both$span_Age.RegionHIPPO) > 10)
table('brainseq' = abs(pcheck_both$Age.RegionHIPPO) > 20, 'span' = abs(pcheck_both$span_Age.RegionHIPPO) > 20)
table('brainseq' = abs(pcheck_both$Age.RegionHIPPO) > 100, 'span' = abs(pcheck_both$span_Age.RegionHIPPO) > 100)
table('brainseq' = abs(pcheck_both$Age.RegionHIPPO) > 1000, 'span' = abs(pcheck_both$span_Age.RegionHIPPO) > 1000)

weird <- which(abs(pcheck_both$Age.RegionHIPPO) > 10 | abs(pcheck_both$span_Age.RegionHIPPO) > 10)
table(pcheck_both$type[weird])

summary(c(pcheck_both$Age.RegionHIPPO[pcheck_both$type == 'jxn'], pcheck_both$span_Age.RegionHIPPO[pcheck_both$type == 'jxn']))
max(abs(c(pcheck_both$Age.RegionHIPPO[pcheck_both$type == 'jxn'], pcheck_both$span_Age.RegionHIPPO[pcheck_both$type == 'jxn'])))

table(pcheck_both$type[which(abs(pcheck_both$Age.RegionHIPPO) > 16 | abs(pcheck_both$span_Age.RegionHIPPO) > 16)])

png('pdf/compare_with_span_logFC.png', type = 'cairo')
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


pdf('pdf/compare_with_span_logFC_density.pdf', useDingbats = FALSE)
ggplot(pcheck_both, aes(x = logFC, y = span_logFC)) +
    facet_grid(age ~ type) + ylab('BrainSpan log FC') +
    xlab('BrainSeq log FC') + stat_density2d() + xlim(-16, 16) + ylim(-16, 16) +
    geom_smooth(method=lm, se=FALSE)
dev.off()

pdf('pdf/compare_with_span_logFC_density_noTx.pdf', useDingbats = FALSE)
ggplot(subset(pcheck_both, type != 'tx'), aes(x = logFC, y = span_logFC)) +
    facet_grid(age ~ type) + ylab('BrainSpan log FC') +
    xlab('BrainSeq log FC') + stat_density2d() + xlim(-7.5, 7.5) +
    ylim(-7.5, 7.5) +
    geom_smooth(method=lm, se=FALSE)
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

ggplot(subset(rep_span, pvar == 'P.Bonf' & type != 'tx'), aes(x = factor(paste0('p<', cutoff), paste0('p<', c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001))), y = replicated / number_de)) + facet_grid(type ~ term_clean) + ylab('Replication rate') + xlab('p-threshold') + geom_point() + theme_bw(base_size = 18)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(c(0, 1))

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

if(!file.exists('rda/go_case_genes_fdr1.Rdata')) {
    system.time( go_case_genes_fdr1 <- run_go(case_genes_fdr1[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_case_genes_fdr1.Rdata'))
    save(go_case_genes_fdr1, file = 'rda/go_case_genes_fdr1.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_case_genes_fdr1.Rdata'))
    load('rda/go_case_genes_fdr1.Rdata', verbose = TRUE)
}
sapply(go_case_genes_fdr1, class)

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

if(!file.exists('rda/go_me_genes_fdr05.Rdata')) {
    system.time( go_me_genes_fdr05 <- run_go(me_genes_fdr05[c('gene', 'exon', 'jxn')]) )
    message(paste(Sys.time(), 'saving rda/go_me_genes_fdr05.Rdata'))
    save(go_me_genes_fdr05, file = 'rda/go_me_genes_fdr05.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_me_genes_fdr05.Rdata'))
    load('rda/go_me_genes_fdr05.Rdata', verbose = TRUE)
}
sapply(go_me_genes_fdr05, class)


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


## Volcano plots



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


## Re-loading
f <- dir('rda', full.names = TRUE)
f <- f[!grepl('limma', f)]
for(ff in f) load(ff, verbose = TRUE)
rm(ff, f)

