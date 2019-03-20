library('clusterProfiler')
library('gplots')
library('GenomicRanges')
library('devtools')
library('VennDiagram')
library('RColorBrewer')
library('ggplot2')
library('jaffelab')
library('limma')
library('edgeR')
library('SummarizedExperiment')

dir.create('pdf', showWarnings = FALSE)
dir.create('rdas', showWarnings = FALSE)

## Load case-control results
files <- c(
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV_matchHIPPO.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV_matchDLPFC.rda'
)


outGene0 <- lapply(files, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(outGene0)
})

names(outGene0) <- c('HIPPO_matchQSV', 'DLPFC_matchQSV')



## Get number of DE genes at different FDR cutoffs
n_de <- do.call(rbind, lapply(c(0.05, 0.1, 0.15, 0.2), function(cut) {
    xx <- sapply(outGene0, function(x) {
        table(factor(x$adj.P.Val < cut, levels = c('FALSE', 'TRUE')))
    })
    cbind(xx, cutoff = cut)
}))

options(width = 200)
n_de
#       HIPPO_matchQSV DLPFC_matchQSV cutoff
# FALSE          24589          23568   0.05
# TRUE              63           1084   0.05
# FALSE          24518          22813   0.10
# TRUE             134           1839   0.10
# FALSE          24404          22127   0.15
# TRUE             248           2525   0.15
# FALSE          24216          21434   0.20
# TRUE             436           3218   0.20

data.frame(colnames(n_de))
#   colnames.n_de.
# 1 HIPPO_matchQSV
# 2 DLPFC_matchQSV
# 3         cutoff

n_de_sign <- do.call(rbind, lapply(c(0.05, 0.1, 0.15, 0.2), function(cut) {
    xx <- lapply(outGene0, function(x) {
        y <- table(factor(x$adj.P.Val < cut, levels = c('FALSE', 'TRUE')), factor(sign(x$logFC), levels = c(-1, 0, 1)))
        data.frame(de_status = rownames(y), n = as.vector(y), sign = rep(colnames(y), each = 2), cutoff = cut)
    })
    for(i in seq_len(length(xx))) { xx[[i]]$model = names(xx)[i] }
    names(xx) <- NULL
    do.call(rbind, xx)
}))
n_de_sign$group <- ifelse(n_de_sign$sign == 0, 'none', ifelse(n_de_sign$sign == -1, 'Control', 'SCZD'))
n_de_sign
# 1      FALSE 13136   -1   0.05 HIPPO_matchQSV Control
# 2       TRUE    31   -1   0.05 HIPPO_matchQSV Control
# 3      FALSE     0    0   0.05 HIPPO_matchQSV    none
# 4       TRUE     0    0   0.05 HIPPO_matchQSV    none
# 5      FALSE 11453    1   0.05 HIPPO_matchQSV    SCZD
# 6       TRUE    32    1   0.05 HIPPO_matchQSV    SCZD
# 7      FALSE 11969   -1   0.05 DLPFC_matchQSV Control
# 8       TRUE   507   -1   0.05 DLPFC_matchQSV Control
# 9      FALSE     0    0   0.05 DLPFC_matchQSV    none
# 10      TRUE     0    0   0.05 DLPFC_matchQSV    none
# 11     FALSE 11599    1   0.05 DLPFC_matchQSV    SCZD
# 12      TRUE   577    1   0.05 DLPFC_matchQSV    SCZD
# 13     FALSE 13110   -1   0.10 HIPPO_matchQSV Control
# 14      TRUE    57   -1   0.10 HIPPO_matchQSV Control
# 15     FALSE     0    0   0.10 HIPPO_matchQSV    none
# 16      TRUE     0    0   0.10 HIPPO_matchQSV    none
# 17     FALSE 11408    1   0.10 HIPPO_matchQSV    SCZD
# 18      TRUE    77    1   0.10 HIPPO_matchQSV    SCZD
# 19     FALSE 11606   -1   0.10 DLPFC_matchQSV Control
# 20      TRUE   870   -1   0.10 DLPFC_matchQSV Control
# 21     FALSE     0    0   0.10 DLPFC_matchQSV    none
# 22      TRUE     0    0   0.10 DLPFC_matchQSV    none
# 23     FALSE 11207    1   0.10 DLPFC_matchQSV    SCZD
# 24      TRUE   969    1   0.10 DLPFC_matchQSV    SCZD
# 25     FALSE 13063   -1   0.15 HIPPO_matchQSV Control
# 26      TRUE   104   -1   0.15 HIPPO_matchQSV Control
# 27     FALSE     0    0   0.15 HIPPO_matchQSV    none
# 28      TRUE     0    0   0.15 HIPPO_matchQSV    none
# 29     FALSE 11341    1   0.15 HIPPO_matchQSV    SCZD
# 30      TRUE   144    1   0.15 HIPPO_matchQSV    SCZD
# 31     FALSE 11258   -1   0.15 DLPFC_matchQSV Control
# 32      TRUE  1218   -1   0.15 DLPFC_matchQSV Control
# 33     FALSE     0    0   0.15 DLPFC_matchQSV    none
# 34      TRUE     0    0   0.15 DLPFC_matchQSV    none
# 35     FALSE 10869    1   0.15 DLPFC_matchQSV    SCZD
# 36      TRUE  1307    1   0.15 DLPFC_matchQSV    SCZD
# 37     FALSE 12972   -1   0.20 HIPPO_matchQSV Control
# 38      TRUE   195   -1   0.20 HIPPO_matchQSV Control
# 39     FALSE     0    0   0.20 HIPPO_matchQSV    none
# 40      TRUE     0    0   0.20 HIPPO_matchQSV    none
# 41     FALSE 11244    1   0.20 HIPPO_matchQSV    SCZD
# 42      TRUE   241    1   0.20 HIPPO_matchQSV    SCZD
# 43     FALSE 10877   -1   0.20 DLPFC_matchQSV Control
# 44      TRUE  1599   -1   0.20 DLPFC_matchQSV Control
# 45     FALSE     0    0   0.20 DLPFC_matchQSV    none
# 46      TRUE     0    0   0.20 DLPFC_matchQSV    none
# 47     FALSE 10557    1   0.20 DLPFC_matchQSV    SCZD
# 48      TRUE  1619    1   0.20 DLPFC_matchQSV    SCZD

pdf('pdf/n_de_sign.pdf', useDingbats = FALSE, width = 10, height = 18)
ggplot(n_de_sign, aes(x = de_status, y = n, fill = group)) + geom_bar(stat = 'identity', width = 0.5, position = 'dodge') + facet_grid(model ~ cutoff) + theme_bw(base_size = 18)
ggplot(subset(n_de_sign, de_status == 'TRUE'), aes(x = de_status, y = n, fill = group)) + geom_bar(stat = 'identity', width = 0.5, position = 'dodge') + facet_grid(model ~ cutoff) + theme_bw(base_size = 18)
ggplot(subset(n_de_sign, de_status == 'TRUE' & sign != 0), aes(x = de_status, y = n, fill = group)) + geom_bar(stat = 'identity', width = 0.5, position = 'dodge') + facet_grid(model ~ cutoff) + theme_bw(base_size = 18)
dev.off()

## Explore number of DE genes across models
make_venn <- function(i, txt = 'HIPPO_', cut = 0.1) {
    vinfo <- lapply(outGene0[i], function(x) {
        x$ensemblID[x$adj.P.Val < cut]
    })
    names(vinfo) <- gsub('_match', '.gene', gsub(txt, '', names(vinfo)))
    venn(vinfo) + title(paste('FDR cutoff:', cut))
}

make_venn2 <- function(i, txt = 'QSV|IPPO|LPFC') {
    make_venn(i, txt = txt, cut = 0.05)
    make_venn(i, txt = txt)
}

data.frame(name = names(outGene0))
#             name
# 1 HIPPO_matchQSV
# 2 DLPFC_matchQSV

pdf('pdf/venn_across_models.pdf', useDingbats = FALSE)
make_venn2(1:2)
dev.off()

## DE genes by sign
de_genes_sign <- mapply(function(x, cut, sign) {
    x$ensemblID[x$adj.P.Val < cut & sign(x$logFC) == sign]
}, outGene0[rep(1:2, each = 2)], rep(c(0.2, 0.1), each = 2), sign = c(-1, 1))
names(de_genes_sign) <- paste0(gsub('HIPPO', 'H', gsub('DLPFC', 'D', gsub('_matchQSV', '.gene', names(de_genes_sign)))), '_', c('Control', 'SCZD'))

xx <- as.data.frame(sapply(de_genes_sign, function(x) length(unique(x[!is.na(x)]))))
colnames(xx) <- 'unique_gene_id'
xx$i <- seq_len(nrow(xx))
xx$n_noNA <- sapply(de_genes_sign, function(x) length(x[!is.na(x)]))
xx
#                unique_gene_id i n_noNA
# H.gene_Control            195 1    195
# H.gene_SCZD               241 2    241
# D.gene_Control            870 3    870
# D.gene_SCZD               969 4    969


## Top DE genes by sign
de_genes_sign_top <- mapply(function(x, sign) {
    x <- x[sign(x$logFC) == sign, ]
    head(x$ensemblID[order(x$adj.P.Val, decreasing = FALSE)], 400)
}, outGene0[c(1, 1, 2, 2)], sign = rep(c(-1, 1), 2), SIMPLIFY = FALSE)
names(de_genes_sign_top) <- names(de_genes_sign)
as.data.frame(sapply(de_genes_sign_top, length))

## Pretty venn code
#venn_cols <- brewer.pal('Set1', n = 4)
venn_cols <- c('skyblue3', 'dark orange', 'red', 'purple')
names(venn_cols) <- paste0(rep(c('DLPFC', 'HIPPO'), each = 2), '_', c('control', 'SCZD'))
make_venn_pretty <- function(genes, title = 'DE features grouped by gene id') {
    genes <- lapply(genes, function(x) x[!is.na(x)])
    v <- venn.diagram(genes, filename = NULL,
        main = title,
        col = 'transparent', fill = venn_cols[seq_len(length(genes))],
        alpha = 0.5, margin = 0,
        main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
        cat.col = venn_cols[seq_len(length(genes))], scaled = FALSE, cat.pos = 0)
    grid.newpage()
    grid.draw(v)
}

pdf('pdf/venn_de_genes_by_sign.pdf', useDingbats = FALSE)
## gene
make_venn_pretty(de_genes_sign[1:4], 'Gene DLPFC FDR10%, HIPPO FDR20%')
make_venn_pretty(list('HIPPO' = de_genes_sign[[1]], 'DLPFC' = de_genes_sign[[3]]), 'Gene Control (DLPFC FDR10%, HIPPO FDR20%)')
make_venn_pretty(list('HIPPO' = de_genes_sign[[2]], 'DLPFC' = de_genes_sign[[4]]), 'Gene SCZD (DLPFC FDR10%, HIPPO FDR20%)')

## top 200 gene
make_venn_pretty(lapply(de_genes_sign_top[1:4], head, n = 200), 'gene: top 200 in each group')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[1]], 200), 'DLPFC' = head(de_genes_sign_top[[3]], 200)), 'gene: Control (top 200)')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[2]], 200), 'DLPFC' = head(de_genes_sign_top[[4]], 200)), 'gene: SCZD (top 200)')

make_venn_pretty(lapply(de_genes_sign_top[1:4], head, n = 150), 'gene: top 150 in each group')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[1]], 150), 'DLPFC' = head(de_genes_sign_top[[3]], 150)), 'gene: Control (top 150)')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[2]], 150), 'DLPFC' = head(de_genes_sign_top[[4]], 150)), 'gene: SCZD (top 150)')

make_venn_pretty(lapply(de_genes_sign_top[1:4], head, n = 100), 'gene: top 100 in each group')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[1]], 100), 'DLPFC' = head(de_genes_sign_top[[3]], 100)), 'gene: Control (top 100)')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[2]], 100), 'DLPFC' = head(de_genes_sign_top[[4]], 100)), 'gene: SCZD (top 100)')

make_venn_pretty(lapply(de_genes_sign_top[1:4], head, n = 50), 'gene: top 50 in each group')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[1]], 50), 'DLPFC' = head(de_genes_sign_top[[3]], 50)), 'gene: Control (top 50)')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[2]], 50), 'DLPFC' = head(de_genes_sign_top[[4]], 50)), 'gene: SCZD (top 50)')
dev.off()
system('rm VennDiagram*.log')


## Create outFeat just for the universe
## so the universe matches exactly was used in the original analysis
outFeat <- lapply(c('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda'), function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    outTx$ensemblID <- gsub('\\..*', '', outTx$gene_id)
    return(list('exon' = outExon, 'jxn' = outJxn, 'tx' = outTx))
})
names(outFeat) <- c('DLPFC')

## Gene ontology analysis
uni <- c(outGene0[[1]]$ensemblID, outFeat[[1]][[1]]$ensemblID, outFeat[[1]][[2]]$ensemblID)
uni <- unique(uni[!is.na(uni)])
stopifnot(length(uni) == 27514)

run_go <- function(genes, ont = c('BP', 'MF', 'CC')) {
    ## Change to ENSEMBL ids and remove NAs
    genes_ens <- lapply(lapply(genes, function(x) { gsub('\\..*', '', x) }), function(y) y[!is.na(y)])

    #genes_venn <- venn(genes_ens, show.plot = FALSE)

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


## DE genes
if(!file.exists('rdas/go_de_genes.Rdata')) {
    system.time( go_de_genes <- run_go(de_genes_sign) )
    message(paste(Sys.time(), 'saving rdas/go_de_genes.Rdata'))
    save(go_de_genes, file = 'rdas/go_de_genes.Rdata')
} else {
    message(paste(Sys.time(), 'loading rdas/go_de_genes.Rdata'))
    load('rdas/go_de_genes.Rdata', verbose = TRUE)
}
sapply(go_de_genes, class)

if(!file.exists('rdas/go_de_genes_top.Rdata')) {
    go_de_genes_top <- lapply(c(50, 100, 150, 200), function(n) {
        run_go(lapply(de_genes_sign_top, head, n = n))
    })
    names(go_de_genes_top) <- c(50, 100, 150, 200)
    message(paste(Sys.time(), 'saving rdas/go_de_genes_top.Rdata'))
    save(go_de_genes_top, file = 'rdas/go_de_genes_top.Rdata')
} else {
    message(paste(Sys.time(), 'loading rdas/go_de_genes_top.Rdata'))
    load('rdas/go_de_genes_top.Rdata', verbose = TRUE)
}
lapply(go_de_genes_top, function(x) sapply(x, class))


simplify_go <- function(x) {
    #gsub('QSV|IPPO|LPFC', '', x)
    #gsub('_matchQSV', '', x)
    gsub('IPPO|LPFC|ontrol|chizo|CZD|xon_|xn_|x_|ene_|\\.', '', x)
}

plot_go <- function(go_cluster, cat = 10) {
    lapply(names(go_cluster), function(bp) {
        go <- go_cluster[[bp]]
        if(is.null(go)) {
            message(paste(Sys.time(), 'found no results for', bp))
            return(NULL)
        }

        ## Simplify names
        go@compareClusterResult$Cluster <- factor(simplify_go(go@compareClusterResult$Cluster), levels = simplify_go(levels(go@compareClusterResult$Cluster)))
        names(go@geneClusters) <- simplify_go(names(go@geneClusters))

        print(plot(go, title = paste('ontology:', bp), font.size = 18, showCategory = cat, includeAll = TRUE))
        return(NULL)
    })
}


pdf('pdf/go_de_genes.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_go(go_de_genes)
dev.off()

pdf('pdf/go_all_de_genes.pdf', width = 22, height = 70, useDingbats = FALSE)
plot_go(go_de_genes, cat = NULL)
dev.off()

for(i in names(go_de_genes_top)) {
    pdf(paste0('pdf/go_de_genes_top', i, '.pdf'), width = 22, height = 9, useDingbats = FALSE)
    plot_go(go_de_genes_top[[i]])
    dev.off()

    pdf(paste0('pdf/go_all_de_genes_top', i, '.pdf'), width = 22, height = 25, useDingbats = FALSE)
    plot_go(go_de_genes_top[[i]], cat = NULL)
    dev.off()
}

run_gse <-  function(region, ont = c('BP', 'MF', 'CC')) {

    genes <- outGene0[[paste0(region, '_matchQSV')]]$t
    names(genes) <- outGene0[[paste0(region, '_matchQSV')]]$ensemblID
    genes <- genes[order(genes, decreasing = TRUE)]

    ## Run gseGO analysis
    go_cluster <- lapply(ont, function(bp) {
        message(paste(Sys.time(), 'running gseGO analysis for', bp))
        tryCatch(gseGO(genes, OrgDb = 'org.Hs.eg.db',
            ont = bp, pAdjustMethod = "BH",
            keyType = 'ENSEMBL', verbose = TRUE),
            error = function(e) { return(NULL) })
    })
    names(go_cluster) <- ont
    
    
    
    genes_tab <- bitr(names(genes), fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
    genes_tab <- genes_tab[!duplicated(genes_tab$ENTREZID), ]
        
    genes_m <- match(genes_tab$ENSEMBL, names(genes))
    genes_ncbi <- genes[genes_m]
    names(genes_ncbi) <- genes_tab$ENTREZID

    go_cluster$KEGG <- tryCatch(gseKEGG(genes_ncbi, organism = 'hsa',
        pAdjustMethod = "BH",
        keyType = 'ncbi-geneid', verbose = TRUE),
        error = function(e) { return(NULL) })
        
    return(go_cluster)
}


system.time( gse_hippo <- run_gse('HIPPO') )
system.time( gse_dlpfc <- run_gse('DLPFC') )
save(gse_hippo, gse_dlpfc, file = 'rdas/gse.Rdata')

plot_gse <- function(gse) {
    lapply(names(gse), function(bp) {
        go <- gse[[bp]]
        if(is.null(go)) {
            message(paste(Sys.time(), 'found no results for', bp))
            return(NULL)
        }

        print(dotplot(go, title = paste('ontology:', bp), font.size = 18))
        return(NULL)
    })
}


pdf('pdf/gse_hippo.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_gse(gse_hippo)
dev.off()

pdf('pdf/gse_dlpfc.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_gse(gse_dlpfc)
dev.off()

## For DE_qual plots
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
                     function(x)
                       rgb(x[1], x[2], x[3], alpha=alpha))
}

## Load BrainSeq Phase 1 and Common Mind results
load('/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/caseControl/rdas/expressed_de_features.rda', verbose = TRUE)

## Add BSP1 and CMC
prev <- list(
    'BSP1' = data.frame(
        ensemblID = names(outStatsExprs$Gene),
        adj.P.Val = outStatsExprs$Gene$fdr_qsva,
        logFC = outStatsExprs$Gene$log2FC_qsva,
        t = outStatsExprs$Gene$tstat_qsva
        ),
    'CMC' = data.frame(
        ensemblID = names(outStatsExprs$Gene),
        adj.P.Val = p.adjust(outStatsExprs$Gene$CMC_pval_qsva, method = 'fdr'),
        logFC = outStatsExprs$Gene$CMC_log2FC_qsva,
        t = outStatsExprs$Gene$CMC_tstat_qsva
        )
)




## Scatter of logFC
comp_log <- function(x, y, xlab, ylab, var = 'logFC', de = FALSE, n = 150, onlyx = FALSE) {
    if(de) {
        if(!onlyx) {
            common <- unique(c(
                head(x$ensemblID[order(x$adj.P.Val, decreasing = FALSE)], n),
                head(y$ensemblID[order(y$adj.P.Val, decreasing = FALSE)], n)
            ))
        } else {
            common <- unique(c(
                head(x$ensemblID[order(x$adj.P.Val, decreasing = FALSE)], n)
            ))
        }
        
    } else {
        common <- intersect(x$ensemblID, y$ensemblID)
    }
    x <- x[match(common, x$ensemblID), ]
    y <- y[match(common, y$ensemblID), ]
    corr = signif(cor(x[, var], y[, var], use = 'pairwise.complete.obs'), 3)

    plot(x = x[, var], y = y[, var],
         xlab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), xlab),
         ylab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), ylab),
         col = add.alpha('black', ifelse(de, 1/2, 1/10)), pch = 16)
    legend('topleft', legend = paste('r =', corr))
    # lines(loess.smooth(y = y[, var], x = x[, var]), col = 'red')
    # abline(lm(y[, var] ~ x[, var]), col = 'blue')
    abline(h = 0, col = 'grey20')
    abline(v = 0, col = 'grey20')
}

pdf('pdf/scatter_models.pdf', useDingbats = FALSE)
comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO', 'DLPFC')
comp_log(outGene0[[1]], outGene[[1]], 'HIPPO', 'BSP1')
comp_log(outGene0[[1]], outGene[[2]], 'HIPPO', 'CMC')
comp_log(outGene0[[2]], outGene[[1]], 'DLPFC', 'BSP1')
comp_log(outGene0[[2]], outGene[[2]], 'DLPFC', 'CMC')

comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO', 'HIPPO with qSVA')

comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO', 'DLPFC', var = 't')
comp_log(outGene0[[1]], outGene[[1]], 'HIPPO', 'BSP1', var = 't')
comp_log(outGene0[[1]], outGene[[2]], 'HIPPO', 'CMC', var = 't')
comp_log(outGene0[[2]], outGene[[1]], 'DLPFC', 'BSP1', var = 't')
comp_log(outGene0[[2]], outGene[[2]], 'DLPFC', 'CMC', var = 't')
dev.off()

pdf('pdf/scatter_models_top150de.pdf', useDingbats = FALSE)
comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO', 'DLPFC', de = TRUE)
comp_log(outGene0[[1]], outGene[[1]], 'HIPPO', 'BSP1', de = TRUE)
comp_log(outGene0[[1]], outGene[[2]], 'HIPPO', 'CMC', de = TRUE)
comp_log(outGene0[[2]], outGene[[1]], 'DLPFC', 'BSP1', de = TRUE)
comp_log(outGene0[[2]], outGene[[2]], 'DLPFC', 'CMC', de = TRUE)

comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO', 'DLPFC', var = 't', de = TRUE)
comp_log(outGene0[[1]], outGene[[1]], 'HIPPO', 'BSP1', var = 't', de = TRUE)
comp_log(outGene0[[1]], outGene[[2]], 'HIPPO', 'CMC', var = 't', de = TRUE)
comp_log(outGene0[[2]], outGene[[1]], 'DLPFC', 'BSP1', var = 't', de = TRUE)
comp_log(outGene0[[2]], outGene[[2]], 'DLPFC', 'CMC', var = 't', de = TRUE)
dev.off()

pdf('pdf/scatter_models_top400de.pdf', useDingbats = FALSE)
comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO', 'DLPFC', de = TRUE, n = 400)
comp_log(outGene0[[1]], outGene[[1]], 'HIPPO', 'BSP1', de = TRUE, n = 400)
comp_log(outGene0[[1]], outGene[[2]], 'HIPPO', 'CMC', de = TRUE, n = 400)
comp_log(outGene0[[2]], outGene[[1]], 'DLPFC', 'BSP1', de = TRUE, n = 400)
comp_log(outGene0[[2]], outGene[[2]], 'DLPFC', 'CMC', de = TRUE, n = 400)

comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO', 'DLPFC', var = 't', de = TRUE, n = 400)
comp_log(outGene0[[1]], outGene[[1]], 'HIPPO', 'BSP1', var = 't', de = TRUE, n = 400)
comp_log(outGene0[[1]], outGene[[2]], 'HIPPO', 'CMC', var = 't', de = TRUE, n = 400)
comp_log(outGene0[[2]], outGene[[1]], 'DLPFC', 'BSP1', var = 't', de = TRUE, n = 400)
comp_log(outGene0[[2]], outGene[[2]], 'DLPFC', 'CMC', var = 't', de = TRUE, n = 400)
dev.off()


pdf('pdf/scatter_models_top400de_onlyX.pdf', useDingbats = FALSE)
comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO', 'DLPFC', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[1]], outGene[[1]], 'HIPPO', 'BSP1', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[1]], outGene[[2]], 'HIPPO', 'CMC', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[2]], outGene[[1]], 'DLPFC', 'BSP1', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[2]], outGene[[2]], 'DLPFC', 'CMC', de = TRUE, n = 400, onlyx = TRUE)

comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO', 'DLPFC', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[1]], outGene[[1]], 'HIPPO', 'BSP1', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[1]], outGene[[2]], 'HIPPO', 'CMC', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[2]], outGene[[1]], 'DLPFC', 'BSP1', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[2]], outGene[[2]], 'DLPFC', 'CMC', var = 't', de = TRUE, n = 400, onlyx = TRUE)
dev.off()




## Load expression data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata', verbose = TRUE)
rse_gene_HIPPO <- rse_gene[, keepIndex]
mod_HIPPO <- mod
modQsva_HIPPO <- modQsva
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_DLPFC.Rdata', verbose = TRUE)
rse_gene_DLPFC <- rse_gene[, keepIndex]
mod_DLPFC <- mod
modQsva_DLPFC <- modQsva

## Get normalized expression and clean it
cleaned <- mapply(function(expr, model) {
    dge = DGEList(counts = assays(expr)$counts,
	genes = rowData(expr))
    #calculate library-size adjustment
    dge = calcNormFactors(dge)
    vGene = voom(dge, model, plot=FALSE)
    list('norm' = vGene$E, 'cleaned' = cleaningY(vGene$E, model, P = 2))
},
    list('HIPPO' = rse_gene_HIPPO, 'DLPFC' = rse_gene_DLPFC),
    list('HIPPO' = mod_HIPPO, 'DLPFC' = mod_DLPFC),
    SIMPLIFY = FALSE
)

## Adapted some functions from
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/region_specific/explore_reg_specific.R
# and
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/region_specific/explore_reg_specific_top.R
get_ylab <- function(type, cleaned = FALSE) {
    if(type %in% c('gene', 'exon', 'jxn')) {
        res <- 'log2(CPM + 0.5)'
    } else {
        res <- 'log2(TPM + 0.5)'
    }

    if(cleaned) res <- paste(res, '- covariate effects removed')
    return(res)
}

get_name <- function(region, group) {
    paste0(substr(region, 1, 1), '.gene_', ifelse(group == 'schizo', 'SCZD', 'Control'))
}

get_main <- function(i, region, group = 'schizo', type = 'gene') {
    if(type %in% c('gene', 'exon')) {
        var <- 'gencodeID'
        vars <- 'Symbol'
    } else if (type == 'jxn') {
        var <- 'gencodeGeneID'
        vars <- 'Symbol'
    } else {
        var <- 'gene_id'
        vars <- 'gene_name'
    }

    rse <- if(region == 'HIPPO') rse_gene_HIPPO else rse_gene_DLPFC

    topnow <- outGene0[[paste0(region, '_matchQSV')]]

    j <- which(topnow$ensemblID == de_genes_sign_top[[get_name(region, group)]][i])
    k <- which(names(rowRanges(rse)) == rownames(topnow)[j])

    paste(if(type != 'jxn') mcols(rowRanges(rse))[, var][k] else rownames(topnow)[j],
          if(is.na(mcols(rowRanges(rse))[, vars][k])) '' else mcols(rowRanges(rse))[, vars][k],
          'FDR',
          signif(topnow$adj.P.Val[j], 3),
          ifelse(group == 'schizo', 'SCZD', 'Control'))
}

get_ylim_mult <- function(rang) {
    c(
        ifelse(sign(rang[1]) == 1, 0.95, 1.05),
        ifelse(sign(rang[2]) == 1, 1.05, 0.95)
    )
}


plot_top_sign <- function(i, region, group = 'schizo', normtype = 'norm') {
    g <- de_genes_sign_top[[get_name(region, group)]][i]
    j <- which(gsub('\\..*', '', rownames(cleaned[[region]][[normtype]])) == g)

    set.seed(20180426)
    dx <- if(region == 'HIPPO') colData(rse_gene_HIPPO)$Dx else colData(rse_gene_DLPFC)$Dx
    dx[dx == 'Schizo'] <- 'SCZD'
    dx <- factor(dx, levels = c('Control', 'SCZD'))

    y <- cleaned[[region]][[normtype]][j, ]
    boxplot(y ~ dx,
            ylab = get_ylab('gene', cleaned = normtype == 'cleaned'),
            main = get_main(i, region, group),
            col = c('orchid1', 'aquamarine1'),
            ylim = abs(range(y)) * get_ylim_mult(range(y)) * sign(range(y)),
            outline = FALSE)
    points(y ~ jitter(as.integer(dx), 1), pch = 21,
           bg = c('orchid4', 'aquamarine4')[as.integer(dx)])
}

for(reg in c('HIPPO', 'DLPFC')) {
    pdf(paste0('pdf/top_50each_', reg, '_norm.pdf'))
    for(i in 1:50) {
        plot_top_sign(i, reg)
        plot_top_sign(i, reg, 'control')
    }
    dev.off()
    pdf(paste0('pdf/top_50each_', reg, '_cleaned.pdf'))
    for(i in 1:50) {
        plot_top_sign(i, reg, normtype = 'cleaned')
        plot_top_sign(i, reg, 'control', normtype = 'cleaned')
    }
    dev.off()
}

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
