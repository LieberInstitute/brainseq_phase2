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

outGene <- lapply(files, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(outGene)
})
names(outGene) <- c('HIPPO_matchQSV', 'DLPFC_matchQSV')
## Make sure the order is correct
stopifnot(all(sapply(gsub('_.*', '', names(outGene)), grep, x = files) - 1:2 == 0))

## Not the same order as in
## https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/explore_case_control.R
## But that way I don't have to change the code that
## I started last week
outGene <- c(prev, outGene)


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
    
    ## Use colors by FDR < 0.05
    colgroup <- rep('none', nrow(x))
    colgroup[x$adj.P.Val < 0.05] <- 'x'
    colgroup[y$adj.P.Val < 0.05] <- 'y'
    colgroup[x$adj.P.Val < 0.05 & y$adj.P.Val < 0.05] <- 'xy'
    cols <- dplyr::case_when(
        colgroup == 'none' ~ add.alpha('black', ifelse(de, 1/2, 1/10)),
        colgroup == 'x' ~ add.alpha('magenta', ifelse(de, 1/2, 1/4)),
        colgroup == 'y' ~ add.alpha('magenta', ifelse(de, 1/2, 1/4)),
        colgroup == 'xy' ~ add.alpha('purple', 1/2)
    )

    par(cex.axis = 2, cex.lab = 2, mar = c(5, 5, 4, 2) + 0.1)
    plot(x = x[, var], y = y[, var],
         xlab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), xlab),
         ylab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), ylab),
         col = cols, pch = 16)
    legend('topleft', legend = paste('r =', corr), cex = 1.8)
    # lines(loess.smooth(y = y[, var], x = x[, var]), col = 'red')
    # abline(lm(y[, var] ~ x[, var]), col = 'blue')
    abline(h = 0, col = 'grey20')
    abline(v = 0, col = 'grey20')
}

pdf('pdf/scatter_models.pdf', useDingbats = FALSE)
# comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO no qSVA', 'DLPFC no qSVA')
# comp_log(outGene0[[1]], outGene[[1]], 'HIPPO no qSVA', 'BSP1')
# comp_log(outGene0[[1]], outGene[[2]], 'HIPPO no qSVA', 'CMC')
# comp_log(outGene0[[2]], outGene[[1]], 'DLPFC no qSVA', 'BSP1')
# comp_log(outGene0[[2]], outGene[[2]], 'DLPFC no qSVA', 'CMC')
#
# comp_log(outGene0[[1]], outGene[[3]], 'HIPPO no qSVA', 'HIPPO with qSVA')
# comp_log(outGene0[[2]], outGene[[4]], 'DLFPC no qSVA', 'DLPFC with qSVA')

comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO no qSVA', 'DLPFC no qSVA', var = 't')
comp_log(outGene0[[1]], outGene[[1]], 'HIPPO no qSVA', 'BSP1', var = 't')
comp_log(outGene0[[1]], outGene[[2]], 'HIPPO no qSVA', 'CMC', var = 't')
comp_log(outGene0[[2]], outGene[[1]], 'DLPFC no qSVA', 'BSP1', var = 't')
comp_log(outGene0[[2]], outGene[[2]], 'DLPFC no qSVA', 'CMC', var = 't')

comp_log(outGene0[[1]], outGene[[3]], 'HIPPO no qSVA', 'HIPPO with qSVA', var = 't')
comp_log(outGene0[[2]], outGene[[4]], 'DLFPC no qSVA', 'DLPFC with qSVA', var = 't')
dev.off()

pdf('pdf/scatter_models_top150de.pdf', useDingbats = FALSE)
# comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO no qSVA', 'DLPFC no qSVA', de = TRUE)
# comp_log(outGene0[[1]], outGene[[1]], 'HIPPO no qSVA', 'BSP1', de = TRUE)
# comp_log(outGene0[[1]], outGene[[2]], 'HIPPO no qSVA', 'CMC', de = TRUE)
# comp_log(outGene0[[2]], outGene[[1]], 'DLPFC no qSVA', 'BSP1', de = TRUE)
# comp_log(outGene0[[2]], outGene[[2]], 'DLPFC no qSVA', 'CMC', de = TRUE)
#
# comp_log(outGene0[[1]], outGene[[3]], 'HIPPO no qSVA', 'HIPPO with qSVA', de = TRUE)
# comp_log(outGene0[[2]], outGene[[4]], 'DLFPC no qSVA', 'DLPFC with qSVA', de = TRUE)

comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO no qSVA', 'DLPFC no qSVA', var = 't', de = TRUE)
comp_log(outGene0[[1]], outGene[[1]], 'HIPPO no qSVA', 'BSP1', var = 't', de = TRUE)
comp_log(outGene0[[1]], outGene[[2]], 'HIPPO no qSVA', 'CMC', var = 't', de = TRUE)
comp_log(outGene0[[2]], outGene[[1]], 'DLPFC no qSVA', 'BSP1', var = 't', de = TRUE)
comp_log(outGene0[[2]], outGene[[2]], 'DLPFC no qSVA', 'CMC', var = 't', de = TRUE)

comp_log(outGene0[[1]], outGene[[3]], 'HIPPO no qSVA', 'HIPPO with qSVA', var = 't', de = TRUE)
comp_log(outGene0[[2]], outGene[[4]], 'DLFPC no qSVA', 'DLPFC with qSVA', var = 't', de = TRUE)
dev.off()

pdf('pdf/scatter_models_top400de.pdf', useDingbats = FALSE)
# comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO no qSVA', 'DLPFC no qSVA', de = TRUE, n = 400)
# comp_log(outGene0[[1]], outGene[[1]], 'HIPPO no qSVA', 'BSP1', de = TRUE, n = 400)
# comp_log(outGene0[[1]], outGene[[2]], 'HIPPO no qSVA', 'CMC', de = TRUE, n = 400)
# comp_log(outGene0[[2]], outGene[[1]], 'DLPFC no qSVA', 'BSP1', de = TRUE, n = 400)
# comp_log(outGene0[[2]], outGene[[2]], 'DLPFC no qSVA', 'CMC', de = TRUE, n = 400)
#
# comp_log(outGene0[[1]], outGene[[3]], 'HIPPO no qSVA', 'HIPPO with qSVA', de = TRUE, n = 400)
# comp_log(outGene0[[2]], outGene[[4]], 'DLFPC no qSVA', 'DLPFC with qSVA', de = TRUE, n = 400)

comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO no qSVA', 'DLPFC no qSVA', var = 't', de = TRUE, n = 400)
comp_log(outGene0[[1]], outGene[[1]], 'HIPPO no qSVA', 'BSP1', var = 't', de = TRUE, n = 400)
comp_log(outGene0[[1]], outGene[[2]], 'HIPPO no qSVA', 'CMC', var = 't', de = TRUE, n = 400)
comp_log(outGene0[[2]], outGene[[1]], 'DLPFC no qSVA', 'BSP1', var = 't', de = TRUE, n = 400)
comp_log(outGene0[[2]], outGene[[2]], 'DLPFC no qSVA', 'CMC', var = 't', de = TRUE, n = 400)

comp_log(outGene0[[1]], outGene[[3]], 'HIPPO no qSVA', 'HIPPO with qSVA', var = 't', de = TRUE, n = 400)
comp_log(outGene0[[2]], outGene[[4]], 'DLFPC no qSVA', 'DLPFC with qSVA', var = 't', de = TRUE, n = 400)
dev.off()


pdf('pdf/scatter_models_top400de_onlyX.pdf', useDingbats = FALSE)
# comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO no qSVA', 'DLPFC no qSVA', de = TRUE, n = 400, onlyx = TRUE)
# comp_log(outGene0[[1]], outGene[[1]], 'HIPPO no qSVA', 'BSP1', de = TRUE, n = 400, onlyx = TRUE)
# comp_log(outGene0[[1]], outGene[[2]], 'HIPPO no qSVA', 'CMC', de = TRUE, n = 400, onlyx = TRUE)
# comp_log(outGene0[[2]], outGene[[1]], 'DLPFC no qSVA', 'BSP1', de = TRUE, n = 400, onlyx = TRUE)
# comp_log(outGene0[[2]], outGene[[2]], 'DLPFC no qSVA', 'CMC', de = TRUE, n = 400, onlyx = TRUE)
#
# comp_log(outGene0[[1]], outGene[[3]], 'HIPPO no qSVA', 'HIPPO with qSVA', de = TRUE, n = 400, onlyx = TRUE)
# comp_log(outGene0[[2]], outGene[[4]], 'DLFPC no qSVA', 'DLPFC with qSVA', de = TRUE, n = 400, onlyx = TRUE)

comp_log(outGene0[[1]], outGene0[[2]], 'HIPPO no qSVA', 'DLPFC no qSVA', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[1]], outGene[[1]], 'HIPPO no qSVA', 'BSP1', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[1]], outGene[[2]], 'HIPPO no qSVA', 'CMC', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[2]], outGene[[1]], 'DLPFC no qSVA', 'BSP1', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[2]], outGene[[2]], 'DLPFC no qSVA', 'CMC', var = 't', de = TRUE, n = 400, onlyx = TRUE)

comp_log(outGene0[[1]], outGene[[3]], 'HIPPO no qSVA', 'HIPPO with qSVA', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene0[[2]], outGene[[4]], 'DLFPC no qSVA', 'DLPFC with qSVA', var = 't', de = TRUE, n = 400, onlyx = TRUE)
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

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.5.1 Patched (2018-10-29 r75535)
#  os       Red Hat Enterprise Linux Server release 6.9 (Santiago)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-03-20
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  AnnotationDbi          1.44.0    2018-10-30 [1] Bioconductor
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  backports              1.1.3     2018-12-14 [2] CRAN (R 3.5.1)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.6    2019-02-10 [1] Bioconductor
#  bit                    1.1-14    2018-05-29 [2] CRAN (R 3.5.1)
#  bit64                  0.9-7     2017-05-08 [2] CRAN (R 3.5.0)
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  blob                   1.1.1     2018-03-25 [2] CRAN (R 3.5.0)
#  callr                  3.1.1     2018-12-21 [2] CRAN (R 3.5.1)
#  caTools                1.17.1.2  2019-03-06 [2] CRAN (R 3.5.1)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  clusterProfiler      * 3.10.1    2018-12-20 [1] Bioconductor
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  cowplot                0.9.4     2019-01-08 [1] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  data.table             1.12.0    2019-01-13 [1] CRAN (R 3.5.1)
#  DBI                    1.0.0     2018-05-02 [2] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  desc                   1.2.0     2018-05-01 [2] CRAN (R 3.5.1)
#  devtools             * 2.0.1     2018-10-26 [1] CRAN (R 3.5.1)
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  DO.db                  2.9       2018-05-03 [1] Bioconductor
#  DOSE                   3.8.2     2019-01-14 [1] Bioconductor
#  dplyr                  0.8.0.1   2019-02-15 [1] CRAN (R 3.5.1)
#  edgeR                * 3.24.3    2019-01-02 [1] Bioconductor
#  enrichplot             1.2.0     2018-10-30 [2] Bioconductor
#  europepmc              0.3       2018-04-20 [2] CRAN (R 3.5.1)
#  farver                 1.1.0     2018-11-20 [1] CRAN (R 3.5.1)
#  fastmatch              1.1-0     2017-01-28 [1] CRAN (R 3.5.0)
#  fgsea                  1.8.0     2018-10-30 [1] Bioconductor
#  formatR                1.6       2019-03-05 [1] CRAN (R 3.5.1)
#  fs                     1.2.6     2018-08-23 [2] CRAN (R 3.5.1)
#  futile.logger        * 1.4.3     2016-07-10 [1] CRAN (R 3.5.0)
#  futile.options         1.0.1     2018-04-20 [2] CRAN (R 3.5.0)
#  gdata                  2.18.0    2017-06-06 [2] CRAN (R 3.5.0)
#  GenomeInfoDb         * 1.18.2    2019-02-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggforce                0.2.1     2019-03-12 [2] CRAN (R 3.5.1)
#  ggplot2              * 3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  ggplotify              0.0.3     2018-08-03 [2] CRAN (R 3.5.1)
#  ggraph                 1.0.2     2018-07-07 [2] CRAN (R 3.5.1)
#  ggrepel                0.8.0     2018-05-09 [1] CRAN (R 3.5.0)
#  ggridges               0.5.1     2018-09-27 [1] CRAN (R 3.5.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
#  GO.db                  3.7.0     2018-11-02 [1] Bioconductor
#  GOSemSim               2.8.0     2018-10-30 [1] Bioconductor
#  gplots               * 3.0.1.1   2019-01-27 [1] CRAN (R 3.5.1)
#  gridExtra              2.3       2017-09-09 [2] CRAN (R 3.5.0)
#  gridGraphics           0.3-0     2018-04-03 [2] CRAN (R 3.5.1)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  gtools                 3.8.1     2018-06-26 [2] CRAN (R 3.5.1)
#  hms                    0.4.2     2018-03-10 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  httr                   1.4.0     2018-12-11 [1] CRAN (R 3.5.1)
#  igraph                 1.2.4     2019-02-13 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab             * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  KernSmooth             2.23-15   2015-06-29 [3] CRAN (R 3.5.1)
#  labeling               0.3       2014-08-23 [2] CRAN (R 3.5.0)
#  lambda.r               1.2.3     2018-05-17 [1] CRAN (R 3.5.0)
#  later                  0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma                * 3.38.3    2018-12-02 [1] Bioconductor
#  locfit                 1.5-9.1   2013-04-20 [2] CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  MASS                   7.3-51.1  2018-11-01 [3] CRAN (R 3.5.1)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  memoise                1.1.0     2017-04-21 [2] CRAN (R 3.5.0)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgbuild               1.0.2     2018-10-16 [2] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  pkgload                1.0.2     2018-10-29 [2] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  polyclip               1.10-0    2019-03-14 [2] CRAN (R 3.5.1)
#  prettyunits            1.0.2     2015-07-13 [1] CRAN (R 3.5.0)
#  processx               3.3.0     2019-03-10 [1] CRAN (R 3.5.1)
#  progress               1.2.0     2018-06-14 [1] CRAN (R 3.5.1)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  ps                     1.3.0     2018-12-21 [2] CRAN (R 3.5.1)
#  purrr                  0.3.1     2019-03-03 [2] CRAN (R 3.5.1)
#  qvalue                 2.14.1    2019-01-10 [1] Bioconductor
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer         * 1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  remotes                2.0.2     2018-10-30 [1] CRAN (R 3.5.1)
#  reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.5.0)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  rprojroot              1.3-2     2018-01-03 [2] CRAN (R 3.5.0)
#  RSQLite                2.1.1     2018-05-06 [2] CRAN (R 3.5.0)
#  rvcheck                0.1.3     2018-12-06 [1] CRAN (R 3.5.1)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr                  0.13      2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo            1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  stringi                1.4.3     2019-03-12 [2] CRAN (R 3.5.1)
#  stringr                1.4.0     2019-02-10 [1] CRAN (R 3.5.1)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  testthat               2.0.1     2018-10-13 [1] CRAN (R 3.5.1)
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyr                  0.8.3     2019-03-01 [2] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  triebeard              0.3.0     2016-08-04 [1] CRAN (R 3.5.0)
#  tweenr                 1.0.1     2018-12-14 [1] CRAN (R 3.5.1)
#  UpSetR                 1.3.3     2017-03-21 [1] CRAN (R 3.5.0)
#  urltools               1.7.2     2019-02-04 [1] CRAN (R 3.5.1)
#  usethis              * 1.4.0     2018-08-14 [2] CRAN (R 3.5.1)
#  VennDiagram          * 1.6.20    2018-03-28 [1] CRAN (R 3.5.0)
#  viridis                0.5.1     2018-03-29 [2] CRAN (R 3.5.0)
#  viridisLite            0.3.0     2018-02-01 [2] CRAN (R 3.5.0)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.5       2019-02-20 [1] CRAN (R 3.5.1)
#  xml2                   1.2.0     2018-01-24 [2] CRAN (R 3.5.0)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
