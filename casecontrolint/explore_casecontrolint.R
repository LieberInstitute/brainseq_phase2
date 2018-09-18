library('limma')
library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('ggplot2')
library('gplots')
library('VennDiagram')
library('RColorBrewer')
library('clusterProfiler')

dir.create('rda', showWarnings = FALSE)
dir.create('pdf', showWarnings = FALSE)

## Load qSVs (n=712)
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold.Rdata', verbose = TRUE)

## Load expr data: function
load_foo <- function(type) {
    load_file <- file.path(
        '/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff',
        paste0('rse_', type, '.Rdata'))
    stopifnot(file.exists(load_file))
    message(paste(Sys.time(), 'loading', load_file))
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

    ## Keep samples in qSVs
    rse <- rse[, keepIndex]

    ## Set as factor
    colData(rse)$Region <- relevel(factor(colData(rse)$Region), 'DLPFC')
    colData(rse)$Race <- relevel(factor(colData(rse)$Race), ref = 'CAUC')
    colData(rse)$Sex <- relevel(factor(colData(rse)$Sex), ref = 'F')

    ## Add means
    colData(rse)$mean_mitoRate <- mean(colData(rse)$mitoRate)
    colData(rse)$mean_totalAssignedGene <- mean(colData(rse)$totalAssignedGene)
    colData(rse)$mean_rRNA_rate <- mean(colData(rse)$rRNA_rate)
    colData(rse)$mean_RIN <- mean(colData(rse)$RIN)

    return(rse)
}

## Load expr data
# types <- c('gene', 'exon', 'jxn', 'tx')
# rses <- lapply(types, load_foo)
# names(rses) <- types
rse <- load_foo('gene')

## Load BrainSeq model results
if(!file.exists('rda/raw.Rdata')) {
    raw <- lapply(c('gene', 'exon', 'jxn', 'tx'), function(type) {
        f <- paste0('rda/limma_casecontrol_interaction_', type, '.Rdata')
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
for(i in seq_len(length(top))) colnames(top[[i]])[1:6] <- paste0('interaction_', colnames(top[[i]])[1:6])
fit <- lapply(raw, '[[', 'fit')
exprsNorm <- lapply(raw, '[[', 'exprsNorm')

sapply(top, function(x) {
    mean(x$interaction_adj.P.Val < 0.05) * 100
})
#        gene        exon         jxn          tx
# 0.450267727 0.200968786 0.059223167 0.007548635

sapply(top, function(x) {
    sum(x$interaction_adj.P.Val < 0.05)
})
# gene exon  jxn   tx
#  111  797  176    7

sapply(top, dim)
#       gene   exon    jxn    tx
# [1,] 24652 396579 297181 92732
# [2,]     7      7      7     7

## Load SCZD case-control results
message(paste(Sys.time(), 'loading SCZD case-control results'))
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/correlation/rda/out_info.Rdata', verbose = TRUE)

## Check that the data is in the same order (should be)
stopifnot(identical(rownames(outGene[['HIPPO_matchQSV']]), rownames(top$gene)))
stopifnot(identical(rownames(outGene[['DLPFC_matchQSV']]), rownames(top$gene)))

## Exons almost match, here I dropped a few
m <- match(rownames(top$exon), rownames(outFeat[['DLPFC']]$exon))
stopifnot(any(!is.na(m)))
outFeat[['DLPFC']]$exon <- outFeat[['DLPFC']]$exon[m, ]
m <- match(rownames(top$exon), rownames(outFeat[['HIPPO']]$exon))
stopifnot(any(!is.na(m)))
outFeat[['HIPPO']]$exon <- outFeat[['HIPPO']]$exon[m, ]

## now they do
stopifnot(identical(rownames(outFeat[['DLPFC']]$exon), rownames(top$exon)))
stopifnot(identical(rownames(outFeat[['HIPPO']]$exon), rownames(top$exon)))

## jxn and tx match
stopifnot(identical(rownames(outFeat[['DLPFC']]$jxn), rownames(top$jxn)))
stopifnot(identical(rownames(outFeat[['DLPFC']]$tx), rownames(top$tx)))
stopifnot(identical(rownames(outFeat[['HIPPO']]$jxn), rownames(top$jxn)))
stopifnot(identical(rownames(outFeat[['HIPPO']]$tx), rownames(top$tx)))


## Paste together with SCZD case-controls results
message(paste(Sys.time(), 'pasting with SCZD case-control results'))

outGene[['HIPPO_matchQSV']] <- cbind(outGene[['HIPPO_matchQSV']], top$gene)
outGene[['DLPFC_matchQSV']] <- cbind(outGene[['DLPFC_matchQSV']], top$gene)

outFeat[['HIPPO']]$exon <- cbind(outFeat[['HIPPO']]$exon, top$exon)
outFeat[['HIPPO']]$jxn <- cbind(outFeat[['HIPPO']]$jxn, top$jxn)
outFeat[['HIPPO']]$tx <- cbind(outFeat[['HIPPO']]$tx, top$tx)
outFeat[['DLPFC']]$exon <- cbind(outFeat[['DLPFC']]$exon, top$exon)
outFeat[['DLPFC']]$jxn <- cbind(outFeat[['DLPFC']]$jxn, top$jxn)
outFeat[['DLPFC']]$tx <- cbind(outFeat[['DLPFC']]$tx, top$tx)

message(paste(Sys.time(), 'saving DE info'))
save(outGene, outFeat, file = 'rda/out_info_withInteraction.Rdata')

comp <- function(out) {
    table('SCZD' = factor(out$adj.P.Val < 0.05, levels = c('TRUE', 'FALSE')), 'Interaction' = factor(out$interaction_adj.P.Val < 0.05, levels = c('TRUE', 'FALSE')))
}

comp_hippo <- list(
    'gene' = comp(outGene[['HIPPO_matchQSV']]),
    'exon' = comp(outFeat[['HIPPO']]$exon),
    'jxn' = comp(outFeat[['HIPPO']]$jxn),
    'tx' = comp(outFeat[['HIPPO']]$tx)
)
comp_hippo
# > comp_hippo
# $gene
#        Interaction
# SCZD     TRUE FALSE
#   TRUE      4    44
#   FALSE   107 24497
#
# $exon
#        Interaction
# SCZD      TRUE  FALSE
#   TRUE      25    172
#   FALSE    772 395610
#
# $jxn
#        Interaction
# SCZD      TRUE  FALSE
#   TRUE       4     37
#   FALSE    172 296968
#
# $tx
#        Interaction
# SCZD     TRUE FALSE
#   TRUE      0     0
#   FALSE     7 92725


comp_dlpfc <- list(
    'gene' = comp(outGene[['DLPFC_matchQSV']]),
    'exon' = comp(outFeat[['DLPFC']]$exon),
    'jxn' = comp(outFeat[['DLPFC']]$jxn),
    'tx' = comp(outFeat[['DLPFC']]$tx)
)
comp_dlpfc
# > comp_dlpfc
# $gene
#        Interaction
# SCZD     TRUE FALSE
#   TRUE      5   240
#   FALSE   106 24301
#
# $exon
#        Interaction
# SCZD      TRUE  FALSE
#   TRUE      12    428
#   FALSE    785 395354
#
# $jxn
#        Interaction
# SCZD      TRUE  FALSE
#   TRUE       6     31
#   FALSE    170 296974
#
# $tx
#        Interaction
# SCZD     TRUE FALSE
#   TRUE      0     6
#   FALSE     7 92719

lapply(comp_hippo, chisq.test)
# $gene
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 50.219, df = 1, p-value = 1.375e-12
#
#
# $exon
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 1471.2, df = 1, p-value < 2.2e-16
#
#
# $jxn
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 497.89, df = 1, p-value < 2.2e-16
#
#
# $tx
#
#     Pearson's Chi-squared test
#
# data:  X[[i]]
# X-squared = NaN, df = 1, p-value = NA

lapply(comp_dlpfc, chisq.test)
# $gene
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 10.612, df = 1, p-value = 0.001123
#
#
# $exon
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 127.84, df = 1, p-value < 2.2e-16
#
#
# $jxn
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 1370.5, df = 1, p-value < 2.2e-16
#
#
# $tx
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  X[[i]]
# X-squared = 1.2966e-28, df = 1, p-value = 1

## Interaction (Dx by region) DE genes
get_de <- function(x) {
    x$ensemblID[x$interaction_adj.P.Val < 0.05]
}
genes <- list(
    'gene' = get_de(outGene[[1]]),
    'exon' = get_de(outFeat[[1]]$exon),
    'jxn' = get_de(outFeat[[1]]$jxn),
    'tx' = get_de(outFeat[[1]]$tx)
)
stopifnot(identical(
    sapply(top, function(x) {
        sum(x$interaction_adj.P.Val < 0.05)
    }),
    sapply(genes, length)
))

## Remove NAs: un-annotated junctions
genes <- lapply(genes, function(x) { x[!is.na(x)] })
sapply(genes, length)
# gene exon  jxn   tx
#  111  797  147    7

## Pretty venn code
venn_cols <- brewer.pal('Set1', n = 4)
names(venn_cols) <- names(genes)
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
make_venn(genes)
make_venn(genes[c('gene', 'exon', 'jxn')])
dev.off()
system('rm VennDiagram*.log')


## Gene ontology analysis
uni <- c(outGene[[1]]$ensemblID, outFeat[[1]]$exon$ensemblID, outFeat[[1]]$jxn$ensemblID)
uni <- unique(uni[!is.na(uni)])
length(uni)
# [1] 27514

run_go <- function(genes, ont = c('BP', 'MF', 'CC')) {
    ## Change to ENSEMBL ids and remove NAs
    #genes_ens <- lapply(lapply(genes, function(x) { gsub('\\..*', '', x) }), function(y) y[!is.na(y)])
    ## They are already ENSEMBL ids in this case
    genes_ens <- genes

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

## GO for the DE genes
if(!file.exists('rda/go_de_genes.Rdata')) {
    system.time( go_de_genes <- run_go(genes) )
    message(paste(Sys.time(), 'saving rda/go_de_genes.Rdata'))
    save(go_de_genes, file = 'rda/go_de_genes.Rdata')
} else {
    message(paste(Sys.time(), 'loading rda/go_de_genes.Rdata'))
    load('rda/go_de_genes.Rdata', verbose = TRUE)
}
sapply(go_de_genes, class)
#                     BP                     MF                     CC
# "compareClusterResult" "compareClusterResult" "compareClusterResult"
#                   KEGG
# "compareClusterResult"

plot_go <- function(go_cluster, cat = 10) {
    lapply(names(go_cluster), function(bp) {
        go <- go_cluster[[bp]]
        if(is.null(go)) {
            message(paste(Sys.time(), 'found no results for', bp))
            return(NULL)
        }
        
        print(plot(go, title = paste('ontology:', bp), font.size = 18, showCategory = cat, includeAll = TRUE))
        return(NULL)
    })
}

## Made with conda_R/3.4.x
pdf('pdf/go_de_genes.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_go(go_de_genes)
dev.off()


pdf('pdf/go_all_de_genes.pdf', width = 16, height = 70, useDingbats = FALSE)
plot_go(go_de_genes, cat = NULL)
dev.off()


get_ylab <- function(type, cleaned = FALSE) {
    if(type %in% c('gene', 'exon', 'jxn')) {
        res <- 'log2(CPM + 0.5)'
    } else {
        res <- 'log2(TPM + 0.5)'
    }

    if(cleaned) res <- paste(res, '- covariate effects removed')
    return(res)
}

get_main <- function(i, region = 'HIPPO', type = 'gene') {
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

    # rse <- if(region == 'HIPPO') rse_gene_HIPPO else rse_gene_DLPFC

    topnow <- outGene[[paste0(region, '_matchQSV')]]

    j <- which(topnow$ensemblID == genes$gene[i])
    k <- which(names(rowRanges(rse)) == rownames(topnow)[j])

    paste(if(type != 'jxn') mcols(rowRanges(rse))[, var][k] else rownames(topnow)[j],
          if(is.na(mcols(rowRanges(rse))[, vars][k])) '' else mcols(rowRanges(rse))[, vars][k],
          paste0('FDR=',
          signif(topnow$interaction_adj.P.Val[j], 3), '\n', 
          'DLPFC: LFC=', signif(outGene[['DLPFC_matchQSV']]$logFC[j], 3), ', FDR=', signif(outGene[['DLPFC_matchQSV']]$adj.P.Val[j], 3), '; ',
          'HIPPO: LFC=', signif(outGene[['HIPPO_matchQSV']]$logFC[j], 3), ', FDR=', signif(outGene[['HIPPO_matchQSV']]$adj.P.Val[j], 3)
          ))
}

get_ylim_mult <- function(rang) {
    c(
        ifelse(sign(rang[1]) == 1, 0.95, 1.05),
        ifelse(sign(rang[2]) == 1, 1.05, 0.95)
    )
}

## Add an interaction by region for Dx
design <- cbind(with(colData(rse), model.matrix(~ Dx * Region)),
    modQsva[, -grep('Dx|RegionHIPPO|Intercept', colnames(modQsva))]
)
stopifnot(is.fullrank(design))
maincol <- grep(':', colnames(design))
mod <- design[, c(1, maincol, 2:(maincol - 1), (maincol + 1):ncol(design))]

exprsCleaned <- cleaningY(exprsNorm$gene, mod, P = 2)

plot_top <- function(i, normtype = 'cleaned') {
    g <- genes$gene[i]
    j <- which(outGene[['HIPPO_matchQSV']]$ensemblID == g)

    set.seed(20180917)
    dx <- factor(ifelse(colData(rse)$Dx == 'Control', 'Control', 'SCZD'), levels = c('Control', 'SCZD'))
    region <- colData(rse)$Region

    y <- if(normtype == 'expr') exprsNorm$gene[j, ] else exprsCleaned[j, ]
    boxplot(y ~ dx + region,
            ylab = get_ylab('gene', cleaned = normtype == 'cleaned'),
            main = get_main(i),
            col = rep(c('orchid1', 'aquamarine1'), 2),
            ylim = abs(range(y)) * get_ylim_mult(range(y)) * sign(range(y)),
            outline = FALSE)
    points(y ~ jitter(as.integer(dx) + ifelse(region == 'HIPPO', 2, 0), 1), pch = 21, #pch = ifelse(dx == 'Control', 21, 22),
           bg = ifelse(region == 'DLPFC', 'dark orange', 'skyblue3'))
}

## Order by FDR
genes$gene <- genes$gene[order(outGene[['HIPPO_matchQSV']]$interaction_adj.P.Val[match(genes$gene, outGene[['HIPPO_matchQSV']]$ensemblID)])]

## Make the plots
pdf('pdf/de_gene_expr.pdf')
for(i in seq_len(length(genes$gene))) {
    plot_top(i, 'expr')
}
dev.off()

pdf('pdf/de_gene_cleaned.pdf')
for(i in seq_len(length(genes$gene))) {
    plot_top(i, 'cleaned')
}
dev.off()

## Get the SCZD case-control DE genes among these 111 DE genes by Dx * Region
which(outGene[['HIPPO_matchQSV']]$adj.P.Val[match(genes$gene, outGene[['HIPPO_matchQSV']]$ensemblID)] < 0.05)
# [1] 57 63 83 95
which(outGene[['DLPFC_matchQSV']]$adj.P.Val[match(genes$gene, outGene[['DLPFC_matchQSV']]$ensemblID)] < 0.05)
# [1] 13 18 53 61 80

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# Session info ----------------------------------------------------------------------------------------------------------
#  setting  value
#  version  R version 3.5.0 Patched (2018-04-30 r74679)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  tz       US/Eastern
#  date     2018-09-17
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  AnnotationDbi          1.42.1    2018-05-17 Bioconductor
#  assertthat             0.2.0     2017-04-11 CRAN (R 3.5.0)
#  base                 * 3.5.0     2018-05-02 local
#  bindr                  0.1.1     2018-03-13 CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 CRAN (R 3.5.0)
#  Biobase              * 2.40.0    2018-05-02 Bioconductor
#  BiocGenerics         * 0.26.0    2018-05-03 Bioconductor
#  BiocParallel         * 1.14.2    2018-07-08 Bioconductor
#  bit                    1.1-14    2018-05-29 CRAN (R 3.5.0)
#  bit64                  0.9-7     2017-05-08 CRAN (R 3.5.0)
#  bitops                 1.0-6     2013-08-17 CRAN (R 3.5.0)
#  blob                   1.1.1     2018-03-25 CRAN (R 3.5.0)
#  caTools                1.17.1.1  2018-07-20 CRAN (R 3.5.0)
#  clusterProfiler      * 3.8.1     2018-06-13 Bioconductor
#  colorout             * 1.2-0     2018-05-02 Github (jalvesaq/colorout@c42088d)
#  colorspace             1.3-2     2016-12-14 CRAN (R 3.5.0)
#  compiler               3.5.0     2018-05-02 local
#  cowplot                0.9.3     2018-07-15 CRAN (R 3.5.0)
#  crayon                 1.3.4     2017-09-16 CRAN (R 3.5.0)
#  data.table             1.11.4    2018-05-27 cran (@1.11.4)
#  datasets             * 3.5.0     2018-05-02 local
#  DBI                    1.0.0     2018-05-02 CRAN (R 3.5.0)
#  DelayedArray         * 0.6.5     2018-08-15 Bioconductor
#  devtools             * 1.13.6    2018-06-27 CRAN (R 3.5.0)
#  digest                 0.6.15    2018-01-28 CRAN (R 3.5.0)
#  DO.db                  2.9       2018-05-03 Bioconductor
#  DOSE                   3.6.1     2018-06-20 Bioconductor
#  dplyr                  0.7.6     2018-06-29 CRAN (R 3.5.0)
#  enrichplot             1.0.2     2018-09-11 Bioconductor
#  fastmatch              1.1-0     2017-01-28 CRAN (R 3.5.0)
#  fgsea                  1.6.0     2018-05-03 Bioconductor
#  formatR                1.5       2017-04-25 CRAN (R 3.5.0)
#  futile.logger        * 1.4.3     2016-07-10 CRAN (R 3.5.0)
#  futile.options         1.0.1     2018-04-20 CRAN (R 3.5.0)
#  gdata                  2.18.0    2017-06-06 CRAN (R 3.5.0)
#  GenomeInfoDb         * 1.16.0    2018-05-03 Bioconductor
#  GenomeInfoDbData       1.1.0     2018-04-17 Bioconductor
#  GenomicRanges        * 1.32.6    2018-07-20 Bioconductor
#  ggforce                0.1.3     2018-07-07 CRAN (R 3.5.0)
#  ggplot2              * 3.0.0     2018-07-03 CRAN (R 3.5.0)
#  ggraph                 1.0.2     2018-07-07 CRAN (R 3.5.0)
#  ggrepel                0.8.0     2018-05-09 CRAN (R 3.5.0)
#  ggridges               0.5.0     2018-04-05 CRAN (R 3.5.0)
#  glue                   1.3.0     2018-07-17 CRAN (R 3.5.0)
#  GO.db                  3.6.0     2018-05-03 Bioconductor
#  GOSemSim               2.6.0     2018-05-03 Bioconductor
#  gplots               * 3.0.1     2016-03-30 CRAN (R 3.5.0)
#  graphics             * 3.5.0     2018-05-02 local
#  grDevices            * 3.5.0     2018-05-02 local
#  grid                 * 3.5.0     2018-05-02 local
#  gridExtra              2.3       2017-09-09 CRAN (R 3.5.0)
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.5.0)
#  gtools                 3.8.1     2018-06-26 CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 CRAN (R 3.5.0)
#  htmlwidgets            1.2       2018-04-19 CRAN (R 3.5.0)
#  httpuv                 1.4.5     2018-07-19 CRAN (R 3.5.0)
#  igraph                 1.2.2     2018-07-27 CRAN (R 3.5.0)
#  IRanges              * 2.14.10   2018-05-17 Bioconductor
#  jaffelab             * 0.99.21   2018-05-03 Github (LieberInstitute/jaffelab@7ed0ab7)
#  KernSmooth             2.23-15   2015-06-29 CRAN (R 3.5.0)
#  lambda.r               1.2.3     2018-05-17 CRAN (R 3.5.0)
#  later                  0.7.4     2018-08-31 CRAN (R 3.5.0)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.5.0)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.5.0)
#  limma                * 3.36.2    2018-06-21 Bioconductor
#  magrittr               1.5       2014-11-22 CRAN (R 3.5.0)
#  MASS                   7.3-50    2018-04-30 CRAN (R 3.5.0)
#  Matrix                 1.2-14    2018-04-13 CRAN (R 3.5.0)
#  matrixStats          * 0.54.0    2018-07-23 CRAN (R 3.5.0)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.5.0)
#  methods              * 3.5.0     2018-05-02 local
#  mime                   0.5       2016-07-07 CRAN (R 3.5.0)
#  munsell                0.5.0     2018-06-12 CRAN (R 3.5.0)
#  parallel             * 3.5.0     2018-05-02 local
#  pillar                 1.3.0     2018-07-14 CRAN (R 3.5.0)
#  pkgconfig              2.0.1     2017-03-21 CRAN (R 3.5.0)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 CRAN (R 3.5.0)
#  purrr                  0.2.5     2018-05-29 CRAN (R 3.5.0)
#  qvalue                 2.12.0    2018-05-03 Bioconductor
#  R6                     2.2.2     2017-06-17 CRAN (R 3.5.0)
#  rafalib              * 1.0.0     2015-08-09 CRAN (R 3.5.0)
#  RColorBrewer         * 1.1-2     2014-12-07 CRAN (R 3.5.0)
#  Rcpp                   0.12.18   2018-07-23 CRAN (R 3.5.0)
#  RCurl                  1.95-4.11 2018-07-15 CRAN (R 3.5.0)
#  reshape2               1.4.3     2017-12-11 CRAN (R 3.5.0)
#  rlang                  0.2.1     2018-05-30 cran (@0.2.1)
#  rmote                * 0.3.4     2018-05-02 deltarho (R 3.5.0)
#  RSQLite                2.1.1     2018-05-06 CRAN (R 3.5.0)
#  rvcheck                0.1.0     2018-05-23 CRAN (R 3.5.0)
#  S4Vectors            * 0.18.3    2018-06-13 Bioconductor
#  scales                 1.0.0     2018-08-09 CRAN (R 3.5.0)
#  segmented              0.5-3.0   2017-11-30 CRAN (R 3.5.0)
#  servr                  0.10      2018-05-30 CRAN (R 3.5.0)
#  splines                3.5.0     2018-05-02 local
#  stats                * 3.5.0     2018-05-02 local
#  stats4               * 3.5.0     2018-05-02 local
#  stringi                1.2.4     2018-07-20 CRAN (R 3.5.0)
#  stringr                1.3.1     2018-05-10 CRAN (R 3.5.0)
#  SummarizedExperiment * 1.10.1    2018-05-17 Bioconductor
#  tibble                 1.4.2     2018-01-22 CRAN (R 3.5.0)
#  tidyr                  0.8.1     2018-05-18 CRAN (R 3.5.0)
#  tidyselect             0.2.4     2018-02-26 CRAN (R 3.5.0)
#  tools                  3.5.0     2018-05-02 local
#  tweenr                 0.1.5     2016-10-10 CRAN (R 3.5.0)
#  units                  0.6-0     2018-06-09 CRAN (R 3.5.0)
#  UpSetR                 1.3.3     2017-03-21 CRAN (R 3.5.0)
#  utils                * 3.5.0     2018-05-02 local
#  VennDiagram          * 1.6.20    2018-03-28 CRAN (R 3.5.0)
#  viridis                0.5.1     2018-03-29 CRAN (R 3.5.0)
#  viridisLite            0.3.0     2018-02-01 CRAN (R 3.5.0)
#  withr                  2.1.2     2018-03-15 CRAN (R 3.5.0)
#  xfun                   0.3       2018-07-06 CRAN (R 3.5.0)
#  XVector                0.20.0    2018-05-03 Bioconductor
#  zlibbioc               1.26.0    2018-05-02 Bioconductor
#
