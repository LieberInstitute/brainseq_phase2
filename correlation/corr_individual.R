## Usage:
# qrsh -l mem_free=100G,h_vmem=100G,h_fsize=100G
# Rscript corr_individual.R  > logs/corr_individual.txt 2>&1

library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('scales')
library('derfinder')
library('clusterProfiler')

dir.create('pdf', showWarnings = FALSE)
dir.create('rda', showWarnings = FALSE)

## Load expression data (cleaned and uncleaned) + rse objects for the pheno data
load('rda/expr_and_cleaned.Rdata', verbose = TRUE)
load('rda/rse_and_modQsva.Rdata', verbose = TRUE)

## Function for computing paired cors
paircor <- function(x, y) {
    stopifnot(nrow(x) == nrow(y))
    stopifnot(ncol(x) == ncol(y))
    sapply(seq_len(ncol(x)), function(i) {
        cor(x[, i], y[, i])
    })
}

## Test paircor
test <- cor(cleaned[['DLPFC']][['geneRpkm']][1:100, ], cleaned[['HIPPO']][['geneRpkm']][1:100, ])
test2 <- paircor(cleaned[['DLPFC']][['geneRpkm']][1:100, ], cleaned[['HIPPO']][['geneRpkm']][1:100, ])
stopifnot(all(diag(test) - test2 == 0))

computecor <- function(exp) {
    sets <- names(exp[[1]])
    res <- lapply(sets, function(feature) {
        message(paste(Sys.time(), 'processing feature', feature))
        paircor(exp[['DLPFC']][[feature]], exp[['HIPPO']][[feature]])
    })
    names(res) <- sets
    return(res)
}

## Compute individual level correlations across features
indv_expr <- computecor(expr)
indv_cleaned <- computecor(cleaned)
save(indv_expr, indv_cleaned, file = 'rda/indv_corr.Rdata')


pdf('pdf/indv_hist.pdf', useDingbats = FALSE)
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = 'light blue')
}, indv_expr, names(indv_expr), 'expr')
mapply(function(dat, set, type) {
    hist(dat, main = paste(type, '-', set), col = "#009E73")
}, indv_cleaned, names(indv_cleaned), 'cleaned expr (keeping Dx)')
dev.off()

pdf('pdf/indv_expr_vs_cleaned.pdf', useDingbats = FALSE)
mapply(function(x, y, set) {
    m <- !is.na(x) & !is.na(y)
    plot(x[m], y[m], xlab = 'expr', ylab = 'cleaned', main = set, pch = 19, col = alpha("#D55E00", 1))
    abline(a = 0, b = 1, col = 'blue')
}, indv_expr, indv_cleaned, names(indv_expr))
dev.off()

pdf('pdf/indv_expr_vs_cleaned_sameLims.pdf', useDingbats = FALSE)
mapply(function(x, y, set) {
    m <- !is.na(x) & !is.na(y)
    lims <- c(min(c(x[m], y[m])), max(c(x[m], y[m])))
    plot(x[m], y[m], xlab = 'expr', ylab = 'cleaned', main = set, pch = 19, col = alpha("#D55E00", 1),
        xlim = lims, ylim = lims
    )
    abline(a = 0, b = 1, col = 'blue')
}, indv_expr, indv_cleaned, names(indv_expr))
dev.off()


pdf('pdf/indv_box_sczd.pdf', useDingbats = FALSE)
mapply(function(x, y, set, type) {
    m <- !is.na(x)
    dx <- colData(y)$Dx[m]
    dx <- factor(ifelse(dx == 'Schizo', 'SCZD', ifelse(dx == 'Control', 'Control', 'hmm')))
    f <- lm(x[m] ~ dx)
    p <- summary(f)$coef[2, 4]
    ylim <- range(x[m])
    boxplot(x[m] ~ dx, main = paste(type, '-', set, '\n p-value:', signif(p, 3)),
        xlab = 'SCZD diagnosis', outline = FALSE, yab = 'Correlation', ylim = ylim)
    points(x[m] ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = 'light blue')
}, indv_expr, simple_rse[[1]], names(indv_expr), 'expr')
mapply(function(x, y, set, type) {
    m <- !is.na(x)
    dx <- colData(y)$Dx[m]
    f <- lm(x[m] ~ dx)
    p <- summary(f)$coef[2, 4]
    ylim <- range(x[m])
    dx <- factor(ifelse(dx == 'Schizo', 'SCZD', ifelse(dx == 'Control', 'Control', 'hmm')))
    boxplot(x[m] ~ dx, main = paste(type, '-', set, '\n p-value:', signif(p, 3)),
        xlab = 'SCZD diagnosis', outline = FALSE, ylab = 'Correlation', ylim = ylim)
    points(x[m] ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = '#009E73')
}, indv_cleaned, simple_rse[[1]], names(indv_cleaned), 'cleaned expr (keeping Dx)')
dev.off()



## Explore covariates to see if low correlation is related to any main covariate
pdf('pdf/indv_corr_vs_covariates.pdf', useDingbats = FALSE)
lapply(c('Age'), function(var) {
    message(paste(Sys.time(), 'processing', var))
        pd <- colData(simple_rse[[1]][[1]])
        y <- pd[, var]
        if(is(y, 'CompressedNumericList') | is(y, 'CompressedIntegerList')) {
            y <- mean(y)
        } else if (is.character(y)) {
            y <- as.factor(y)
        }
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - expr', pch = 19, col = 'deepskyblue3')
        }, indv_expr, names(indv_expr))
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - cleaned', pch = 19, col = '#009E73')
        }, indv_cleaned, names(indv_cleaned))
        return(NULL)
})
lapply(c('Sex', 'Race'), function(var) {
    message(paste(Sys.time(), 'processing', var))
        pd <- colData(simple_rse[[1]][[1]])
        y <- pd[, var]
        if(is(y, 'CompressedNumericList') | is(y, 'CompressedIntegerList')) {
            y <- mean(y)
        } else if (is.character(y)) {
            y <- as.factor(y)
        }
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - expr')
        }, indv_expr, names(indv_expr))
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = var, xlab = 'Corr - cleaned')
        }, indv_cleaned, names(indv_cleaned))
        return(NULL)
})

lapply(c('snpPC1', 'snpPC2', 'RIN', 'rRNA_rate', 'totalAssignedGene', 'mitoRate', 'totalMapped', 'mitoMapped', 'overallMapRate'), function(var) {
    message(paste(Sys.time(), 'processing', var))
    lapply(names(simple_rse), function(region) {
        pd <- colData(simple_rse[[region]][[1]])
        y <- pd[, var]
        if(is(y, 'CompressedNumericList') | is(y, 'CompressedIntegerList')) {
            y <- mean(y)
        } else if (is.character(y)) {
            y <- as.factor(y)
        }
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = paste0(region, ': ', var), xlab = 'Corr - expr', pch = 19, col = 'deepskyblue3')
        }, indv_expr, names(indv_expr))
        mapply(function(x, set) {
            plot(y ~ x, main = set, ylab = paste0(region, ': ', var), xlab = 'Corr - cleaned', pch = 19, col = '#009E73')
        }, indv_cleaned, names(indv_cleaned))
        return(NULL)
    })
})
dev.off()



## Explore expression across all genes for a few samples with either very low or very high correlation (expr, not cleaned)
pdf('pdf/indv_expr_high_low_corr.pdf', useDingbats = FALSE)
lapply(names(indv_expr)[1], function(set) {
    i <- match(1:5, rank(indv_expr[[set]]))
    j <- match(261:265, rank(indv_expr[[set]]))
    
    mapply(function(k, rk) {
        plot(expr[[1]][[set]][, k], expr[[2]][[set]][, k], xlab = 'DLPFC', ylab = 'HIPPO', main = paste(set, 'rank:', rk, '- expr'), pch = 19, col = alpha('deepskyblue3', 1/10))
        abline(a = 0, b = 1, col = 'grey80')
        
        plot(cleaned[[1]][[set]][, k], cleaned[[2]][[set]][, k], xlab = 'DLPFC', ylab = 'HIPPO', main = paste(set, 'rank:', rk, '- cleaned expr'), pch = 19, col = alpha('#009E73', 1/10))
        abline(a = 0, b = 1, col = 'grey80')
        
    }, c(i, j), c(1:5, 261:265))
    return(NULL)
})
dev.off()


## Get KEGG sets
## code by Andrew Jaffe
KEGG_DATA <- clusterProfiler:::get_data_from_KEGG_db("hsa")
keggList = get("PATHID2EXTID", envir = KEGG_DATA)
keggTerms <- DOSE:::TERM2NAME(names(keggList), KEGG_DATA)

## Split genes by KEGG pathway
kegg_df <- data.frame(
    'entrez' = unlist(keggList),
    'kegg_path' = rep(names(keggList), sapply(keggList, length)),
    'kegg_term' = rep(keggTerms, sapply(keggList, length)),
    stringsAsFactors = FALSE
)

m <- match(rowRanges(simple_rse[['DLPFC']][['gene']])$EntrezID, as.integer(kegg_df$entrez))
m1 <- !is.na(m)
table(m1)
m2 <- m[m1]

## 1,420 KEGG entrez ids are missing
length(unique(kegg_df$entrez))
m_rev <- match(as.integer(unique(kegg_df$entrez)), rowRanges(simple_rse[['DLPFC']][['gene']])$EntrezID)
m1_rev <- !is.na(m_rev)
table(m1_rev)

kegg_df$m_gene <- match(as.integer(kegg_df$entrez), rowRanges(simple_rse[['DLPFC']][['gene']])$EntrezID)
kegg_df_comp <- kegg_df[!is.na(kegg_df$m_gene), ]

kegg_genes <- split(kegg_df_comp$m_gene, kegg_df_comp$kegg_path)
length(kegg_genes)
kegg_genes <- kegg_genes[sapply(kegg_genes, length) > 1]
length(kegg_genes)


computecor_kegg <- function(exp) {
    pathways <- names(kegg_genes)
    res <- lapply(pathways, function(pathway) {
        message(paste(Sys.time(), 'processing pathway', pathway))
        paircor(exp[['DLPFC']][['geneRpkm']][kegg_genes[[pathway]], ], exp[['HIPPO']][['geneRpkm']][kegg_genes[[pathway]], ])
    })
    names(res) <- pathways
    return(res)
}

## Compute individual level correlations across features
indv_kegg_expr <- computecor_kegg(expr)
indv_kegg_cleaned <- computecor_kegg(cleaned)
save(indv_kegg_expr, indv_kegg_cleaned, kegg_df, kegg_df_comp, kegg_genes, file = 'rda/indv_kegg_corr.Rdata')


pdf('pdf/indv_mean_kegg.pdf', useDingbats = FALSE)
plot(sapply(indv_kegg_expr, mean), sapply(indv_kegg_cleaned, mean), xlab = 'expr', ylab = 'cleaned', main = paste('mean correlation for', length(kegg_genes), 'KEGG pathways'), pch = 19, col = 'darkviolet')
hist(sapply(indv_kegg_expr, mean), main = 'mean KEGG corr - expr', col = 'deepskyblue3')
hist(sapply(indv_kegg_cleaned, mean), main = 'mean KEGG corr - cleaned', col = '#009E73')
dev.off()

top_kegg <- function(kegg_cor, n = 10) {
    top <- head(sort(sapply(kegg_cor, mean), decreasing = TRUE), n = n)
    res <- cbind(kegg_df_comp[match(names(top), kegg_df_comp$kegg_path), ], mean_corr = top)
    res$n_genes <- sapply(res$kegg_path, function(x) { sum(kegg_df_comp$kegg_path == x) })
    rownames(res) <- NULL
    return(res[, c('kegg_path', 'kegg_term', 'mean_corr', 'n_genes')])
}
options(width = 120)
print('expr')
top_kegg(indv_kegg_expr)
print('cleaned')
top_kegg(indv_kegg_cleaned)

# > print('expr')
# [1] "expr"
# > top_kegg(indv_kegg_expr)
#    kegg_path                                           kegg_term mean_corr n_genes
# 1   hsa00400 Phenylalanine, tyrosine and tryptophan biosynthesis 0.9949677       3
# 2   hsa03450                          Non-homologous end-joining 0.9868043      11
# 3   hsa00785                              Lipoic acid metabolism 0.9846689       3
# 4   hsa04744                                   Phototransduction 0.9827878      12
# 5   hsa00740                               Riboflavin metabolism 0.9815958       7
# 6   hsa00524                 Butirosin and neomycin biosynthesis 0.9808939       4
# 7   hsa00471              D-Glutamine and D-glutamate metabolism 0.9808905       4
# 8   hsa04114                                      Oocyte meiosis 0.9800333      95
# 9   hsa03010                                            Ribosome 0.9795631      86
# 10  hsa01040             Biosynthesis of unsaturated fatty acids 0.9795592      21
# > print('cleaned')
# [1] "cleaned"
# > top_kegg(indv_kegg_cleaned)
#    kegg_path                                           kegg_term mean_corr n_genes
# 1   hsa00780                                   Biotin metabolism 1.0000000       2
# 2   hsa00785                              Lipoic acid metabolism 0.9813169       3
# 3   hsa00400 Phenylalanine, tyrosine and tryptophan biosynthesis 0.9783005       3
# 4   hsa00524                 Butirosin and neomycin biosynthesis 0.9763295       4
# 5   hsa00740                               Riboflavin metabolism 0.9720471       7
# 6   hsa04614                            Renin-angiotensin system 0.9648249       9
# 7   hsa03450                          Non-homologous end-joining 0.9602165      11
# 8   hsa00860                Porphyrin and chlorophyll metabolism 0.9554085      23
# 9   hsa00190                           Oxidative phosphorylation 0.9496137     122
# 10  hsa00062                               Fatty acid elongation 0.9489531       8


kegg_sczd <- function(kegg_cor, type) {
    dx <- colData(simple_rse[['DLPFC']][['gene']])$Dx
    dx <- factor(ifelse(dx == 'Schizo', 'SCZD', ifelse(dx == 'Control', 'Control', 'hmm')))
    # > levels(dx)
    # [1] "Control" "SCZD"
    
    res <- do.call(rbind, mapply(function(pathway, pathname) {
        f <- lm(pathway ~ dx)
        p <- summary(f)$coef[2, 4]
        ylim <- range(pathway)
        
        boxplot(pathway ~ dx, main = paste0(pathname, ': ',
            kegg_df_comp$kegg_term[which(kegg_df_comp$kegg_path == pathname)[1]],
            '\np-value: ', signif(p, 3)),
            ylim = ylim,
            xlab = 'SCZD diagnosis', outline = FALSE, ylab = paste('Correlation -', type))
        points(pathway ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = ifelse(type == 'expr', 'deepskyblue3', '#009E73'))
                
        c(summary(f)$coef[2, ], mean_sczd = mean(pathway[dx == 'SCZD']), mean_control = mean(pathway[dx == 'Control']))
    }, kegg_cor, names(kegg_cor), SIMPLIFY = FALSE))
    m <- match(rownames(res), kegg_df_comp$kegg_path)
    res <- cbind(res, kegg_df_comp[m, c('kegg_path', 'kegg_term')])
    res$n_genes <- sapply(res$kegg_path, function(x) { sum(kegg_df_comp$kegg_path == x) })
    res$padj <- p.adjust(res[, 'Pr(>|t|)'], method = 'fdr')
    res <- res[order(res$padj, decreasing = FALSE), ]
    rownames(res) <- NULL
    
    return(res)
}

pdf('pdf/indv_box_sczd_kegg.pdf', useDingbats = FALSE)
indv_kegg_expr_sczd <- kegg_sczd(indv_kegg_expr, 'expr')
indv_kegg_cleaned_sczd <- kegg_sczd(indv_kegg_cleaned, 'cleaned')
dev.off()

nrow(indv_kegg_expr_sczd)
# [1] 227
table(indv_kegg_expr_sczd$padj < 0.05)
# FALSE
#   227
table(indv_kegg_cleaned_sczd$padj < 0.05)
# FALSE  TRUE
#   214    13
options(width = 200)
subset(indv_kegg_cleaned_sczd, padj < 0.05)
#        Estimate   Std. Error   t value     Pr(>|t|) mean_sczd mean_control kegg_path                               kegg_term n_genes         padj
# 1  -0.011095422 0.0023079374 -4.807506 2.573804e-06 0.8441098    0.8552052  hsa03320                  PPAR signaling pathway      51 0.0005842534
# 2  -0.004573640 0.0009932781 -4.604591 6.440982e-06 0.8916175    0.8961911  hsa04142                                Lysosome     109 0.0007310515
# 3  -0.014778589 0.0037914968 -3.897825 1.232383e-04 0.7534047    0.7681833  hsa04640              Hematopoietic cell lineage      42 0.0093250285
# 4  -0.009263283 0.0024614478 -3.763347 2.067723e-04 0.8204990    0.8297623  hsa05145                           Toxoplasmosis     105 0.0117343254
# 5  -0.004902039 0.0013356126 -3.670255 2.933653e-04 0.8842804    0.8891824  hsa04210                               Apoptosis      76 0.0133187844
# 6  -0.011181312 0.0032033927 -3.490459 5.652704e-04 0.7469369    0.7581182  hsa00511                Other glycan degradation      16 0.0169275367
# 7   0.009304367 0.0026565159  3.502470 5.414783e-04 0.8257343    0.8164300  hsa04740                  Olfactory transduction      41 0.0169275367
# 8  -0.015689972 0.0045146292 -3.475362 5.965652e-04 0.7248851    0.7405750  hsa05140                           Leishmaniasis      58 0.0169275367
# 9  -0.010150273 0.0029939956 -3.390210 8.056484e-04 0.7362409    0.7463911  hsa00071                  Fatty acid degradation      33 0.0203202440
# 10  0.009464357 0.0028201798  3.355941 9.076858e-04 0.9244389    0.9149745  hsa04744                       Phototransduction      12 0.0206044679
# 11 -0.002924630 0.0008956438 -3.265394 1.238149e-03 0.9067666    0.9096912  hsa00970             Aminoacyl-tRNA biosynthesis      41 0.0239204922
# 12 -0.004382762 0.0013447454 -3.259176 1.264519e-03 0.9349081    0.9392909  hsa01040 Biosynthesis of unsaturated fatty acids      21 0.0239204922
# 13 -0.015603131 0.0049429940 -3.156615 1.781964e-03 0.7893232    0.8049264  hsa00460              Cyanoamino acid metabolism       5 0.0311158305
save(indv_kegg_expr_sczd, indv_kegg_cleaned_sczd, file = 'rda/indv_kegg_sczd.Rdata')



## Now by GO groups
entrez <- rowRanges(simple_rse[['DLPFC']][['gene']])$EntrezID
length(entrez[!is.na(entrez)])
# [1] 16677
ggo_cc <- groupGO(as.character(entrez[!is.na(entrez)]), 'org.Hs.eg.db', ont = 'CC', level = 3)
ggo_bp <- groupGO(as.character(entrez[!is.na(entrez)]), 'org.Hs.eg.db', ont = 'BP', level = 3)
ggo_mf <- groupGO(as.character(entrez[!is.na(entrez)]), 'org.Hs.eg.db', ont = 'MF', level = 3)
all_ggo <- list('CC' = ggo_cc, 'BP' = ggo_bp, 'MF' = ggo_mf)
save(all_ggo, file = 'rda/all_ggo.Rdata')
(ggo_len <- sapply(all_ggo, function(x) { y <- summary(x); sum(y$Count > 1)}))
#  CC  BP  MF
# 327 338  95

computecor_go <- function(ggo, exp) {
    df <- ggo@result
    df <- df[df$Count > 1, ]
    gos <- as.character(df$ID)
    res <- lapply(gos, function(go) {
        message(paste(Sys.time(), 'processing pathway', go))
        genes <- strsplit(as.character(df$geneID[which(gos == go)]), '/')[[1]]
        m <- match(genes, as.character(rowRanges(simple_rse[['DLPFC']][['gene']])$EntrezID))
        stopifnot(!any(is.na(m)))
        paircor(exp[['DLPFC']][['geneRpkm']][m, ], exp[['HIPPO']][['geneRpkm']][m, ])
    })
    names(res) <- gos
    
    filt <- sapply(res, function(x) { !any(is.na(x))})
    print('Removing those with some NAs (FALSE in the table)')
    print(table(filt))
    res <- res[filt]
    
    return(res)
}

indv_go_expr <- lapply(all_ggo, computecor_go, exp = expr)
indv_go_cleaned <- lapply(all_ggo, computecor_go, exp = cleaned)
save(indv_go_expr, indv_go_cleaned, file = 'rda/indv_go_corr.Rdata')

pdf('pdf/indv_mean_go.pdf', useDingbats = FALSE)
for(i in 1:3) {
    plot(sapply(indv_go_expr[[i]], mean), sapply(indv_go_cleaned[[i]], mean), xlab = 'expr', ylab = 'cleaned', main = paste('mean correlation for', ggo_len[i], names(indv_go_expr)[i], 'GOs'), pch = 19, col = 'darkviolet')
    hist(sapply(indv_go_expr[[i]], mean), main = paste0('mean GO corr - expr -' , names(indv_go_expr)[i]), col = 'deepskyblue3')
    hist(sapply(indv_go_cleaned[[i]], mean), main = paste0('mean GO corr - cleaned - ', names(indv_go_expr)[i]), col = '#009E73')
}
dev.off()


kegg_go <- function(go_cor, go_info, type, go_type) {
    dx <- colData(simple_rse[['DLPFC']][['gene']])$Dx
    dx <- factor(ifelse(dx == 'Schizo', 'SCZD', ifelse(dx == 'Control', 'Control', 'hmm')))
    # > levels(dx)
    # [1] "Control" "SCZD"
    
    go_info <- go_info@result   
    go_info$ID <- as.character(go_info$ID)
    go_info$Description <- as.character(go_info$Description)
    go_info$GeneRatio <- as.character(go_info$GeneRatio)
    go_info$geneID <- as.character(go_info$geneID)
    
    
    
    res <- do.call(rbind, mapply(function(pathway, pathname) {
        f <- lm(pathway ~ dx)
        p <- summary(f)$coef[2, 4]
        ylim <- range(pathway)
        
        boxplot(pathway ~ dx, main = paste0(go_type, ' ', pathname, ': ',
            go_info$Description[which(go_info$ID == pathname)[1]],
            '\np-value: ', signif(p, 3)),
            ylim = ylim,
            xlab = 'SCZD diagnosis', outline = FALSE, ylab = paste('Correlation -', type))
        points(pathway ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = ifelse(type == 'expr', 'deepskyblue3', '#009E73'))
                
        c(summary(f)$coef[2, ], mean_sczd = mean(pathway[dx == 'SCZD']), mean_control = mean(pathway[dx == 'Control']))
    }, go_cor, names(go_cor), SIMPLIFY = FALSE))
    
    m <- match(rownames(res), go_info$ID)
    res <- cbind(res, go_info[m, ])
    res$n_genes <- as.integer(sapply(sapply(lapply(res$GeneRatio, strsplit, '/'), '[', 1), '[', 1))
    res$padj <- p.adjust(res[, 'Pr(>|t|)'], method = 'fdr')
    res <- res[order(res$padj, decreasing = FALSE), ]
    rownames(res) <- NULL
    
    return(res)
}

indv_go_expr_sczd <- indv_go_cleaned_sczd <- vector('list', 3)
names(indv_go_expr_sczd) <- names(indv_go_cleaned_sczd) <- names(indv_go_cleaned)
pdf('pdf/indv_box_sczd_go.pdf', useDingbats = FALSE)
for(i in 1:3) {
    message(paste(Sys.time(), 'processing', names(indv_go_cleaned)[i], '- expr'))
    indv_go_expr_sczd[[i]] <- kegg_go(indv_go_expr[[i]], all_ggo[[i]], 'expr', names(indv_go_expr)[i])
    message(paste(Sys.time(), 'processing', names(indv_go_cleaned)[i], '- cleaned'))
    indv_go_cleaned_sczd[[i]] <- kegg_go(indv_go_cleaned[[i]], all_ggo[[i]], 'cleaned', names(indv_go_cleaned)[i])
}
dev.off()


sapply(indv_go_expr_sczd, nrow)
#sapply(indv_go_cleaned_sczd, nrow)
#  CC  BP  MF
# 327 337  95

sapply(indv_go_expr_sczd, function(x) { sum(x$padj < 0.05 )})
# CC BP MF
#  5  0  2
sapply(indv_go_cleaned_sczd, function(x) { sum(x$padj < 0.05 )})
# CC BP MF
# 10 14  3
options(width = 200)
lapply(lapply(indv_go_expr_sczd, subset, padj < 0.05), function(x) { x[, -11]})
# $CC
#       Estimate   Std. Error   t value     Pr(>|t|) mean_sczd mean_control         ID                                        Description Count GeneRatio n_genes       padj
# 1 -0.002933467 0.0007621624 -3.848874 0.0001490311 0.9899106    0.9928441 GO:0002199                    zona pellucida receptor complex     8   8/16644       8 0.04837948
# 2  0.280051753 0.0820167908  3.414566 0.0007397474 0.9082569    0.6282051 GO:0005602                    complement component C1 complex     2   2/16644       2 0.04837948
# 3  0.040855594 0.0116233227  3.514967 0.0005177219 0.9262735    0.8854179 GO:0036454                              growth factor complex     4   4/16644       4 0.04837948
# 4  0.040855594 0.0116233227  3.514967 0.0005177219 0.9262735    0.8854179 GO:0016942 insulin-like growth factor binding protein complex     4   4/16644       4 0.04837948
# 5 -0.159255207 0.0465609938 -3.420357 0.0007248384 0.5741280    0.7333832 GO:0030906                  retromer, cargo-selective complex     4   4/16644       4 0.04837948
#
# $BP
#  [1] Estimate     Std. Error   t value      Pr(>|t|)     mean_sczd    mean_control ID           Description  Count        GeneRatio    n_genes      padj
# <0 rows> (or 0-length row.names)
#
# $MF
#     Estimate  Std. Error   t value     Pr(>|t|)  mean_sczd mean_control         ID             Description Count GeneRatio n_genes       padj
# 1  0.4249588 0.119323251  3.561408 0.0004377293 0.02752294   -0.3974359 GO:0050436     microfibril binding     2   2/16644       2 0.04158428
# 2 -0.0103953 0.003117011 -3.335022 0.0009757660 0.97852616    0.9889215 GO:0008384 IkappaB kinase activity     3   3/16644       3 0.04634888
lapply(lapply(indv_go_cleaned_sczd, subset, padj < 0.05), function(x) { x[, -11]})
# $CC
#        Estimate   Std. Error   t value     Pr(>|t|) mean_sczd mean_control         ID                              Description Count GeneRatio n_genes       padj
# 1  -0.002589671 0.0006587388 -3.931256 0.0001081201 0.9422399    0.9448296 GO:0005875           microtubule associated complex   135 135/16644     135 0.02870965
# 2  -0.006130560 0.0016981598 -3.610120 0.0003663869 0.9593958    0.9655264 GO:1990131                 Gtr1-Gtr2 GTPase complex     4   4/16644       4 0.02870965
# 3  -0.036624693 0.0101221684 -3.618266 0.0003555804 0.8154615    0.8520862 GO:0044326                     dendritic spine neck     4   4/16644       4 0.02870965
# 4  -0.011914673 0.0033403774 -3.566864 0.0004291340 0.8128534    0.8247681 GO:0045178                       basal part of cell    42  42/16644      42 0.02870965
# 5  -0.005732611 0.0016100037 -3.560619 0.0004389855 0.8879753    0.8937079 GO:0090543                            Flemming body    23  23/16644      23 0.02870965
# 6  -0.018027721 0.0053059401 -3.397649 0.0007849614 0.7793642    0.7973919 GO:0099240 intrinsic component of synaptic membrane     9   9/16644       9 0.04278039
# 7  -0.002330112 0.0007092071 -3.285517 0.0011562681 0.9672466    0.9695768 GO:0000974                            Prp19 complex    13  13/16644      13 0.04377292
# 8   0.056213484 0.0173373911  3.242327 0.0013386213 0.8003857    0.7441722 GO:0030689                              Noc complex     3   3/16644       3 0.04377292
# 9  -0.323453305 0.0996645652 -3.245419 0.0013247249 0.3944954    0.7179487 GO:0043511                          inhibin complex     2   2/16644       2 0.04377292
# 10 -0.007370197 0.0022702302 -3.246454 0.0013201080 0.8291516    0.8365218 GO:0072562                      blood microparticle    62  62/16644      62 0.04377292
#
# $BP
#        Estimate  Std. Error   t value     Pr(>|t|) mean_sczd mean_control         ID                                                      Description Count GeneRatio n_genes         padj
# 1  -0.008584444 0.001791930 -4.790613 2.781277e-06 0.9514168    0.9600012 GO:0044110                         growth involved in symbiotic interaction     6   6/16644       6 0.0009372903
# 2  -0.006666111 0.001477834 -4.510731 9.743652e-06 0.8921225    0.8987886 GO:0007585                                     respiratory gaseous exchange    50  50/16644      50 0.0016418054
# 3   0.057537157 0.014098412  4.081109 5.949020e-05 0.5501997    0.4926625 GO:0061744                                                   motor behavior     5   5/16644       5 0.0066827319
# 4  -0.008929998 0.002251035 -3.967063 9.388695e-05 0.9396881    0.9486181 GO:0008340                                  determination of adult lifespan    10  10/16644      10 0.0079099754
# 5  -0.005343728 0.001386367 -3.854482 1.458351e-04 0.8372702    0.8426139 GO:0001816                                              cytokine production   514 514/16644     514 0.0081722943
# 6  -0.007919468 0.002070762 -3.824423 1.637541e-04 0.8214242    0.8293436 GO:0044706                             multi-multicellular organism process   148 148/16644     148 0.0081722943
# 7  -0.005717320 0.001498620 -3.815056 1.697509e-04 0.9614569    0.9671742 GO:1900454             positive regulation of long term synaptic depression     4   4/16644       4 0.0081722943
# 8  -0.004628983 0.001272581 -3.637477 3.312693e-04 0.8771629    0.8817919 GO:0006457                                                  protein folding   224 224/16644     224 0.0139547173
# 9  -0.005932432 0.001661894 -3.569682 4.247576e-04 0.8449342    0.8508666 GO:0002440              production of molecular mediator of immune response   122 122/16644     122 0.0159048111
# 10  0.012302585 0.003605407  3.412260 7.457652e-04 0.8050615    0.7927589 GO:0045494                                   photoreceptor cell maintenance    25  25/16644      25 0.0251322885
# 11 -0.005301388 0.001609693 -3.293415 1.125530e-03 0.8521564    0.8574578 GO:0002253                                    activation of immune response   423 423/16644     423 0.0344821600
# 12 -0.007411402 0.002285888 -3.242242 1.339005e-03 0.9411178    0.9485292 GO:0060033                                  anatomical structure regression     7   7/16644       7 0.0376037119
# 13 -0.004705034 0.001490353 -3.156993 1.779746e-03 0.9213972    0.9261023 GO:0035036                                            sperm-egg recognition    24  24/16644      24 0.0428410170
# 14 -0.004894814 0.001548555 -3.160891 1.756965e-03 0.8541901    0.8590849 GO:0051091 positive regulation of DNA binding transcription factor activity   209 209/16644     209 0.0428410170
#
# $MF
#       Estimate  Std. Error   t value     Pr(>|t|) mean_sczd mean_control         ID                Description Count GeneRatio n_genes       padj
# 1 -0.010348218 0.003022956 -3.423211 0.0007175930 0.9006746    0.9110228 GO:0098748 endocytic adaptor activity    11  11/16644      11 0.02966208
# 2 -0.003167418 0.000905006 -3.499887 0.0005465153 0.8768780    0.8800454 GO:0048037           cofactor binding   401 401/16644     401 0.02966208
# 3 -0.005409429 0.001616273 -3.346854 0.0009366973 0.9463041    0.9517135 GO:0051920     peroxiredoxin activity     8   8/16644       8 0.02966208
save(indv_go_expr_sczd, indv_go_cleaned_sczd, file = 'rda/indv_go_sczd.Rdata')


## Load SCZD case-control results
load('rda/out_info.Rdata', verbose = TRUE)
compare_corr_de <- function(, cutvar = cutde) {
    out_comp <- out[ !is.na(out$adj.P.Val) & !is.na(out$expr.pval), ]
    res <- lapply(colnames(gene_pinfo)[1:6], function(var) {
        table(DE = out$adj.P.Val < cutde, Corr = out[, var] < cutvar)
    })
    names(res) <- colnames(gene_pinfo)[1:6]
    return(res)
}


computecor_sczd <- function(out, exp, cutde = 0.05) {
    
    genes <- unique(out$ensemblID[out$adj.P.Val < cutde])
    m <- match(genes, as.character(rowRanges(simple_rse[['DLPFC']][['gene']])$ensemblID))
    print(table(is.na(m)))
    m <- m[!is.na(m)]
    stopifnot(!any(is.na(m)))
    paircor(exp[['DLPFC']][['geneRpkm']][m, ], exp[['HIPPO']][['geneRpkm']][m, ])
}

indv_corrsczd_expr <- lapply(c(outGene, list('combined' = rbind(outGene[['HIPPO_matchQSV']], outGene[['DLPFC_matchQSV']]))), computecor_sczd, exp = expr)
#
# FALSE
#    48
#
# FALSE
#   245
#
# FALSE  TRUE
#   171    12
#
# FALSE  TRUE
#   339    37
#
# FALSE
#   293
indv_corrsczd_cleaned <- lapply(c(outGene, list('combined' = rbind(outGene[['HIPPO_matchQSV']], outGene[['DLPFC_matchQSV']]))), computecor_sczd, exp = cleaned)


corrsczd_info <- function(corrsczd, type) {
    dx <- colData(simple_rse[['DLPFC']][['gene']])$Dx
    dx <- factor(ifelse(dx == 'Schizo', 'SCZD', ifelse(dx == 'Control', 'Control', 'hmm')))
    # > levels(dx)
    # [1] "Control" "SCZD"
    
    res <- do.call(rbind, mapply(function(pathway, pathname) {
        f <- lm(pathway ~ dx)
        p <- summary(f)$coef[2, 4]
        ylim <- range(pathway)
        
        boxplot(pathway ~ dx, main = paste0(pathname, ': ',
            '\np-value: ', signif(p, 3)),
            ylim = ylim,
            xlab = 'SCZD diagnosis', outline = FALSE, ylab = paste('Correlation -', type))
        points(pathway ~ jitter(as.numeric(dx), amount = 0.15), cex = 1.5, pch = 21, bg = ifelse(type == 'expr', 'deepskyblue3', '#009E73'))
                
        c(summary(f)$coef[2, ], mean_sczd = mean(pathway[dx == 'SCZD']), mean_control = mean(pathway[dx == 'Control']))
    }, corrsczd, names(corrsczd), SIMPLIFY = FALSE))
    
    res <- as.data.frame(res)
    res$n_genes <- c(48, 245, 171, 339, 293)
    res$padj <- p.adjust(res[, 'Pr(>|t|)'], method = 'fdr')
    res <- res[order(res$padj, decreasing = FALSE), ]
    rownames(res) <- NULL
    
    return(res)
}

pdf('pdf/indv_box_corrsczd.pdf', useDingbats = FALSE)
indv_corrsczd_info_expr <- corrsczd_info(indv_corrsczd_expr, 'expr')
indv_corrsczd_info_cleaned <- corrsczd_info(indv_corrsczd_cleaned, 'cleaned')
dev.off()
options(width = 120)
indv_corrsczd_info_expr
#        Estimate  Std. Error    t value  Pr(>|t|) mean_sczd mean_control n_genes      padj
# 1 -0.0008518767 0.005349746 -0.1592368 0.8736045 0.9389143    0.9397661      48 0.8736045
# 2 -0.0014287886 0.005414014 -0.2639056 0.7920593 0.9124227    0.9138515     245 0.8736045
# 3 -0.0006810401 0.003531634 -0.1928399 0.8472332 0.9264122    0.9270933     171 0.8736045
# 4 -0.0026563724 0.005823791 -0.4561242 0.6486771 0.9146195    0.9172759     339 0.8736045
# 5 -0.0011000972 0.005267718 -0.2088375 0.8347367 0.9157445    0.9168446     293 0.8736045
indv_corrsczd_info_cleaned
#        Estimate  Std. Error    t value     Pr(>|t|) mean_sczd mean_control n_genes         padj
# 1 -0.0187940225 0.003019046 -6.2251529 1.891553e-09 0.7856714    0.8044655      48 9.457767e-09
# 2  0.0013918067 0.001547340  0.8994834 3.692180e-01 0.8623785    0.8609867     171 6.153634e-01
# 3 -0.0018085495 0.001784012 -1.0137539 3.116319e-01 0.8150913    0.8168999     293 6.153634e-01
# 4  0.0007077981 0.001791474  0.3950927 6.930947e-01 0.8221908    0.8214830     245 6.930947e-01
# 5 -0.0005068032 0.001219738 -0.4155016 6.781135e-01 0.8477722    0.8482790     339 6.930947e-01

save(indv_corrsczd_expr, indv_corrsczd_cleaned, indv_corrsczd_info_expr, indv_corrsczd_info_cleaned, file = 'rda/indv_corrsczd.Rdata')

## Re-loading if necessary
if(FALSE) {
    f <- dir('rda', full.names = TRUE)
    f <- f[!grepl('limma', f)]
    for(ff in f) load(ff, verbose = TRUE)
    rm(ff, f)
}

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
#  date     2018-08-23
#
# Packages --------------------------------------------------------------------------------------------------------------
#  package              * version   date       source
#  acepack                1.4.1     2016-10-29 CRAN (R 3.5.0)
#  AnnotationDbi          1.42.1    2018-05-17 Bioconductor
#  assertthat             0.2.0     2017-04-11 CRAN (R 3.5.0)
#  backports              1.1.2     2017-12-13 CRAN (R 3.5.0)
#  base                 * 3.5.0     2018-05-02 local
#  base64enc              0.1-3     2015-07-28 CRAN (R 3.5.0)
#  bibtex                 0.4.2     2017-06-30 CRAN (R 3.5.0)
#  bindr                  0.1.1     2018-03-13 CRAN (R 3.5.0)
#  bindrcpp               0.2.2     2018-03-29 CRAN (R 3.5.0)
#  Biobase              * 2.40.0    2018-05-02 Bioconductor
#  BiocGenerics         * 0.26.0    2018-05-03 Bioconductor
#  BiocParallel         * 1.14.2    2018-07-08 Bioconductor
#  biomaRt                2.36.1    2018-06-13 Bioconductor
#  Biostrings             2.48.0    2018-05-03 Bioconductor
#  bit                    1.1-14    2018-05-29 CRAN (R 3.5.0)
#  bit64                  0.9-7     2017-05-08 CRAN (R 3.5.0)
#  bitops                 1.0-6     2013-08-17 CRAN (R 3.5.0)
#  blob                   1.1.1     2018-03-25 CRAN (R 3.5.0)
#  BSgenome               1.48.0    2018-05-03 Bioconductor
#  bumphunter             1.22.0    2018-05-03 Bioconductor
#  checkmate              1.8.5     2017-10-24 CRAN (R 3.5.0)
#  cluster                2.0.7-1   2018-04-13 CRAN (R 3.5.0)
#  clusterProfiler      * 3.8.1     2018-06-13 Bioconductor
#  codetools              0.2-15    2016-10-05 CRAN (R 3.5.0)
#  colorout             * 1.2-0     2018-05-02 Github (jalvesaq/colorout@c42088d)
#  colorspace             1.3-2     2016-12-14 CRAN (R 3.5.0)
#  compiler               3.5.0     2018-05-02 local
#  cowplot                0.9.3     2018-07-15 CRAN (R 3.5.0)
#  crayon                 1.3.4     2017-09-16 CRAN (R 3.5.0)
#  data.table             1.11.4    2018-05-27 cran (@1.11.4)
#  datasets             * 3.5.0     2018-05-02 local
#  DBI                    1.0.0     2018-05-02 CRAN (R 3.5.0)
#  DelayedArray         * 0.6.2     2018-07-23 Bioconductor
#  derfinder            * 1.14.0    2018-05-03 Bioconductor
#  derfinderHelper        1.14.0    2018-05-03 Bioconductor
#  devtools             * 1.13.6    2018-06-27 CRAN (R 3.5.0)
#  digest                 0.6.15    2018-01-28 CRAN (R 3.5.0)
#  DO.db                  2.9       2018-05-03 Bioconductor
#  doRNG                  1.7.1     2018-06-22 CRAN (R 3.5.0)
#  DOSE                   3.6.1     2018-06-20 Bioconductor
#  dplyr                  0.7.6     2018-06-29 CRAN (R 3.5.0)
#  enrichplot             1.0.2     2018-06-13 Bioconductor
#  fastmatch              1.1-0     2017-01-28 CRAN (R 3.5.0)
#  fgsea                  1.6.0     2018-05-03 Bioconductor
#  foreach                1.4.4     2017-12-12 CRAN (R 3.5.0)
#  foreign                0.8-70    2017-11-28 CRAN (R 3.5.0)
#  Formula                1.2-3     2018-05-03 CRAN (R 3.5.0)
#  GenomeInfoDb         * 1.16.0    2018-05-03 Bioconductor
#  GenomeInfoDbData       1.1.0     2018-04-17 Bioconductor
#  GenomicAlignments      1.16.0    2018-05-03 Bioconductor
#  GenomicFeatures        1.32.0    2018-05-03 Bioconductor
#  GenomicFiles           1.16.0    2018-05-03 Bioconductor
#  GenomicRanges        * 1.32.6    2018-07-20 Bioconductor
#  ggforce                0.1.3     2018-07-07 CRAN (R 3.5.0)
#  ggplot2                3.0.0     2018-07-03 CRAN (R 3.5.0)
#  ggraph                 1.0.2     2018-07-07 CRAN (R 3.5.0)
#  ggrepel                0.8.0     2018-05-09 CRAN (R 3.5.0)
#  ggridges               0.5.0     2018-04-05 CRAN (R 3.5.0)
#  glue                   1.3.0     2018-07-17 CRAN (R 3.5.0)
#  GO.db                  3.6.0     2018-05-03 Bioconductor
#  GOSemSim               2.6.0     2018-05-03 Bioconductor
#  graphics             * 3.5.0     2018-05-02 local
#  grDevices            * 3.5.0     2018-05-02 local
#  grid                   3.5.0     2018-05-02 local
#  gridExtra              2.3       2017-09-09 CRAN (R 3.5.0)
#  gtable                 0.2.0     2016-02-26 CRAN (R 3.5.0)
#  Hmisc                  4.1-1     2018-01-03 CRAN (R 3.5.0)
#  hms                    0.4.2     2018-03-10 CRAN (R 3.5.0)
#  htmlTable              1.12      2018-05-26 CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 CRAN (R 3.5.0)
#  htmlwidgets            1.2       2018-04-19 CRAN (R 3.5.0)
#  httpuv                 1.4.5     2018-07-19 CRAN (R 3.5.0)
#  httr                   1.3.1     2017-08-20 CRAN (R 3.5.0)
#  igraph                 1.2.1     2018-03-10 CRAN (R 3.5.0)
#  IRanges              * 2.14.10   2018-05-17 Bioconductor
#  iterators              1.0.10    2018-07-13 CRAN (R 3.5.0)
#  jaffelab             * 0.99.21   2018-05-03 Github (LieberInstitute/jaffelab@7ed0ab7)
#  knitr                  1.20      2018-02-20 CRAN (R 3.5.0)
#  later                  0.7.3     2018-06-08 CRAN (R 3.5.0)
#  lattice                0.20-35   2017-03-25 CRAN (R 3.5.0)
#  latticeExtra           0.6-28    2016-02-09 CRAN (R 3.5.0)
#  lazyeval               0.2.1     2017-10-29 CRAN (R 3.5.0)
#  limma                  3.36.2    2018-06-21 Bioconductor
#  locfit                 1.5-9.1   2013-04-20 CRAN (R 3.5.0)
#  magrittr               1.5       2014-11-22 CRAN (R 3.5.0)
#  MASS                   7.3-50    2018-04-30 CRAN (R 3.5.0)
#  Matrix                 1.2-14    2018-04-13 CRAN (R 3.5.0)
#  matrixStats          * 0.54.0    2018-07-23 CRAN (R 3.5.0)
#  memoise                1.1.0     2017-04-21 CRAN (R 3.5.0)
#  methods              * 3.5.0     2018-05-02 local
#  munsell                0.5.0     2018-06-12 CRAN (R 3.5.0)
#  nnet                   7.3-12    2016-02-02 CRAN (R 3.5.0)
#  parallel             * 3.5.0     2018-05-02 local
#  pillar                 1.3.0     2018-07-14 CRAN (R 3.5.0)
#  pkgconfig              2.0.1     2017-03-21 CRAN (R 3.5.0)
#  pkgmaker               0.27      2018-05-25 CRAN (R 3.5.0)
#  plyr                   1.8.4     2016-06-08 CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 CRAN (R 3.5.0)
#  prettyunits            1.0.2     2015-07-13 CRAN (R 3.5.0)
#  progress               1.2.0     2018-06-14 CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 CRAN (R 3.5.0)
#  purrr                  0.2.5     2018-05-29 CRAN (R 3.5.0)
#  qvalue                 2.12.0    2018-05-03 Bioconductor
#  R6                     2.2.2     2017-06-17 CRAN (R 3.5.0)
#  rafalib              * 1.0.0     2015-08-09 CRAN (R 3.5.0)
#  RColorBrewer           1.1-2     2014-12-07 CRAN (R 3.5.0)
#  Rcpp                   0.12.18   2018-07-23 CRAN (R 3.5.0)
#  RCurl                  1.95-4.11 2018-07-15 CRAN (R 3.5.0)
#  registry               0.5       2017-12-03 CRAN (R 3.5.0)
#  reshape2               1.4.3     2017-12-11 CRAN (R 3.5.0)
#  rlang                  0.2.1     2018-05-30 cran (@0.2.1)
#  rmote                * 0.3.4     2018-05-02 deltarho (R 3.5.0)
#  rngtools               1.3.1     2018-05-15 CRAN (R 3.5.0)
#  rpart                  4.1-13    2018-02-23 CRAN (R 3.5.0)
#  Rsamtools              1.32.2    2018-07-03 Bioconductor
#  RSQLite                2.1.1     2018-05-06 CRAN (R 3.5.0)
#  rstudioapi             0.7       2017-09-07 CRAN (R 3.5.0)
#  rtracklayer            1.40.3    2018-06-13 Bioconductor
#  rvcheck                0.1.0     2018-05-23 CRAN (R 3.5.0)
#  S4Vectors            * 0.18.3    2018-06-13 Bioconductor
#  scales               * 0.5.0     2017-08-24 CRAN (R 3.5.0)
#  segmented              0.5-3.0   2017-11-30 CRAN (R 3.5.0)
#  servr                  0.10      2018-05-30 CRAN (R 3.5.0)
#  splines                3.5.0     2018-05-02 local
#  stats                * 3.5.0     2018-05-02 local
#  stats4               * 3.5.0     2018-05-02 local
#  stringi                1.2.4     2018-07-20 CRAN (R 3.5.0)
#  stringr                1.3.1     2018-05-10 CRAN (R 3.5.0)
#  SummarizedExperiment * 1.10.1    2018-05-17 Bioconductor
#  survival               2.42-3    2018-04-16 CRAN (R 3.5.0)
#  tibble                 1.4.2     2018-01-22 CRAN (R 3.5.0)
#  tidyr                  0.8.1     2018-05-18 CRAN (R 3.5.0)
#  tidyselect             0.2.4     2018-02-26 CRAN (R 3.5.0)
#  tools                  3.5.0     2018-05-02 local
#  tweenr                 0.1.5     2016-10-10 CRAN (R 3.5.0)
#  units                  0.6-0     2018-06-09 CRAN (R 3.5.0)
#  UpSetR                 1.3.3     2017-03-21 CRAN (R 3.5.0)
#  utils                * 3.5.0     2018-05-02 local
#  VariantAnnotation      1.26.1    2018-07-04 Bioconductor
#  viridis                0.5.1     2018-03-29 CRAN (R 3.5.0)
#  viridisLite            0.3.0     2018-02-01 CRAN (R 3.5.0)
#  withr                  2.1.2     2018-03-15 CRAN (R 3.5.0)
#  xfun                   0.3       2018-07-06 CRAN (R 3.5.0)
#  XML                    3.98-1.12 2018-07-15 CRAN (R 3.5.0)
#  xtable                 1.8-2     2016-02-05 CRAN (R 3.5.0)
#  XVector                0.20.0    2018-05-03 Bioconductor
#  zlibbioc               1.26.0    2018-05-02 Bioconductor
