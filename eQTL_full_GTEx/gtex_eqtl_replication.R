library('data.table')
library('devtools')

## Load subsets of data
files_sub <- dir('rdas', pattern = '_compare_', full.names = TRUE)
stopifnot(length(files_sub) == 12)
for(f in files_sub) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
}

## For testing
# gtex <- inter_gtex_genes
# brainseq <- i_sig_genes

merge_qtl <- function(gtex, brainseq) {
    ## All should be in the same order
    stopifnot(identical(gtex$snps, brainseq$snps))
    stopifnot(identical(gtex$gene, brainseq$gene))
    
    ## Use the latest Symbol info from GTEx files
    brainseq$Symbol <- gtex$Symbol
    
    ## Keep only a few GTEx columns
    to_add <- gtex[, c('statistic', 'pvalue', 'FDR', 'beta')]
    colnames(to_add) <- paste0('gtex_', colnames(to_add))
    
    ## Combine
    cbind(brainseq, to_add)
}

## Merge
message(paste(Sys.time(), 'merging interaction QTLs'))
interaction <- list(
    'gene' = merge_qtl(inter_gtex_genes, i_sig_genes),
    'exon' = merge_qtl(inter_gtex_exons, i_sig_exons),
    'jxn' = merge_qtl(inter_gtex_jxns, i_sig_jxns),
    'tx' = merge_qtl(inter_gtex_txs, i_sig_txs)
)
#interaction_all <- do.call(rbind, interaction)
message(paste(Sys.time(), 'saving merged QTLs'))
save(interaction, file = 'rdas/merged_GTEx_BrainSeq_QTLs_interaction.Rdata')

message(paste(Sys.time(), 'merging hippo QTLs'))
hippo <- list(),
    'gene' = merge_qtl(hippo_gtex_genes, h_sig_genes),
    'exon' = merge_qtl(hippo_gtex_exons, h_sig_exons),
    'jxn' = merge_qtl(hippo_gtex_jxns, h_sig_jxns),
    'tx' = merge_qtl(hippo_gtex_txs, h_sig_txs)
)
#hippo_all <- do.call(rbind, hippo)
message(paste(Sys.time(), 'saving merged QTLs'))
save(hippo, file = 'rdas/merged_GTEx_BrainSeq_QTLs_hippo.Rdata')


message(paste(Sys.time(), 'merging dlpfc QTLs'))
dlpfc <- list(
    'gene' = merge_qtl(dlpfc_gtex_genes, d_sig_genes),
    'exon' = merge_qtl(dlpfc_gtex_exons, d_sig_exons),
    'jxn' = merge_qtl(dlpfc_gtex_jxns, d_sig_jxns),
    'tx' = merge_qtl(dlpfc_gtex_txs, d_sig_txs)
)
#interaction_all <- do.call(rbind, dlpfc)
message(paste(Sys.time(), 'saving merged QTLs'))
save(dlpfc, file = 'rdas/merged_GTEx_BrainSeq_QTLs_dlpfc.Rdata')


## Explore
comp_qtl <- function(type, dfs, perc = FALSE) {
    df <- dfs[[type]]
    if(any(is.na(df$gtex_statistic)) & !perc) {
        message(paste(Sys.time(), 'removing some NAs from GTEx (TRUEs below) for type', type))
        print(table(is.na(df$gtex_statistic)))
    }
    
    res <- addmargins(table('Equal sign' = sign(df$statistic) == sign(df$gtex_statistic), 'GTEx p<0.01' = df$gtex_pvalue < 0.01))
    if(!perc) return(res)
    
    ## Calculate percent over all of brainseq
    ## the total marginal will not be 100% unless there were no NAs
    res / nrow(df) * 100
}
comp_qtl_short <- function(dfs, perc = FALSE) {
    res <- lapply(names(dfs), comp_qtl, dfs = dfs, perc = perc)
    names(res) <- names(dfs)
    return(res)
}

comp_qtl_short(interaction)
# 2018-09-18 17:13:37 removing some NAs from GTEx (TRUEs below) for type gene
#
# FALSE  TRUE
# 35958  4134
# 2018-09-18 17:13:37 removing some NAs from GTEx (TRUEs below) for type exon
#
# FALSE  TRUE
# 79435 10488
# 2018-09-18 17:13:37 removing some NAs from GTEx (TRUEs below) for type jxn
#
# FALSE  TRUE
# 66829  8774
# 2018-09-18 17:13:37 removing some NAs from GTEx (TRUEs below) for type tx
#
# FALSE  TRUE
# 19658  2237
# $gene
#           GTEx p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE  8613   334  8947
#      TRUE  20852  6159 27011
#      Sum   29465  6493 35958
#
# $exon
#           GTEx p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE 17026   759 17785
#      TRUE  39580 22070 61650
#      Sum   56606 22829 79435
#
# $jxn
#           GTEx p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE 24361   639 25000
#      TRUE  32179  9650 41829
#      Sum   56540 10289 66829
#
# $tx
#           GTEx p<0.01
# Equal sign FALSE  TRUE   Sum
#      FALSE  2774   234  3008
#      TRUE  11077  5573 16650
#      Sum   13851  5807 19658

comp_qtl_short(interaction, perc = TRUE)
# $gene
#           GTEx p<0.01
# Equal sign      FALSE       TRUE        Sum
#      FALSE 21.4830889  0.8330839 22.3161728
#      TRUE  52.0103761 15.3621670 67.3725432
#      Sum   73.4934650 16.1952509 89.6887160
#
# $exon
#           GTEx p<0.01
# Equal sign      FALSE       TRUE        Sum
#      FALSE 18.9339768  0.8440555 19.7780323
#      TRUE  44.0154354 24.5432203 68.5586557
#      Sum   62.9494123 25.3872758 88.3366881
#
# $jxn
#           GTEx p<0.01
# Equal sign      FALSE       TRUE        Sum
#      FALSE 32.2222663  0.8452046 33.0674709
#      TRUE  42.5631258 12.7640438 55.3271696
#      Sum   74.7853921 13.6092483 88.3946404
#
# $tx
#           GTEx p<0.01
# Equal sign     FALSE      TRUE       Sum
#      FALSE 12.669559  1.068737 13.738296
#      TRUE  50.591459 25.453300 76.044759
#      Sum   63.261018 26.522037 89.783055


comp_qtl_short(hippo)
comp_qtl_short(hippo, perc = TRUE)
comp_qtl_short(dlpfc)
comp_qtl_short(dlpfc, perc = TRUE)


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
