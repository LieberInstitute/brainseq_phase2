library('limma')
library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('ggplot2')
library('getopt')

## Specify parameters
spec <- matrix(c(
    'type', 't', 1, 'character', 'Either gene, exon, tx or jxn',
    'age', 'a', 1, 'character', 'Either adult or fetal',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## For testing
if(FALSE){
    opt <- list('type' = 'gene', 'age' = 'fetal')
}

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

stopifnot(opt$age %in% c('adult', 'fetal'))

source('load_funs.R')
dir.create('rda', showWarnings = FALSE)
dir.create('pdf', showWarnings = FALSE)


## Reg specific model
rse <- load_foo('gene', 'adult')
design <- get_mods( colData(rse) )$mod

min(top$adult_gene$adj.P.Val)
which.min(top$adult_gene$adj.P.Val)
which(rank(top$adult_gene$adj.P.Val) == 1)

cleanedVoom <- cleaningY(exprsNorm$adult_gene, design, 2)


pdf('pdf/top100_hits_adult_gene.pdf', useDingbats = FALSE)
for(j in seq_len(100)) {
    i <- seq_len(nrow(top$adult_gene))[order(top$adult_gene$adj.P.Val)][j]
    set.seed(20180319)
    boxplot(cleanedVoom[i, ] ~ colData(rse)$Region, ylab = 'Norm. Expr - adj covariates removed', main = paste(rownames(top$adult_gene)[i], rowRanges(rse)$Symbol[rowRanges(rse)$gencodeID == rownames(top$adult_gene)[i]], 'FDR', signif(top$adult_gene$adj.P.Val[i], 3)), col = c('lightgoldenrod', 'light blue'), ylim = abs(range(cleanedVoom[i, ])) * 1.05 * sign(range(cleanedVoom[i, ])), outline = FALSE)
    points(cleanedVoom[i, ] ~ jitter(as.integer(colData(rse)$Region), 1), pch = 21, bg = c('darkgoldenrod2', 'steelblue1')[as.integer(colData(rse)$Region)])
}
dev.off()

rse <- load_foo('gene', 'fetal')
design <- get_mods( colData(rse) )$mod
cleanedVoom <- cleaningY(exprsNorm$fetal_gene, design, 2)

pdf('pdf/top100_hits_fetal_gene.pdf', useDingbats = FALSE)
for(j in seq_len(100)) {
    i <- seq_len(nrow(top$fetal_gene))[order(top$fetal_gene$adj.P.Val)][j]
    set.seed(20180319)
    boxplot(cleanedVoom[i, ] ~ colData(rse)$Region, ylab = 'Norm. Expr - adj covariates removed', main = paste(rownames(top$fetal_gene)[i], rowRanges(rse)$Symbol[rowRanges(rse)$gencodeID == rownames(top$fetal_gene)[i]], 'FDR', signif(top$fetal_gene$adj.P.Val[i], 3)), col = c('lightgoldenrod', 'light blue'), ylim = abs(range(cleanedVoom[i, ])) * 1.05 * sign(range(cleanedVoom[i, ])), outline = FALSE)
    points(cleanedVoom[i, ] ~ jitter(as.integer(colData(rse)$Region), 1), pch = 21, bg = c('darkgoldenrod2', 'steelblue1')[as.integer(colData(rse)$Region)])
}
dev.off()





top$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 1), ]
boxplot(exprsNorm$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 1), ] ~ colData(rse)$Region)
t.test(exprsNorm$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 1), ] ~ colData(rse)$Region)

boxplot(cleanedVoom[which(rank(top$adult_gene$adj.P.Val) == 1), ] ~ colData(rse)$Region)
t.test(cleanedVoom[which(rank(top$adult_gene$adj.P.Val) == 1), ] ~ colData(rse)$Region)

top$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 2), ]
boxplot(exprsNorm$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 2), ] ~ colData(rse)$Region)
t.test(exprsNorm$adult_gene[which(rank(top$adult_gene$adj.P.Val) == 2), ] ~ colData(rse)$Region)

boxplot(cleanedVoom[which(rank(top$adult_gene$adj.P.Val) == 2), ] ~ colData(rse)$Region)
t.test(cleanedVoom[which(rank(top$adult_gene$adj.P.Val) == 2), ] ~ colData(rse)$Region)

boxplot(assays(rse)$rpkm[which(rank(top$adult_gene$adj.P.Val) == 2), ] ~ colData(rse)$Region)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
