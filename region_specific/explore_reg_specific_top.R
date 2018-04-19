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

## Load required data
load('rda/raw.Rdata', verbose = TRUE)
top <- lapply(raw, '[[', 'top')
fit <- lapply(raw, '[[', 'fit')
exprsNorm <- lapply(raw, '[[', 'exprsNorm')

## Reg specific model
rse <- load_foo(opt$type, opt$age)
design <- get_mods( colData(rse) )$mod

cleanedVoom <- cleaningY(exprsNorm[[paste0(opt$age, '_', opt$type)]], design, 2)

get_main <- function(i) {
    if(opt$type %in% c('gene', 'exon')) {
        var <- 'gencodeID'
        vars <- 'Symbol'
    } else if (opt$type == 'jxn') {
        var <- 'gencodeGeneID'
        vars <- 'Symbol'
    } else {
        var <- 'gene_id'
        vars <- 'gene_name'
    }

    topnow <- top[[paste0(opt$age, '_', opt$type)]]

    paste(rownames(topnow)[i],
          mcols(rowRanges(rse))[, vars][mcols(rowRanges(rse))[, var] == rownames(topnow)[i]],
          'FDR',
          signif(topnow$adj.P.Val[i], 3))
}

get_ylim_mult <- function(rang) {
    c(
        ifelse(sign(rang[1]) == 1, 0.95, 1.05),
        ifelse(sign(rang[2]) == 1, 1.05, 0.95)
    )
}


make_plot <- function(df, i, ylab = 'Norm. Expr - adj covariates removed') {
    set.seed(20180319)
    boxplot(df[i, ] ~ colData(rse)$Region, ylab = ylab, main = get_main(i),
            col = c('darkgoldenrod2', 'steelblue1'),
            ylim = abs(range(df[i, ])) * get_ylim_mult(range(df[i, ])) * sign(range(df[i, ])),
            outline = FALSE)
    points(df[i, ] ~ jitter(as.integer(colData(rse)$Region), 1), pch = 21,
           bg = c('dark orange', 'skyblue3')[as.integer(colData(rse)$Region)])
}


pdf(paste0('pdf/top100_hits_', opt$age, '_', opt$type, '_cleaned.pdf'), useDingbats = FALSE)
for(j in seq_len(100)) {
    i <- seq_len(nrow(top[[paste0(opt$age, '_', opt$type)]]))[order(top[[paste0(opt$age, '_', opt$type)]]$adj.P.Val)][j]
    make_plot(cleanedVoom, i)
}
dev.off()

pdf(paste0('pdf/top100_hits_', opt$age, '_', opt$type, '_norm.pdf'), useDingbats = FALSE)
for(j in seq_len(100)) {
    i <- seq_len(nrow(top[[paste0(opt$age, '_', opt$type)]]))[order(top[[paste0(opt$age, '_', opt$type)]]$adj.P.Val)][j]
    make_plot(exprsNorm[[paste0(opt$age, '_', opt$type)]], i, 'Norm. Expr')
}
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
