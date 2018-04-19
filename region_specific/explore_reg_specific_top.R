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
    opt <- list('type' = 'exon', 'age' = 'fetal')
    opt <- list('type' = 'jxn', 'age' = 'fetal')
    opt <- list('type' = 'tx', 'age' = 'fetal')
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
load('rda/pcheck_both.Rdata', verbose = TRUE)
load('rda/raw.Rdata', verbose = TRUE)
top <- lapply(raw, '[[', 'top')
fit <- lapply(raw, '[[', 'fit')
exprsNorm <- lapply(raw, '[[', 'exprsNorm')

## Reg specific model
rse <- load_foo(opt$type, opt$age)
design <- get_mods( colData(rse) )$mod

cleanedVoom <- cleaningY(exprsNorm[[paste0(opt$age, '_', opt$type)]], design, 2)


pinfo <- subset(pcheck_both, type == opt$type & age == opt$age)
pinfo <- pinfo[order(pinfo$P.Bonf), ]
pinfo <- pinfo[sign(pinfo$t) == sign(pinfo$span_t) & pinfo$span_P.Value < 0.05 & pinfo$P.Bonf < 0.01, ]

tocheck <- match(gsub(paste0(opt$age, '_', opt$type, '.'), '', head(rownames(pinfo), 100)), rownames(rse))
tocheck

if(length(tocheck) == 0) {
    print('no results to plot')
    ## Reproducibility information
    print('Reproducibility information:')
    print(Sys.time())
    print(proc.time())
    options(width = 120)
    print(session_info())
    q('n')
}

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

    k <- which(names(rowRanges(rse)) == rownames(topnow)[i])

    paste(if(opt$type != 'jxn') mcols(rowRanges(rse))[, var][k] else rownames(topnow)[i],
          if(is.na(mcols(rowRanges(rse))[, vars][k])) '' else mcols(rowRanges(rse))[, vars][k],
          'FDR',
          signif(topnow$adj.P.Val[i], 3))
}

get_ylim_mult <- function(rang) {
    c(
        ifelse(sign(rang[1]) == 1, 0.95, 1.05),
        ifelse(sign(rang[2]) == 1, 1.05, 0.95)
    )
}

get_ylab <- function(type, cleaned = FALSE) {
    if(type %in% c('gene', 'exon', 'jxn')) {
        res <- 'log2(CPM + 0.5)'
    } else {
        res <- 'log2(TPM + 0.5)'
    }

    if(cleaned) res <- paste(res, '- covariate effects removed')
    return(res)
}


make_plot <- function(df, i, clean = TRUE) {
    set.seed(20180319)
    boxplot(df[i, ] ~ colData(rse)$Region,
            ylab = get_ylab(opt$type, clean), main = get_main(i),
            col = c('darkgoldenrod2', 'steelblue1'),
            ylim = abs(range(df[i, ])) * get_ylim_mult(range(df[i, ])) * sign(range(df[i, ])),
            outline = FALSE)
    points(df[i, ] ~ jitter(as.integer(colData(rse)$Region), 1), pch = 21,
           bg = c('dark orange', 'skyblue3')[as.integer(colData(rse)$Region)])
}


pdf(paste0('pdf/tophits_', opt$age, '_', opt$type, '_cleaned.pdf'), useDingbats = FALSE)
for(j in seq_len(length(tocheck))) {
    make_plot(cleanedVoom, i = tocheck[j])
}
dev.off()

pdf(paste0('pdf/tophits_', opt$age, '_', opt$type, '_norm.pdf'), useDingbats = FALSE)
for(j in seq_len(length(tocheck))) {
    make_plot(exprsNorm[[paste0(opt$age, '_', opt$type)]], i = tocheck[j], clean = FALSE)
}
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
