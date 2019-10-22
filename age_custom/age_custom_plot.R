library('limma')
library('SummarizedExperiment')
library('jaffelab')
library('sessioninfo')
library('ggplot2')
library('getopt')
library('data.table')

## Specify parameters
spec <- matrix(c(
    'type', 't', 1, 'character', 'Either gene, exon, tx or jxn',
    'symbol', 's', 1, 'character', 'Symbol for the gene',
    'outdir', 'o', 1, 'character', 'Output directory: pdf',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## For testing
if(FALSE){
    opt <- list(type = 'gene', symbol = 'BDNF', outdir = 'pdf')
    opt <- list(type = 'exon', symbol = 'BDNF', outdir = 'pdf')
}

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}


stopifnot(opt$type %in% c('gene', 'exon', 'tx', 'jxn'))




source('/dcl01/lieber/ajaffe/lab/brainseq_phase2/development/load_funs.R')
dir.create(opt$outdir, showWarnings = FALSE)

## Load required data
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/development/rda/pcheck_both.Rdata', verbose = TRUE)
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/development/rda/raw.Rdata', verbose = TRUE)
top <- lapply(raw, '[[', 'top')
fit <- lapply(raw, '[[', 'fit')
exprsNorm <- lapply(raw, '[[', 'exprsNorm')

rse <- load_foo(opt$type)


if(opt$type == 'exon') {
    load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/browser/rda/exon_name_map.Rdata', verbose = TRUE)
}

## Reg specific model
design <- get_mods( colData(rse), int = TRUE)$mod


pinfo <- subset(pcheck_both, type == opt$type)

# pinfo <- pinfo[order(pinfo$P.Bonf), ]
# pinfo <- pinfo[sign(pinfo$F) == sign(pinfo$span_F) & pinfo$span_P.Value < 0.05 & pinfo$P.Bonf < 0.01, ]
#
#

find_feature <- function() {
    if(opt$type %in% c('gene', 'exon', 'jxn')) {
        res <- which(opt$symbol == rowRanges(rse)$Symbol)
    } else {
        res <- which(opt$symbol == rowRanges(rse)$gene_name)
    }
    if(length(res) == 0) stop('Could not find', opt$symbol)
    return(res)
}

## Find the exact feature(s) to plot
to_plot <- paste0(opt$type, '.', rownames(rse)[find_feature()])
pinfo <- pinfo[to_plot, ]

tocheck <- match(gsub(paste0(opt$type, '.'), '', head(rownames(pinfo), 100)), rownames(rse))
tocheck


## Re-organize design before cleaningY
protect <- grepl(':RegionHIPPO|RegionHIPPO:|Intercept', colnames(design))
design <- cbind(design[, protect], design[, !protect])

cleanedVoom <- cleaningY(exprsNorm[[opt$type]], design, sum(protect))


get_main <- function(i) {
    j <- tocheck[i]
    if(opt$type %in% c('gene', 'exon')) {
        res <- with(rowRanges(rse)[j, ], paste(gencodeID, Symbol))
    } else if(opt$type == 'jxn') {
        res <- with(rowRanges(rse)[j, ], paste(gencodeGeneID, Symbol))
    } else {
        res <- with(rowRanges(rse)[j, ], paste(gene_id, gene_name))
    }
    res <- paste(res, 'p-bonf', signif(pinfo$P.Bonf[i], 3))
    if(opt$type != 'gene') {
        res <- paste(rownames(rse)[j], res)
    } 
    if (opt$type == 'exon') {
        ## Add Gencode exon id
        res <- paste(exon_name_map$gencode[exon_name_map$libd_bsp2 == rownames(rse)[j]], res)
    }
    return(res)
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

plot_age_mod <-
    design[, c(
        '(Intercept)',
        'Age',
        'RegionHIPPO',
        'fetal',
        'birth',
        'infant',
        'child',
        'teen',
        'adult'
    )]
p_cols <- ifelse(colData(rse)$Region == 'HIPPO', 'skyblue3', 'dark orange')
l_cols <- c('lightgoldenrod', 'light blue')
age_brks <- c(-1, 0, 1, 10, 20, 50, 100)

pdf(file.path(opt$outdir, paste0('bsp2_age_', opt$symbol, '_', opt$type, '_norm.pdf')), width = 14, useDingbats = FALSE)
for(i in seq_len(length(tocheck))) {
    set.seed(20180419)
    agePlotter(
        y = exprsNorm[[opt$type]][tocheck[i],],
        age = colData(rse)$Age,
        pointColor = p_cols,
        ageBreaks = age_brks,
        mainText = get_main(i),
        lineColor = l_cols,
        mod = plot_age_mod,
        ylab = get_ylab(opt$type)
    )
    legend('bottom', c('DLPFC', 'HIPPO'), col = l_cols, lwd = 3, bty = 'n', ncol = 2)
}
dev.off()

pdf(file.path(opt$outdir, paste0('bsp2_age_', opt$symbol, '_', opt$type, '_cleaned.pdf')), width = 14, useDingbats = FALSE)
for(i in seq_len(length(tocheck))) {
    set.seed(20180419)
    agePlotter(
        cleanedVoom[tocheck[i],],
        colData(rse)$Age,
        pointColor = p_cols,
        ageBreaks = age_brks,
        mainText = get_main(i),
        lineColor = l_cols,
        mod = plot_age_mod,
        ylab = get_ylab(opt$type, cleaned = TRUE)
    )
    legend('bottom', c('DLPFC', 'HIPPO'), col = l_cols, lwd = 3, bty = 'n', ncol = 2)
}
dev.off()


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
