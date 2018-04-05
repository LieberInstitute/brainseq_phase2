library('limma')
library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('ggplot2')
library('getopt')

## Specify parameters
spec <- matrix(c(
    'type', 't', 1, 'character', 'Either gene, exon, tx or jxn',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## For testing
if(FALSE){
    opt <- list('type' = 'gene')
}

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

source('load_funs.R')
dir.create('rda', showWarnings = FALSE)
dir.create('pdf', showWarnings = FALSE)

## Load required data
load('rda/pcheck_both.Rdata', verbose = TRUE)
load('rda/raw.Rdata', verbose = TRUE)
top <- lapply(raw, '[[', 'top')
fit <- lapply(raw, '[[', 'fit')
exprsNorm <- lapply(raw, '[[', 'exprsNorm')

rse <- load_foo(opt$type)

## Reg specific model
design <- get_mods( colData(rse), int = TRUE)$mod


pinfo <- subset(pcheck_both, type == opt$type)
pinfo <- pinfo[order(pinfo$P.Bonf), ]
pinfo <- pinfo[sign(pinfo$F) == sign(pinfo$span_F) & pinfo$span_P.Value < 0.05 & pinfo$P.Bonf < 0.01, ]


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
    return(res)
}

get_ylab <- function(type, cleaned = FALSE) {
    if(type %in% c('gene', 'exon', 'jxn')) {
        res <- 'Voom-normalized expression'
    } else {
        res <- 'TPM'
    }

    if(cleaned) res <- paste(res, 'with covariate effects removed')
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
p_cols <- ifelse(colData(rse)$Region == 'HIPPO', 'blue', 'orange')
l_cols <- c('lightgoldenrod', 'light blue')
age_brks <- c(-1, 0, 1, 10, 20, 50, 100)

pdf('pdf/top_gene_replicated_exprNorm.pdf', width = 14, useDingbats = FALSE)
for(i in 1:100) {
    agePlotter(
        exprsNorm$gene[tocheck[i],],
        colData(rse)$Age,
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

pdf('pdf/top_gene_replicated_cleanedVoom.pdf', width = 14, useDingbats = FALSE)
for(i in 1:100) {
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
