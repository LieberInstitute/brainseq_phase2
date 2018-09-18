library('GenomicRanges')
library('data.table')
library('devtools')

## Load subsets of data
files_sub <- dir('rdas', pattern = '_compare_', full.names = TRUE)
stopifnot(length(files_sub) == 12)
for(f in files_sub) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
}

x <- paste0(d_sig_genes$snps, '_', d_sig_genes$gene)
y <- paste0(dlpfc_gtex_genes$snps, '_', dlpfc_gtex_genes$gene)

x2 <- paste0(hmm$snps, '_', hmm$gene)

Sys.time()

length(x)
length(y)
identical(x, y)

identical(x2, y)

m <- match(x, y)
stopifnot(all(!is.na(m)))

dlpfc_gtex_genes[m]


head(x)
head(y)

head(dlpfc_gtex_genes)
head(dlpfc_gtex_genes$snps)

'10:100104794:ACTAAAATAAGT:A_ENSG00000107521.18' %in% x

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
