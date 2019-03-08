library('SummarizedExperiment')
library('data.table')
library('devtools')


## Load and extract BSP1 full eQTL results for the browser
dir.create('rda', showWarnings = FALSE)
load('rda/exon_name_map.Rdata', verbose = TRUE)
setkey(exon_name_map, libd_bsp2)


f <- c('dlpfc_replication_bsp1' = '/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/eqtl/full/eqtl_tables/mergedEqtl_output_dlpfc_4features_in_progress.rda')
message(paste(Sys.time(), 'loading', f))
load(f, verbose=TRUE)

## Keep only some columns: leave the annotation to the other tables
message(paste(Sys.time(), 'processing gene'))
DT_gene <- data.table(as.data.frame(geneEqtl[, c('snps', 'gene', 'statistic', 'pvalue', 'FDR', 'beta', 'Type')]))
rm(geneEqtl)

message(paste(Sys.time(), 'processing exon'))
DT_exon <- data.table(as.data.frame(exonEqtl[, c('snps', 'gene', 'statistic', 'pvalue', 'FDR', 'beta', 'Type')]))
rm(exonEqtl)

message(paste(Sys.time(), 'fixing exon ids'))
## Use data.table syntax ^^
DT_exon$gene <- exon_name_map[.(DT_exon$gene), gencode]

message(paste(Sys.time(), 'processing jxn'))
DT_jxn <- data.table(as.data.frame(jxnEqtl[, c('snps', 'gene', 'statistic', 'pvalue', 'FDR', 'beta', 'Type')]))
rm(jxnEqtl)

message(paste(Sys.time(), 'processing tx'))
DT_tx <- data.table(as.data.frame(txEqtl[, c('snps', 'gene', 'statistic', 'pvalue', 'FDR', 'beta', 'Type')]))
rm(txEqtl)

DT <- rbind(DT_gene, DT_exon, DT_jxn, DT_tx)
rm(DT_gene, DT_exon, DT_jxn, DT_tx)

## So the column matches the snpAnno column name
colnames(DT)[1] <- 'snp'

## So the column matches the feature annotation column
colnames(DT)[2] <- 'feature_id'

## Make Type lowercase to match file names from other tables
DT$Type <- tolower(DT$Type)

f_new <- paste0('BrainSeqPhaseII_eQTL_', names(f), '.txt')
message(paste(Sys.time(), 'writing', f_new))
fwrite(DT, file = f_new, sep = '\t', row.names = FALSE)

message(paste(Sys.time(), 'lines for each file'))
system('wc -l BrainSeqPhaseII*')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
