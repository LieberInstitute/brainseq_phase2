# Usage:
# qrsh -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
# module load conda_R/3.4.x
# Rscript subset_span.R  > subset_span_log.txt 2>&1

library('SummarizedExperiment')
library('devtools')

## Locate and load all rse files
rse_files <- c(dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff',
    pattern = 'rse_[[:alpha:]]*.Rdata', full.names = TRUE),
    dir('/dcl01/lieber/ajaffe/lab/brainspan_analysis/rse', pattern = 'Rdata$',
    full.names = TRUE))
stopifnot(length(rse_files) == 8)

for(f in rse_files) system.time(load(f, verbose = TRUE))
rm(f)

## Subset each rse_span object by the genes expressed in brainseq phase 2
span_subset <- function(rse, rse_span, exon = FALSE) {
    ## Match by feature name
    if(exon) {
        ov <- findOverlaps(rowRanges(rse_exon), rowRanges(rse_span_exon), type = 'equal', ignore.strand = FALSE)
        m <- subjectHits(ov)
    } else {
        m <- match(rownames(rse), rownames(rse_span))
    }
    print('Missing features')
    print(table(is.na(m)))
    
    ## Subset and keep only the DLPFC and HIPPO regions
    rse_span <- rse_span[m, rse_span$Regioncode %in% c('DFC', 'HIP')]
    
    ## Set as factor
    colData(rse_span)$Region <- relevel(factor(colData(rse_span)$Regioncode), 'DFC')
    race <- colData(rse_span)$Ethnicity
    race[race == 'European'] <- 'CAUC'
    colData(rse_span)$Race <- relevel(factor(race), ref = 'CAUC')
    colData(rse_span)$Sex <- relevel(factor(colData(rse_span)$Sex), ref = 'F')
    
    ## Add means (these are just the same numbers, but for consistency with
    ## the rest of the BrainSeq code I'll call them means)
    colData(rse_span)$mean_mitoRate <- colData(rse_span)$mitoRate
    colData(rse_span)$mean_totalAssignedGene <- colData(rse_span)$totalAssignedGene
    colData(rse_span)$mean_rRNA_rate <- colData(rse_span)$rRNA_rate
    colData(rse_span)$mean_RIN <- colData(rse_span)$RIN
    
    print('Dimensions of the data')
    print(dim(rse_span))
    return(rse_span)
}

## Actually perform the subset
rse_span_gene <- span_subset(rse_gene, rse_span_gene)
rse_span_exon <- span_subset(rse_exon, rse_span_exon, exon = TRUE)
## Note here that the BrainSpan data set is unstranded while the BrainSeq phase 2 data is stranded...
table(strand(rse_jxn))
table(strand(rse_span_jxn))
rse_span_jxn <- span_subset(rse_jxn, rse_span_jxn, exon = TRUE)
table(strand(rse_span_jxn))
rse_span_tx <- span_subset(rse_tx, rse_span_tx)

## Save filtered results
save(rse_span_gene, file = 'rse_span_gene.Rdata')
save(rse_span_exon, file = 'rse_span_exon.Rdata')
save(rse_span_jxn, file = 'rse_span_jxn.Rdata')
save(rse_span_tx, file = 'rse_span_tx.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
