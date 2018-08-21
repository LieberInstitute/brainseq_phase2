# Usage:
# qrsh -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
# module load conda_R/3.5
# Rscript subset_gtex.R  > logs/subset_gtex.txt 2>&1

library('SummarizedExperiment')
library('devtools')
library('recount')
library('jaffelab')

## Locate and load all rse files
rse_files <- c(dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff',
    pattern = 'rse_[[:alpha:]]*.Rdata', full.names = TRUE),
    dir(pattern = 'Rdata$', full.names = TRUE))
stopifnot(length(rse_files) == 8)

## Load GTEx first
for(f in rse_files[5:8]) system.time(load(f, verbose = TRUE))

## Add GTEx metadata
meta <- all_metadata('gtex')

## Downlaoded subject info from https://www.gtexportal.org/home/datasets
meta_subj <- read.table('../gtex/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt', header = TRUE, sep = '\t')
colnames(meta_subj) <- tolower(colnames(meta_subj))

## Simplify the age var
meta_subj$agegrp <- meta_subj$age
meta_subj$age <- mapply(function(x, y) { mean(c(x, y)) },
    as.integer(ss(as.character(meta_subj$agegrp), '-', 1)),
    as.integer(ss(as.character(meta_subj$agegrp), '-', 2)))

## Simplify gender using info from GTEx_Data_V6_Annotations_SubjectPhenotypes_DD.xlsx
stopifnot(all(meta_subj$gender %in% 1:2))
meta_subj$gender <- relevel(factor(ifelse(meta_subj$gender == 1, 'M', 'F')), ref = 'F')

## Add to recount2 GTEx metadata
m <- match(paste0(ss(meta$sampid, '-', 1), '-', ss(meta$sampid, '-', 2)), meta_subj$subjid)
meta <- cbind(meta, meta_subj[m, ])

## Add all the data
add_gtex_meta <- function(rse) {
    m <- match(rse$SAMPLE_ID, meta$run)
    colData(rse) <- cbind(colData(rse), meta[m, ])
    return(add_predictions(rse))
}
rse_gtex_gene <- add_gtex_meta(rse_gene)
rse_gtex_exon <- add_gtex_meta(rse_exon)
rse_gtex_jxn <- add_gtex_meta(rse_jx)
rse_gtex_tx <- add_gtex_meta(rse_tx)

## Now load BrainSeq Phase II
for(f in rse_files[1:4]) system.time(load(f, verbose = TRUE))
rm(f, meta)

## Subset each rse_gtex object by the genes expressed in brainseq phase 2
gtex_subset <- function(rse, rse_gtex, exon = FALSE, jxn = FALSE) {
    ## Match by feature name
    if(exon | jxn) {
        ov <- findOverlaps(rowRanges(rse), rowRanges(rse_gtex), type = 'equal', ignore.strand = FALSE)
        m <- subjectHits(ov)
        if(jxn) {
            jxn_df <- DataFrame(lapply(seq_len(ncol(rse_gtex)),
                function(x) Rle(0, nrow(rse)) ))
            colnames(jxn_df) <- colnames(rse_gtex)
            rle_df <- DataFrame(apply(assays(rse_gtex)$counts[m, ], 2, Rle))
            jxn_df[queryHits(ov), ] <- rle_df
            rownames(jxn_df) <- rownames(rse)
            colnames(jxn_df) <- colnames(rle_df)
            rse_gtex <- SummarizedExperiment(assays = list(counts = jxn_df),
                rowRanges = rowRanges(rse), colData = colData(rse_gtex))
        }
    } else {
        m <- match(rownames(rse), rownames(rse_gtex))
    }
    print('Missing features')
    print(table(is.na(m)))

    ## Subset
    if (!jxn) {
        rse_gtex <- rse_gtex[m, ]
    }

    ## Set as factor
    colData(rse_gtex)$Region <- factor('DLPFC', levels = c('DLPFC', 'HIPPO'))
    colData(rse_gtex)$Race <- NA

    # table(grepl('_female', rse_gtex_gene$bigwig_file), grepl('_male', rse_gtex_gene$bigwig_file))
    colData(rse_gtex)$Sex <- relevel(factor(ifelse(grepl('_female', rse_gtex$bigwig_file), 'F', 'M')), ref = 'F')

    ## Add means (these are just the same numbers, but for consistency with
    ## the rest of the BrainSeq code I'll call them means)
    colData(rse_gtex)$mean_mitoRate <- colData(rse_gtex)$mitoRate
    colData(rse_gtex)$mean_totalAssignedGene <- colData(rse_gtex)$totalAssignedGene
    colData(rse_gtex)$mean_rRNA_rate <- colData(rse_gtex)$rRNA_rate
    colData(rse_gtex)$mean_RIN <- colData(rse_gtex)$RIN

    print('Dimensions of the data')
    print(dim(rse_gtex))
    return(rse_gtex)
}

## Actually perform the subset
rse_gtex_gene <- gtex_subset(rse_gene, rse_gtex_gene)
rse_gtex_exon <- gtex_subset(rse_exon, rse_gtex_exon, exon = TRUE)
## Note here that the GTEx data set is unstranded while the BrainSeq phase 2 data is stranded...
table(strand(rse_jxn))
table(strand(rse_gtex_jxn))
rse_gtex_jxn <- gtex_subset(rse_jxn, rse_gtex_jxn, jxn = TRUE)
table(strand(rse_gtex_jxn))
rse_gtex_tx <- gtex_subset(rse_tx, rse_gtex_tx)

## Save filtered results
save(rse_gtex_gene, file = 'rse_gtex_gene.Rdata')
save(rse_gtex_exon, file = 'rse_gtex_exon.Rdata')
save(rse_gtex_jxn, file = 'rse_gtex_jxn.Rdata')
save(rse_gtex_tx, file = 'rse_gtex_tx.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
