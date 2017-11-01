library('SummarizedExperiment')
library('recount')
library('jaffelab')
library('devtools')

files <- dir('../count_data', pattern = 'hg38_rse', full.names = TRUE)
names(files) <- dir('../count_data', pattern = 'hg38_rse')

types <- gsub('_.*', '', gsub('.*_hg38_rse', '', files))
regions <- toupper(gsub('_.*', '', names(files)))

## Load the raw data and calculate RPKMs & RP10M when necessary
all <- mapply(function(f, type, region) {
    load(f)
    if(type == 'Gene') {
        assays(rse_gene)$rpkm <- recount::getRPKM(rse_gene, 'Length')
        rowRanges(rse_gene)$meanExprs <- NA
        colData(rse_gene)$Region <- region
        return(rse_gene)
    } else if (type == 'Exon') {
        rowRanges(rse_exon)$meanExprs <- NA
        assays(rse_exon)$rpkm <- recount::getRPKM(rse_exon, 'Length')
        colData(rse_exon)$Region <- region
        return(rse_exon)
    } else if (type == 'Jxn') {
        ## Try getRPKM based on https://github.com/LieberInstitute/brainseq_phase2/blob/54c73b2b4cd65af93a254ff8f38eed6a8d5c362a/caseControl_analysis_hippo.R#L42
        rowRanges(rse_jxn)$Length <- 100
        #assays(rse_jxn)$rp10m <- recount::getRPKM(rse_jxn, 'Length')
        rowRanges(rse_jxn)$meanExprs <- NA
        colData(rse_jxn)$Region <- region
        return(rse_jxn)
    } else if (type == 'Tx') {
        rowRanges(rse_tx)$meanExprs <- NA
        colData(rse_tx)$Region <- region
        return(rse_tx)
    }
}, files, types, regions)


rse_merge <- function(rses) {
    cols <- sapply(rses, function(x) colnames(colData(x)))
    common <- intersect(cols[[1]], cols[[2]])
    do.call(cbind, lapply(rses, function(r) {
        m <- match(common, colnames(colData(r)))
        colData(r) <- colData(r)[, m[!is.na(m)]]
        return(r)
    }))
}

## Combine across regions
rse_gene <- rse_merge(all[types == 'Gene'])
rse_exon <- rse_merge(all[types == 'Exon'])
rse_tx <- rse_merge(all[types == 'Tx'])

## More complicated case since not all junctions are in both sets
## First merge the counts
jxn_raw <- lapply(all[types == 'Jxn'], function(x) { assays(x)$counts })
all_jxn <- unique(unlist(sapply(jxn_raw, rownames)))
jxn <- matrix(0, nrow = length(all_jxn), ncol = sum(sapply(jxn_raw, ncol)))

m1 <- match(all_jxn, rownames(jxn_raw[[1]]))
jxn[!is.na(m1), seq_len(ncol(jxn_raw[[1]]))] <- jxn_raw[[1]][m1[!is.na(m1)], ]
m2 <- match(all_jxn, rownames(jxn_raw[[2]]))
jxn[!is.na(m2), seq_len(ncol(jxn_raw[[2]])) + ncol(jxn_raw[[1]])] <- jxn_raw[[2]][m2[!is.na(m2)], ]

## Merge jx range info
jxn_gr <- rep(rowRanges(all[[which(types == 'Jxn')[1]]])[1], nrow(jxn))
jxn_gr[!is.na(m1)] <- rowRanges(all[[which(types == 'Jxn')[1]]])[m1[!is.na(m1)]]
jxn_gr[!is.na(m2)] <- rowRanges(all[[which(types == 'Jxn')[2]]])[m2[!is.na(m2)]]
stopifnot(sum(jxn_gr == rowRanges(all[[which(types == 'Jxn')[1]]])[1]) == 1)

## Now merge the columns
cols <- sapply(all[types == 'Jxn'], function(x) colnames(colData(x)))
common <- intersect(cols[[1]], cols[[2]])
jxn_col <- do.call(rbind, lapply(all[types == 'Jxn'], function(r) {
    m <- match(common, colnames(colData(r)))
    colData(r)[, m[!is.na(m)]]
}))

rse_jxn <- SummarizedExperiment(assays = list(counts = jxn), rowRanges = jxn_gr, colData = jxn_col)
# fix junction row names
rownames(rse_jxn) <- paste0(seqnames(rse_jxn), ":", start(rse_jxn), "-",
    end(rse_jxn), "(", strand(rse_jxn), ")")
assays(rse_jxn)$rp10m <- recount::getRPKM(rse_jxn, 'Length')


exprs <- list(
    'Gene' = assays(rse_gene)$rpkm,
    'Exon' = assays(rse_exon)$rpkm,
    'Jxn' = assays(rse_jxn)$rp10m,
    'Tx' = assay(rse_tx)
)


## Identify potential cutoffs
seed <- 20171026
seeds <- seed + 0:3
names(seeds) <- names(exprs)

cutoffs <- sapply(names(exprs), function(type) {
    
    message(type)
    pdf(paste0('suggested_expr_cutoffs_', tolower(type), '.pdf'), width = 12)
    cuts <- expression_cutoff(exprs[[type]], seed = seeds[type])
    message(paste(cuts, collapse = ' '))
    cut <- max(cuts)
    dev.off()
    
    return(cut)

})
# Gene
# 2017-10-27 16:18:28 the suggested expression cutoff is 0.21
# 0.25 0.17
# Exon
# 2017-10-27 16:20:14 the suggested expression cutoff is 0.26
# 0.3 0.22
# Jxn
# 2017-10-27 16:21:28 the suggested expression cutoff is 0.4
# 0.33 0.46
# Tx
# 2017-10-27 16:22:42 the suggested expression cutoff is 0.32
# 0.38 0.25

cutoffs
# Gene Exon  Jxn   Tx
# 0.25 0.30 0.46 0.38

means <- lapply(exprs, rowMeans)

## Add the mean expressions, whether it passes the expression cutoff
## and save the data
rowRanges(rse_gene)$meanExprs <- means[['Gene']]
rowRanges(rse_gene)$passExprsCut <- means[['Gene']] > cutoffs['Gene']
save(rse_gene, file = 'rse_gene_unfiltered.Rdata')
rse_gene <- rse_gene[rowRanges(rse_gene)$passExprsCut]
save(rse_gene, file = 'rse_gene.Rdata')

rowRanges(rse_exon)$meanExprs <- means[['Exon']]
rowRanges(rse_exon)$passExprsCut <- means[['Exon']] > cutoffs['Exon']
save(rse_exon, file = 'rse_exon_unfiltered.Rdata')
rse_exon <- rse_exon[rowRanges(rse_exon)$passExprsCut]
save(rse_exon, file = 'rse_exon.Rdata')

rowRanges(rse_jxn)$meanExprs <- means[['Jxn']]
rowRanges(rse_jxn)$passExprsCut <- means[['Jxn']] > cutoffs['Jxn']
save(rse_jxn, file = 'rse_jxn_unfiltered.Rdata')
rse_jxn <- rse_jxn[rowRanges(rse_jxn)$passExprsCut]
save(rse_jxn, file = 'rse_jxn.Rdata')

rowRanges(rse_tx)$meanExprs <- means[['Tx']]
rowRanges(rse_tx)$passExprsCut <- means[['Tx']] > cutoffs['Tx']
save(rse_tx, file = 'rse_tx_unfiltered.Rdata')
rse_tx <- rse_tx[rowRanges(rse_tx)$passExprsCut]
save(rse_tx, file = 'rse_tx.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
