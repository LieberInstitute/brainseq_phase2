# Usage:
# qrsh -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
# module load conda_R/3.5
# Rscript combine_gtex.R  > logs/combine_gtex.txt 2>&1

library('SummarizedExperiment')
library('devtools')
library('recount')
library('jaffelab')

## Locate and load all rse files
rse_files <- c(dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex',
    pattern = 'rse_gtex_[[:alpha:]]*.Rdata', full.names = TRUE),
    dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_dlpfc',
        pattern = 'rse_gtex_[[:alpha:]]*.Rdata', full.names = TRUE))
stopifnot(length(rse_files) == 8)


reg <- rep(c('hippo', 'dlpfc'), each = 4)
feat <- gsub('.*rse_gtex_', '', gsub('.Rdata', '', rse_files))

names(rse_files) <- sets <- paste0(reg, '_', feat)


load_gtex <- function(feature, rse_file) {
    message(paste(Sys.time(), 'loading', rse_file))
    load(rse_file, verbose = TRUE)
    if(feature == 'gene') {
        return(rse_gtex_gene)
    } else if (feature == 'exon') {
        return(rse_gtex_exon)
    } else if (feature == 'jxn') {
        return(rse_gtex_jxn)
    } else if (feature == 'tx') {
        return(rse_gtex_tx)   
    }
}

xx <- lapply(unique(feat), function(feature) {
    rses <- lapply(unique(reg), function(region) {
        load_gtex(feature, rse_files[which(reg == region & feat == feature)])
    })
    
    if(feature != 'tx') {
        new_mean <- (rowData(rses[[1]])$meanExprs * ncol(rses[[1]]) + rowData(rses[[2]])$meanExprs * ncol(rses[[2]])) / (ncol(rses[[1]]) + ncol(rses[[2]]))
        rowData(rses[[1]])$meanExprs <- rowData(rses[[2]])$meanExprs <- new_mean
    }
    
    ## Commmon columns
    common <- intersect(colnames(colData(rses[[1]])), colnames(colData(rses[[2]])) )
    colData(rses[[1]]) <- colData(rses[[1]])[, common]
    colData(rses[[2]]) <- colData(rses[[2]])[, common]
    
    message(paste(Sys.time(), 'combining rses'))
    new_rse <- do.call(cbind, rses)
    
    
    message(paste(Sys.time(), 'saving new rse file'))
    if(feature == 'gene') {
        rse_gtex_gene <- new_rse
        assays(rse_gtex_gene)$rpkm <- recount::getRPKM(rse_gtex_gene, 'Length')
        save(rse_gtex_gene, file = 'rse_gtex_gene.Rdata')
    } else if (feature == 'exon') {
        rse_gtex_exon <- new_rse
        assays(rse_gtex_exon)$rpkm <- recount::getRPKM(rse_gtex_exon, 'Length')
        save(rse_gtex_exon, file = 'rse_gtex_exon.Rdata')
    } else if (feature == 'jxn') {
        rse_gtex_jxn <- new_rse
        rowRanges(rse_gtex_jxn)$Length <- 100
        jxn <- assays(rse_gtex_jxn)$counts
        rownames(jxn) <- seq_len(nrow(jxn))
        jxn <- as.matrix(as.data.frame(jxn))
        
        jxn_gr <- rowRanges(rse_gtex_jxn)
        names(jxn_gr) <- seq_len(nrow(jxn))
        jxn_col <- colData(rse_gtex_jxn)
        rse_gtex_jxn <- SummarizedExperiment(assays = list(counts = jxn),
            rowRanges = jxn_gr, colData = jxn_col)
        ## Fix rownames
        rownames(jxn) <- paste0(seqnames(rse_gtex_jxn), ":",
            start(rse_gtex_jxn), "-", end(rse_gtex_jxn), "(",
            strand(rse_gtex_jxn), ")")
        assays(rse_gtex_jxn)$rp10m <- recount::getRPKM(rse_gtex_jxn, 'Length')
        save(rse_gtex_jxn, file = 'rse_gtex_jxn.Rdata')
    } else if (feature == 'tx') {
        rse_gtex_tx <- new_rse 
        save(rse_gtex_tx, file = 'rse_gtex_tx.Rdata')
    }
    
    return(NULL)
})


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
