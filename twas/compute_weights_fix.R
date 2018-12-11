library('SummarizedExperiment')
library('jaffelab')
library('data.table')
library('sessioninfo')
library('getopt')


## Fix previous files
for(reg in c('DLPFC', 'HIPPO')) {
    for (feat in c('gene', 'exon', 'jxn', 'tx')) {
        # for(feat in 'gene')
        print(paste(reg, '-', feat))
        
        opt <- list(region = reg, feature = feat, cores = 1, 'pgconly' = TRUE)
        # opt <- list(region = reg, feature = feat, cores = 1, 'pgconly' = FALSE)
        
        rse_file <- file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, '_pgconly', '')), 'working_rse.Rdata')
        if(!file.exists(rse_file)) {
            rse <- load_rse(opt$feature, opt$region)
            message(paste(Sys.time(), 'saving the rse file for later at', rse_file))
            save(rse, file = rse_file)
        } else {
            message(paste(Sys.time(), 'loading previous rse file', rse_file))
            load(rse_file, verbose = TRUE)
        }
        
        bim_file <- paste0(
            '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_LDfiltered_',
            opt$region
        )
        bim <- fread(
            paste0(bim_file, '.bim'),
            col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2')
        )
        bim_gr <- GRanges(
            paste0('chr', bim$chr),
            IRanges(bim$basepair, width = 1)
        )
        mcols(bim_gr) <- bim[, - c('chr', 'basepair')]
        
        
        ## Based on http://gusevlab.org/projects/fusion/#computing-your-own-functional-weights
        ## they restrict to 500kb on each side
        rse_window <- resize(rowRanges(rse), width(rowRanges(rse)) + 500000 * 2, fix = 'center')


        if(opt$pgconly) {
            load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2/scz2_anneal_gr_hg19.Rdata', verbose = TRUE)
    
            ## Which of the best snps did we actually observe?
            table(countOverlaps(scz2_anneal_gr_hg19, bim_gr))
            #  0  1
            # 77 31
            scz2_anneal_gr_hg19 <- scz2_anneal_gr_hg19[countOverlaps(scz2_anneal_gr_hg19, bim_gr) > 0]
    
            ## Keep just the feature windows that include the PCG SNPs that we observed
            keep_feat <- unique(subjectHits(findOverlaps(scz2_anneal_gr_hg19, rse_window)))
            rse <- rse[keep_feat, ]
            rse_window <- rse_window[keep_feat, ]
            stopifnot(nrow(rse) == length(rse_window))
        }

        ## Number of SNPs per feature window
        # table(countOverlaps(rse_window, bim_gr))

        ## Keep only those feature windows with some SNPs nearby
        keep_feat <- which(countOverlaps(rse_window, bim_gr) > 0)
        rse <- rse[keep_feat, ]
        rse_window <- rse_window[keep_feat, ]
        stopifnot(nrow(rse) == length(rse_window))

        print('Final RSE feature dimensions:')
        print(dim(rse))

        rse_file <- file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, '_pgconly', '')), 'subsetted_rse.Rdata')
        message(paste(Sys.time(), 'saving the subsetted rse file for later at', rse_file))
        save(rse, file = rse_file)

        ## Subset to features
        setwd(file.path(opt$region, paste0(opt$feature, ifelse(opt$pgconly, '_pgconly', ''))))
        dir.create('snp_files', showWarnings = FALSE)
        dir.create('bim_files', showWarnings = FALSE)
        dir.create('tmp_files', showWarnings = FALSE)
        dir.create('out_files', showWarnings = FALSE)
        
        output_status <- sapply(seq_len(nrow(rse)), function(i) {
            out_file <- file.path(getwd(), 'out_files', paste0(opt$feature, '_', i))
            return(file.exists(paste0(out_file, '.wgt.RDat')))
        })

        message('*******************************************************************************')
        message(paste(Sys.time(), 'summary output status (TRUE means that there is a file)'))
        print(table(output_status))

        message(paste(Sys.time(), 'saving the output_status.Rdata file'))
        names(output_status) <- paste0(opt$feature, '_', seq_len(nrow(rse)))
        save(output_status, file = 'output_status.Rdata')

        message(paste(Sys.time(), 'creating the wgt profile files'))
        rdat_files <- dir('out_files', '.wgt.RDat', full.names = TRUE)
        stopifnot(length(rdat_files) == sum(output_status))

        wglist <- paste0(opt$region, '_', opt$feature, '.list')
        write.table(rdat_files, file = wglist, row.names = FALSE, col.names = FALSE, quote = FALSE)
        system(paste0('Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/utils/FUSION.profile_wgt.R ', wglist, ' > ', opt$region, '_', opt$feature, '.profile.err 2>&1'))
        ## Rename the wglist_summary.txt file to keep the naming convention consistent
        system(paste0('mv wglist_summary.txt ', opt$region, '_', opt$feature, '.profile'))

        message(paste(Sys.time(), 'creating the .pos file'))
        if(opt$feature %in% c('gene', 'exon')) {
            var <- 'gencodeID'
            # vars <- 'Symbol'
        } else if (opt$feature == 'jxn') {
            var <- 'gencodeGeneID'
            # vars <- 'Symbol'
        } else {
            var <- 'gene_id'
            # vars <- 'gene_name'
        }

        pos_info <- data.frame(
            'WGT' = rdat_files,
            'ID' = mcols(rowRanges(rse))[, var][output_status],
            'CHR' = gsub('chr', '', seqnames(rowRanges(rse)[output_status])),
            'P0' = start(rowRanges(rse[output_status])), 
            'P1' = end(rowRanges(rse[output_status])),
            stringsAsFactors = FALSE
        )
        print(sapply(pos_info, function(x) sum(is.na(x))))
        pos_info$ID[pos_info$ID == ''] <- NA
        write.table(pos_info, file = paste0(opt$region, '_', opt$feature, '.pos'), row.names = FALSE, col.names = TRUE, quote = FALSE)
        save(pos_info, file = 'pos_info.Rdata')
        
        setwd('../..')
        print(getwd())
    }
}
