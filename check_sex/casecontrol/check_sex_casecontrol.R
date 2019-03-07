## Based on https://github.com/LieberInstitute/brainseq_phase2/blob/master/development/sz_effect_devel_overlap.R
library('GenomicRanges')
library('SummarizedExperiment')
library('purrr')
library('sessioninfo')

## For the chr info
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata", verbose = TRUE)


## List the files for each brain region
## either the original ones (DE by Dx) or the
## sex ones (DE by sex)
f_sex <- c(
    'DLPFC' = 'rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda',
    'HIPPO' = 'rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda'
)
stopifnot(all(file.exists(f_sex)))
f_ori <- paste0('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/', f_sex)
names(f_ori) <- names(f_sex)
stopifnot(all(file.exists(f_ori)))

## Create a similar object to:
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/development/sz_effect_devel_overlap.R#L5-L10
load_stats <- function(stats_files) {
    final <- do.call(c, map2(stats_files, names(stats_files), function(f, region) {
        message(paste(Sys.time(), 'loading', f))
        load(f, verbose = TRUE)
        
        stopifnot(identical(rownames(outGene), names(rowRanges(rse_gene))))
        stopifnot(identical(rownames(outExon), names(rowRanges(rse_exon))))
        stopifnot(identical(rownames(outJxn), names(rowRanges(rse_jxn))))
        stopifnot(identical(rownames(outTx), names(rowRanges(rse_tx))))
        
        res <- list(
            'Gene' = cbind(outGene, chr = seqnames(rowRanges(rse_gene))),
            'Exon' = cbind(outExon, chr = seqnames(rowRanges(rse_exon))),
            'Jxn' = cbind(outJxn, chr = seqnames(rowRanges(rse_jxn))),
            'Tx' = cbind(outTx, chr = seqnames(rowRanges(rse_tx)))
        )
        names(res) <- paste0(region, '_', names(res))
        
        ## Add Bonferroni
        res <- map(res, function(x) {
            x$p_bonf = p.adjust(x$P.Value, method = "bonf")
        	return(x)
        })
        return(res)
    }))
    names(final) <- gsub('.*\\.', '', names(final))
    return(final)
}

sz_stats <- load_stats(f_ori)
sex_stats <- load_stats(f_sex)

map_int(sz_stats, nrow)
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#      24652     396583     297181      92732      24652     396583     297181
#   HIPPO_Tx
#      92732

stopifnot(identical(map_int(sz_stats, nrow), map_int(sex_stats, nrow)))

## Line up
sz_stats <- map2(sz_stats, sex_stats, ~ .x[rownames(.y), ])

## Add sex stats
sz_stats <- map2(sz_stats, sex_stats, function(sz, sex) {
    sz$sexReg <- sex$adj.P.Val < 0.05
    sz$sexM <- sex$adj.P.Val < 0.05 & sex$t > 0
    sz$sexF <- sex$adj.P.Val < 0.05 & sex$t < 0
    sz$sext <- sex$t
    return(sz)
})

## Without chrX and Y
sz_stats_nosex <- map(sz_stats, ~ .x[!.x$chr %in% c('chrX', 'chrY'), ])
map_int(sz_stats_nosex, nrow)
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#      23777     384980     287599      90119      23777     384980     287599
#   HIPPO_Tx
#      90119
map_int(sz_stats_nosex, nrow) / map_int(sz_stats, nrow) * 100
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#   96.45059   97.07426   96.77570   97.18220   96.45059   97.07426   96.77570
#   HIPPO_Tx
#   97.18220


map_int(sz_stats, ~ sum(.x$sexReg))
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#        282       1576       1517        340        116       1261       1109
#   HIPPO_Tx
#        356
map_int(sz_stats_nosex, ~ sum(.x$sexReg))
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#        202        537        610        106         45        274        301
#   HIPPO_Tx
#        128

map_int(sz_stats_nosex, ~ sum(.x$sexReg)) / map_int(sz_stats, ~ sum(.x$sexReg)) * 100
# DLPFC_Gene DLPFC_Exon  DLPFC_Jxn   DLPFC_Tx HIPPO_Gene HIPPO_Exon  HIPPO_Jxn
#   71.63121   34.07360   40.21094   31.17647   38.79310   21.72879   27.14157
#   HIPPO_Tx
#   35.95506


foo <- function(stats) {
    table(
        'SCZD DE' = factor(stats$adj.P.Val < 0.05, levels = c('FALSE', 'TRUE')),
        'Sex DE' = factor(stats$sexReg, levels = c('FALSE', 'TRUE'))
    )
}

map(sz_stats, foo)

map_dbl(sz_stats, ~ chisq.test(foo(.x), simulate.p.value = TRUE, B = 1e4)$p.value)

map(sz_stats, ~ chisq.test(foo(.x)))

map(sz_stats, foo)[[7]]



map_dbl(sz_stats_nosex, ~ chisq.test(foo(.x), simulate.p.value = TRUE, B = 1e4)$p.value)
map(sz_stats_nosex, foo)[[7]]
map(sz_stats_nosex, foo)[[1]]

map_dbl(sz_stats, ~ cor(.x$t, .x$sext))
map_dbl(sz_stats_nosex, ~ cor(.x$t, .x$sext))

pdf('pdf/test.pdf')
with(sz_stats[[1]], plot(t, sext))
dev.off()

library('jaffelab')
map_dbl(sz_stats, ~ getOR(foo(.x)))

## find overlaps by features
get_out <- function(stats) {
    o_stats <- map(stats, function(x) {
    	o <- rbind(summary(lm(x$t ~ x$sexM))$coef[2,c(1,4)],
    		summary(lm(x$t ~ x$sexF))$coef[2,c(1,4)])
    	rownames(o) <- c("male", "female")
    	return(o)
    
    })
    
    out = t(sapply(o_stats, function(x) c(x[,1], x[,2])))
    colnames(out) = paste0(colnames(out), "_", rep(c("shift", "pval"),each=2))
    out = as.data.frame(out[c(1,5,2,6,3,7,4,8),])
    
    out$male_bonf <- p.adjust(out$male_pval, 'bonf')
    out$female_bonf <- p.adjust(out$female_pval, 'bonf')
    return(out)
    
}

out <- get_out(sz_stats)
out

out_nosex <- get_out(sz_stats_nosex)
out_nosex


## all plots
make_pdf <- function(stats, out, filename = 'pdf/densityPlots_dxEffects_devRegByDir.pdf') {
    pdf(filename,w=10,h=4, useDingbats = FALSE)
    par(mar=c(5,6,2,2), mfrow = c(1,2), cex.axis=2,cex.lab=2)
    for(i in seq(along=stats)) {
    	x = stats[[i]]
    		plot(density(x$t[x$sexM==0]),
    		col="grey",lwd=2,main=names(stats)[i],
    		xlab="SZ vs Control",xlim=c(-7,7), ylim = c(0, 0.4))
    	lines(density(x$t[x$sexM==1]),
    		col="darkorange",lwd=3)
    	legend("topright", paste0("p=", signif(out$male_pval[i],3)),
    		bty="n", cex=1.5)
		
    	plot(density(x$t[x$sexF==0]),
    		col="grey",lwd=2,main=names(stats)[i],
    		xlab="SZ vs Control",xlim=c(-7,7), 	ylim = c(0, 0.4))
    	lines(density(x$t[x$sexF==1]),
    		col="darkblue",lwd=3)
    	legend("topright", paste0("p=", signif(out$female_pval[i],3)),
    		bty="n", cex=1.5)
		

    }
    dev.off()
}

make_pdf(sz_stats, out, filename = 'pdf/test.pdf')


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
