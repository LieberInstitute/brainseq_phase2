## Load plink before starting R
qrsh -l mem_free=100G,h_vmem=100G,h_fsize=100G
module load plink/1.90b6.6
R

## Now run R code
library('data.table')
library('sessioninfo')

dir.create('LDREF_hg38')

## based on ../filter_data/filter_snps.R
their_bims <- dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LDREF', '.*bim$', full.names = TRUE)
names(their_bims) <- dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LDREF', '.*bim$')

ldref_bim <- lapply(their_bims, function(input_bim) {
    message(paste(Sys.time(), 'reading file', input_bim))
    res <- fread(input_bim,
        col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2'),
        colClasses = c('character', 'character', 'numeric', 'integer', 'character', 'character')
    )
    setkey(res, 'chr', 'basepair')
    return(res)
})

## get the snpMap that contains the hg38 positions
message(paste(Sys.time(), 'loading BSP2 genotype data'))
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda', verbose = TRUE)

colnames(snpMap) <- tolower(colnames(snpMap))
colnames(snpMap)[colnames(snpMap) == 'pos'] <- 'basepair'
snpMap <- data.table(snpMap)
setkey(snpMap, 'chr', 'basepair')

# ref_bim <- ldref_bim[['1000G.EUR.22.bim']]; bim_file <- their_bims['1000G.EUR.22.bim']; chr <- '22'
newbims_exist <- mapply(function(ref_bim, bim_file, chr) {
    message('*****************************************************')
    message(paste(Sys.time(), 'start processing chr', chr))
    
    our <- snpMap[.(ref_bim$chr, ref_bim$basepair)]
    m <- !is.na(our$pos_hg38)
    
    ## Write the list of snps that do show up in our data
    filt_snps <- paste0('1000G.EUR.', chr, '.snps.txt')
    fwrite(as.data.frame(ref_bim$snp[m]),
        file = filt_snps,
        sep = '\t', col.names = FALSE
    )   
    
    newbfile <- file.path('LDREF_hg38', paste0('1000G.EUR.', chr))
    newbfile_bim <- paste0(newbfile, '.bim')
    
    message(paste(Sys.time(), 'subsetting the original LDREF files with plink'))
    system(paste("plink --bfile", gsub('.bim', '', bim_file), '--extract', filt_snps,
         "--make-bed --out", newbfile, " --memory 10000"))
        
    ## Read the new file
    message(paste(Sys.time(), 'read the new bim file'))
    final_bim <- fread(newbfile_bim,
        col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2'),
        colClasses = c('character', 'character', 'numeric', 'integer', 'character', 'character')
    )
    
    ## Complete matching code later
    m_final <- match(with(final_bim, paste0(chr, '-', basepair)), with(our, paste0(chr, '-', basepair)))
    stopifnot(!any(is.na(m_final)))
    final_bim$basepair <- our$pos_hg38[m_final]
    
    message(paste(Sys.time(), 'Write new filtered bim file for chr', chr))
    fwrite(
        final_bim,
        file = newbfile_bim,
        sep = '\t', col.names = FALSE
    )

    return(file.exists(newbfile_bim))
    
}, ldref_bim, their_bims, gsub('1000G.EUR.|.bim', '', names(their_bims)))

## Last bit of the log:

# *****************************************************
# 2018-12-14 15:31:09 start processing chr 9
# 2018-12-14 15:31:11 subsetting the original LDREF files with plink
# PLINK v1.90b6.6 64-bit (10 Oct 2018)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to LDREF_hg38/1000G.EUR.9.log.
# Options in effect:
#   --bfile /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LDREF/1000G.EUR.9
#   --extract 1000G.EUR.9.snps.txt
#   --make-bed
#   --memory 10000
#   --out LDREF_hg38/1000G.EUR.9
#
# 517005 MB RAM detected; reserving 10000 MB for main workspace.
# 55464 variants loaded from .bim file.
# 489 people (0 males, 0 females, 489 ambiguous) loaded from .fam.
# Ambiguous sex IDs written to LDREF_hg38/1000G.EUR.9.nosex .
# --extract: 47883 variants remaining.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 489 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# 47883 variants and 489 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to LDREF_hg38/1000G.EUR.9.bed + LDREF_hg38/1000G.EUR.9.bim +
# LDREF_hg38/1000G.EUR.9.fam ... done.
# 2018-12-14 15:31:12 read the new bim file
# 2018-12-14 15:31:12 Write new filtered bim file for chr 9

table(newbims_exist)
# TRUE
#   22

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.5.1 Patched (2018-10-29 r75535)
#  os       Red Hat Enterprise Linux Server release 6.9 (Santiago)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2018-12-14
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date       lib source
#  assertthat    0.2.0   2017-04-11 [2] CRAN (R 3.5.0)
#  bindr         0.1.1   2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp      0.2.2   2018-03-29 [1] CRAN (R 3.5.0)
#  cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.1)
#  colorout    * 1.2-0   2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace    1.3-2   2016-12-14 [2] CRAN (R 3.5.0)
#  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)
#  data.table  * 1.11.8  2018-09-30 [1] CRAN (R 3.5.1)
#  digest        0.6.18  2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr         0.7.8   2018-11-10 [1] CRAN (R 3.5.1)
#  ggplot2       3.1.0   2018-10-25 [1] CRAN (R 3.5.1)
#  glue          1.3.0   2018-07-17 [1] CRAN (R 3.5.1)
#  gtable        0.2.0   2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools     0.3.6   2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets   1.3     2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv        1.4.5   2018-07-19 [2] CRAN (R 3.5.1)
#  later         0.7.5   2018-09-18 [2] CRAN (R 3.5.1)
#  lattice       0.20-38 2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval      0.2.1   2017-10-29 [2] CRAN (R 3.5.0)
#  magrittr      1.5     2014-11-22 [1] CRAN (R 3.5.0)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 3.5.0)
#  pillar        1.3.0   2018-07-14 [1] CRAN (R 3.5.1)
#  pkgconfig     2.0.2   2018-08-16 [1] CRAN (R 3.5.1)
#  plyr          1.8.4   2016-06-08 [2] CRAN (R 3.5.0)
#  png           0.1-7   2013-12-03 [2] CRAN (R 3.5.0)
#  promises      1.0.1   2018-04-13 [2] CRAN (R 3.5.0)
#  purrr         0.2.5   2018-05-29 [2] CRAN (R 3.5.0)
#  R6            2.3.0   2018-10-04 [2] CRAN (R 3.5.1)
#  Rcpp          1.0.0   2018-11-07 [1] CRAN (R 3.5.1)
#  rlang         0.3.0.1 2018-10-25 [1] CRAN (R 3.5.1)
#  rmote       * 0.3.4   2018-05-02 [1] deltarho (R 3.5.0)
#  scales        1.0.0   2018-08-09 [2] CRAN (R 3.5.1)
#  servr         0.11    2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo * 1.1.1   2018-11-05 [1] CRAN (R 3.5.1)
#  tibble        1.4.2   2018-01-22 [1] CRAN (R 3.5.0)
#  tidyselect    0.2.5   2018-10-11 [2] CRAN (R 3.5.1)
#  withr         2.1.2   2018-03-15 [2] CRAN (R 3.5.0)
#  xfun          0.4     2018-10-23 [1] CRAN (R 3.5.1)
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
