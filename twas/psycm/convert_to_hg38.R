library('data.table')
library('sessioninfo')

## get the snpMap that contains the hg38 positions
message(paste(Sys.time(), 'loading BSP2 genotype data'))
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda', verbose = TRUE)

colnames(snpMap) <- tolower(colnames(snpMap))
colnames(snpMap)[colnames(snpMap) == 'pos'] <- 'basepair'
snpMap <- data.table(snpMap)
setkey(snpMap, 'chr', 'basepair')

## Read the data from http://walters.psycm.cf.ac.uk/

psycm <- fread('clozuk_pgc2.meta.reformatted.sumstats')

## Read the original file that has the chromosome and base pair information
psycm_ori <- fread('clozuk_pgc2.meta.sumstats.txt',
    col.names = c('snp', 'freq.a1', 'chr', 'basepair', 'a1', 'a2', 'or', 'se', 'p')
)
psycm_ori$chr <- as.character(psycm_ori$chr)
setkey(psycm_ori, 'snp')
psycm_ori_sub <- psycm_ori[.(psycm$SNP)]
stopifnot(nrow(psycm_ori_sub) == nrow(psycm))
stopifnot(all(!is.na(psycm_ori_sub$chr)))
stopifnot(identical(psycm$SNP, psycm_ori_sub$snp))

## After passing the checks, assign this info to our table
psycm$chr <- psycm_ori_sub$chr
psycm$basepair <- psycm_ori_sub$basepair

## Change 23 to X; Y is not there
as.character(1:24) %in% unique(psycm$chr)
psycm$chr[psycm$chr == '23'] <- 'X'
stopifnot(all(c(as.character(1:22), 'X') %in% unique(psycm$chr)))


ids <- with(psycm, paste(chr, basepair, sep = ':'))
## Ok, they are all unique
identical(length(unique(ids)), length(ids))


ids_snpmap <- with(snpMap, paste(chr, basepair, sep = ':'))
length(unique(ids_snpmap))
length(ids_snpmap)

table(snpMap$type)
#
# Deletion Insertion       SNV
#   371958    319431   6332471

## Looks like chr + pasepair is the best way to match
table(psycm$SNP %in% snpMap$snp)
#   FALSE    TRUE
# 3145509 3354965
table(ids %in% ids_snpmap)
#   FALSE    TRUE
# 2183881 4316593

table(ids %in% ids_snpmap[snpMap$type == 'SNV'])
#   FALSE    TRUE
# 2184290 4316184

table(is.na(snpMap$pos_hg38))
#
#   FALSE    TRUE
# 7023286     574


snpMap_sub <- subset(snpMap, !is.na(pos_hg38) & type == 'SNV')
ids_snpmap_sub <- with(snpMap_sub, paste(chr, basepair, sep = ':'))
table(psycm$SNP %in% snpMap_sub$snp)
#   FALSE    TRUE
# 3145654 3354820
table(ids %in% ids_snpmap_sub)
#   FALSE    TRUE
# 2184461 4316013

m <- match(ids, ids_snpmap_sub)
table(!is.na(m))
#   FALSE    TRUE
# 2184461 4316013

m_na <- is.na(m)
m2 <- m[!m_na]

psycm_hg38 <- psycm[!m_na, ]

## Check before assigning
stopifnot(all(psycm_hg38$basepair == snpMap_sub$basepair[m2]))
stopifnot(all(psycm_hg38$chr == snpMap_sub$chr[m2]))

psycm_hg38$basepairhg19 <- psycm_hg38$basepair
psycm_hg38$originalSNP <- psycm_hg38$SNP

psycm_hg38$basepair <- snpMap_sub$pos_hg38[m2]
psycm_hg38$SNP <- snpMap_sub$snp[m2]

fwrite(psycm_hg38, file = 'clozuk_pgc2.meta.reformatted.sumstats_hg38_ourname', sep = '\t')


## Check if these SNPs are in the LDref
their_bims_hg38_ourname <- dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38', '.*bim$', full.names = TRUE)
names(their_bims_hg38_ourname) <- dir('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38', '.*bim$')

ldref_bim_hg38_ourname <- do.call(rbind, lapply(their_bims_hg38_ourname, function(input_bim) {
    message(paste(Sys.time(), 'reading file', input_bim))
    res <- fread(input_bim,
        col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2'),
        colClasses = c('character', 'character', 'numeric', 'integer', 'character', 'character')
    )
    #setkey(res, 'chr', 'basepair')
    return(res)
}))

ids_converted <- with(psycm_hg38, paste(chr, basepair, sep = ':'))
ids_ldref <- with(ldref_bim_hg38_ourname, paste(chr, basepair, sep = ':'))
table(ids_converted %in% ids_ldref)
#   FALSE    TRUE
# 3319149  996864

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
#  date     2019-01-28
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date       lib source
#  assertthat    0.2.0   2017-04-11 [2] CRAN (R 3.5.0)
#  bindr         0.1.1   2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp      0.2.2   2018-03-29 [1] CRAN (R 3.5.0)
#  cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.1)
#  colorout    * 1.2-0   2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace    1.4-0   2019-01-13 [2] CRAN (R 3.5.1)
#  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)
#  data.table  * 1.12.0  2019-01-13 [1] CRAN (R 3.5.1)
#  digest        0.6.18  2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr         0.7.8   2018-11-10 [1] CRAN (R 3.5.1)
#  ggplot2       3.1.0   2018-10-25 [1] CRAN (R 3.5.1)
#  glue          1.3.0   2018-07-17 [1] CRAN (R 3.5.1)
#  gtable        0.2.0   2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools     0.3.6   2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets   1.3     2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv        1.4.5.1 2018-12-18 [2] CRAN (R 3.5.1)
#  later         0.7.5   2018-09-18 [2] CRAN (R 3.5.1)
#  lattice       0.20-38 2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval      0.2.1   2017-10-29 [2] CRAN (R 3.5.0)
#  magrittr      1.5     2014-11-22 [1] CRAN (R 3.5.0)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 3.5.0)
#  pillar        1.3.1   2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig     2.0.2   2018-08-16 [1] CRAN (R 3.5.1)
#  plyr          1.8.4   2016-06-08 [2] CRAN (R 3.5.0)
#  png           0.1-7   2013-12-03 [2] CRAN (R 3.5.0)
#  promises      1.0.1   2018-04-13 [2] CRAN (R 3.5.0)
#  purrr         0.2.5   2018-05-29 [2] CRAN (R 3.5.0)
#  R6            2.3.0   2018-10-04 [2] CRAN (R 3.5.1)
#  Rcpp          1.0.0   2018-11-07 [1] CRAN (R 3.5.1)
#  rlang         0.3.1   2019-01-08 [1] CRAN (R 3.5.1)
#  rmote       * 0.3.4   2018-05-02 [1] deltarho (R 3.5.0)
#  scales        1.0.0   2018-08-09 [2] CRAN (R 3.5.1)
#  servr         0.11    2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo * 1.1.1   2018-11-05 [1] CRAN (R 3.5.1)
#  tibble        2.0.1   2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect    0.2.5   2018-10-11 [2] CRAN (R 3.5.1)
#  withr         2.1.2   2018-03-15 [2] CRAN (R 3.5.0)
#  xfun          0.4     2018-10-23 [1] CRAN (R 3.5.1)
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
