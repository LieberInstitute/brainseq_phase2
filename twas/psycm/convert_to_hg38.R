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
# 3326255  989758

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
