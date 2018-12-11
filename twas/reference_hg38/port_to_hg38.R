## Load plink before starting R
module load plink/1.90b6.6

## and copy the LDREF files
cp -R /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LDREF LDREF_hg38

R

## Now run R code
library('data.table')
# library('SummarizedExperiment')
library('sessioninfo')

## based on ../filter_data/filter_snps.R
their_bims <- dir('LDREF_hg38', '.*bim$', full.names = TRUE)
names(their_bims) <- dir('LDREF_hg38', '.*bim$')

ldref_bim <- lapply(their_bims, function(input_bim) {
    message(paste(Sys.time(), 'reading file', input_bim))
    original <- paste0(input_bim, '.original')
    if(file.exists(original)) {
        ## start from the original file if it's around
        system(paste('mv', original, input_bim))
    }
    res <- fread(input_bim,
        col.names = c('chr', 'snp', 'position', 'basepair', 'allele1', 'allele2'),
        colClasses = c('character', 'character', 'numeric', 'integer', 'character', 'character')
    )
    system(paste('mv', input_bim, original))
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

# ref_bim <- ldref_bim[[1]]; bim_file <- their_bims[1]; chr <- '1'
mapply(function(ref_bim, bim_file, chr) {
    our <- snpMap[.(ref_bim$chr, ref_bim$basepair)]
    m <- !is.na(our$pos_hg38)
    new <- ref_bim[m, ]
    new$basepair <- our$basepair[m]
    
    message(paste(Sys.time(), 'Write new filtered bim file for chr', chr))
    fwrite(new,
        file = bim_file,
        sep = '\t', col.names = FALSE
    )
    return(file.exists(bim_file))
    
}, ldref_bim, their_bims, gsub('1000G.EUR.|.bim', '', names(their_bims)))


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
