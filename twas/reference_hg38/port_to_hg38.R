## Load plink before starting R
module load plink/1.90b6.6

## and copy the LDREF files
cp -R /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/filter_data/LDREF LDREF_hg38

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

# ref_bim <- ldref_bim[[1]]; bim_file <- their_bims[1]; chr <- '1'
mapply(function(ref_bim, bim_file, chr) {
    our <- snpMap[.(ref_bim$chr, ref_bim$basepair)]
    m <- !is.na(our$pos_hg38)
    
    ## Write the list of snps that do show up in our data
    filt_snps <- paste0('1000G.EUR.', chr, '.snps.txt')
    fwrite(ref_bim$snp[m],
        file = filter_snps,
        sep = '\t', col.names = FALSE
    )   
    
    newbfile <- paste0('LDREF_hg38/1000G.EUR.', chr)
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
    m_final <- 
    final_bim$basepair <- our$pos_hg38[m_final]
    
    message(paste(Sys.time(), 'Write new filtered bim file for chr', chr))
    fwrite(
        final_bim,
        file = newbfile_bim,
        sep = '\t', col.names = FALSE
    )

    return(file.exists(newbfile_bim))
    
}, ldref_bim, their_bims, gsub('1000G.EUR.|.bim', '', names(their_bims)))



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
