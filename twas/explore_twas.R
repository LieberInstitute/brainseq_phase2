library('readr')
library('purrr')
library('dplyr')
library('stringr')
library('sessioninfo')


regions <- c('DLPFC', 'HIPPO')
types <- c('psycm', 'pgc2')
features <- c('gene', 'exon', 'jxn', 'tx')
file_types <- c('all', 'included', 'dropped')


map()

region <- 'DLPFC'
type <- 'psycm'
feature <- 'gene'



path <- file.path(
    '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas',
    region,
    feature,
    type
)

ftype <- 'all'
ftype <- file_types[3]

fpatt <- case_when(
    ftype == 'all' ~ '',
    ftype == 'included' ~ '\\.analysis\\.joint_included',
    ftype == 'dropped' ~ '\\.analysis\\.joint_dropped'
)



patt <- paste0(type, '\\.[[:digit:]]*', fpatt, '\\.dat$')

patt

dat_files <- dir(path = path, pattern = patt, full.names = TRUE)
## Extract the chromosome
names(dat_files) <- str_extract(basename(dat_files), '[:digit:]+')


testfoo <- function(f, chr, region, feature, type) {
    res <- read_tsv(f)
    res$chr <- chr
    res$region <- region
    res$feature <- feature
    res$type <- type
    return(res)
}

x <- map(file_types, function(ftype) {

    fpatt <- case_when(
        ftype == 'all' ~ '',
        ftype == 'included' ~ '\\.analysis\\.joint_included',
        ftype == 'dropped' ~ '\\.analysis\\.joint_dropped'
    )



    patt <- paste0(type, '\\.[[:digit:]]*', fpatt, '\\.dat$')

    print(patt)

    dat_files <- dir(path = path, pattern = patt, full.names = TRUE)
    ## Extract the chromosome
    names(dat_files) <- str_extract(basename(dat_files), '[:digit:]+')
    
    map2_dfr(dat_files[1:2], names(dat_files)[1:2], testfoo, region = region, feature = feature, type = type)
})



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
