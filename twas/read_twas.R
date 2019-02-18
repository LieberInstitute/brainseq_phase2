library('readr')
library('purrr')
library('dplyr')
library('stringr')
library('testthat')
library('sessioninfo')

## Main options to work through
regions <- c('DLPFC', 'HIPPO')
types <- c('psycm', 'pgc2')
features <- c('gene', 'exon', 'jxn', 'tx')
file_types <- c('all', 'included', 'dropped')

## Function for locating the files
locate_files <- function(ftype, path, type) {
    
    ## Determine the file pattern to use
    fpatt <- case_when(
        ftype == 'all' ~ '',
        ftype == 'included' ~ '\\.analysis\\.joint_included',
        ftype == 'dropped' ~ '\\.analysis\\.joint_dropped'
    )
    patt <- paste0(
        type,
        '\\.[[:digit:]]*',
        fpatt, 
        '\\.dat$'
    )

    ## Find the files
    dat_files <- dir(
        path = path,
        pattern = patt,
        full.names = TRUE
    )
    
    ## Extract the chromosome
    names(dat_files) <- str_extract(
        basename(dat_files),
        '[:digit:]+'
    )
    
    ## Done
    return(dat_files)
}

## Check that the file locator is working
test_files <- locate_files(
    'all',
    '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene/psycm',
    'psycm'
)

test_that('File locator', {
    expect_equal(
        basename(test_files['11']),
        'psycm.11.dat'
    )
    expect_equivalent(
        test_files['1'],
        '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/gene/psycm/psycm.1.dat'
    )
})


## Now read in the data for all file types
## Note that each file type has different numbers of columns
## when why I'm keeping them as different elements of the
## 'twas' list
twas <- map(file_types, function(ftype) {
    
    ## data.frame with all the combinations of arguments
    arg_grid <- expand.grid(
        region = regions,
        type = types,
        feature = features,
        stringsAsFactors = FALSE
    )
    
    
    pmap_dfr(arg_grid, function(region, type, feature) {
        
        ## Construct the path to the files
        path <- file.path(
            '/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas',
            region,
            feature,
            type
        )
        
        ## Now locate the files
        dat_files <- locate_files(
            ftype = ftype,
            path = path,
            type = type
        )
        
        ## Next read the files and add the chromosome info
        result <- map2_dfr(
            dat_files,
            names(dat_files),
            function(f, chr) {
                res <- suppressMessages(read_tsv(f))
                res$chr <- chr
                return(res)
            }
        )
        
        ## Next add the region, feature and type information
        result$region <- region
        result$feature <- feature
        result$type <- type
        
        ## Done
        return(result)
    })
})
names(twas) <- file_types

## Explore the resulting dimensions
map_dfr(twas, dim)
# # A tibble: 2 x 3
#      all included dropped
#    <int>    <int>   <int>
# 1 408043     1572  254181
# 2     24       12      12

## Save the data for later use
dir.create('rda', showWarnings = FALSE)
save(twas, file = 'rda/twas.Rdata')

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
#  date     2019-02-18
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
#  digest        0.6.18  2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr       * 0.7.8   2018-11-10 [1] CRAN (R 3.5.1)
#  fansi         0.4.0   2018-10-05 [1] CRAN (R 3.5.1)
#  ggplot2       3.1.0   2018-10-25 [1] CRAN (R 3.5.1)
#  glue          1.3.0   2018-07-17 [1] CRAN (R 3.5.1)
#  gtable        0.2.0   2016-02-26 [2] CRAN (R 3.5.0)
#  hms           0.4.2   2018-03-10 [2] CRAN (R 3.5.0)
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
#  purrr       * 0.2.5   2018-05-29 [2] CRAN (R 3.5.0)
#  R6            2.3.0   2018-10-04 [2] CRAN (R 3.5.1)
#  Rcpp          1.0.0   2018-11-07 [1] CRAN (R 3.5.1)
#  readr       * 1.3.1   2018-12-21 [1] CRAN (R 3.5.1)
#  rlang         0.3.1   2019-01-08 [1] CRAN (R 3.5.1)
#  rmote       * 0.3.4   2018-05-02 [1] deltarho (R 3.5.0)
#  scales        1.0.0   2018-08-09 [2] CRAN (R 3.5.1)
#  servr         0.11    2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo * 1.1.1   2018-11-05 [1] CRAN (R 3.5.1)
#  stringi       1.2.4   2018-07-20 [2] CRAN (R 3.5.1)
#  stringr     * 1.3.1   2018-05-10 [1] CRAN (R 3.5.0)
#  testthat    * 2.0.1   2018-10-13 [1] CRAN (R 3.5.1)
#  tibble        2.0.1   2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect    0.2.5   2018-10-11 [2] CRAN (R 3.5.1)
#  utf8          1.1.4   2018-05-24 [1] CRAN (R 3.5.0)
#  withr         2.1.2   2018-03-15 [2] CRAN (R 3.5.0)
#  xfun          0.4     2018-10-23 [1] CRAN (R 3.5.1)
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
