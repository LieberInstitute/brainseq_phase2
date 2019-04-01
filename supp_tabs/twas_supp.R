library('sessioninfo')
library('tibble')
library('data.table')
# export _JAVA_OPTIONS="-Xms40g -Xmx60g"
library('xlsx')

load('../twas/rda/tt_objects.Rdata', verbose = TRUE)

## To add the exon ids
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/browser/rda/exon_name_map.Rdata', verbose = TRUE)
setkey(exon_name_map, libd_bsp2)
tt$exonGencodeID <- exon_name_map[.(tt$ID), gencode]
region_twas_z$exonGencodeID <- exon_name_map[.(region_twas_z$ID), gencode]

## Remove unused stuff
tt$FILE <- gsub('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/', '', tt$FILE)
tt <- tt[, -which(colnames(tt) %in% c('PANEL', 'type', 'BEST.GWAS.P', 'BEST.GWAS.OR', 'BEST.GWAS.SE', 'BEST.GWAS.FDR'))]

## Comparison versus Gusev
load('../twas/rda/gusev_gene.Rdata', verbose = TRUE)
colnames(gusev_gene) <- gsub('psycm', 'pgc2_clozuk', colnames(gusev_gene))

table('FDR' = tt$TWAS.FDR < 0.05, 'Bonf' = tt$TWAS.Bonf < 0.05)
#        Bonf
# FDR      FALSE   TRUE
#   FALSE 117410      0
#   TRUE    8495   1346

## Export
write.xlsx2(as.data.frame(tt[tt$TWAS.FDR < 0.05, ]), file = 'BrainSeqPhaseII_TWAS_TableSxx.xlsx', sheetName='TWAS PGC2+CLOZUK FDR<5% results')
write.xlsx2(region_twas_z, file = 'BrainSeqPhaseII_TWAS_TableSxx.xlsx', append = TRUE, sheetName='PGC2+CLOZUK TWAS Z DLPFCvsHIPPO')
write.xlsx2(gusev_gene, file = 'BrainSeqPhaseII_TWAS_TableSxx.xlsx', append = TRUE, sheetName='Versus Gusev et al 2018')


region_twas_z_clozuk <- region_twas_z
tt_clozuk <- tt

## Add the PGC2 TWAS results too
load('../twas/rda/pgc2_tt_objects.Rdata', verbose = TRUE)
tt$exonGencodeID <- exon_name_map[.(tt$ID), gencode]
region_twas_z$exonGencodeID <- exon_name_map[.(region_twas_z$ID), gencode]
tt$FILE <- gsub('/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/', '', tt$FILE)
tt <- tt[, -which(colnames(tt) %in% c('PANEL', 'type', 'BEST.GWAS.P', 'BEST.GWAS.OR', 'BEST.GWAS.SE', 'BEST.GWAS.FDR'))]

write.xlsx2(as.data.frame(tt[tt$TWAS.FDR < 0.05, ]), file = 'BrainSeqPhaseII_TWAS_TableSxx.xlsx', append = TRUE, sheetName='TWAS PGC2 FDR<5% results')
write.xlsx2(region_twas_z, file = 'BrainSeqPhaseII_TWAS_TableSxx.xlsx', append = TRUE, sheetName='PGC2 TWAS Z DLPFC vs HIPPO')

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
#  date     2019-03-29
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date       lib source
#  assertthat    0.2.1   2019-03-21 [2] CRAN (R 3.5.1)
#  cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.1)
#  colorout    * 1.2-0   2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace    1.4-1   2019-03-18 [2] CRAN (R 3.5.1)
#  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)
#  data.table  * 1.12.0  2019-01-13 [1] CRAN (R 3.5.1)
#  digest        0.6.18  2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr         0.8.0.1 2019-02-15 [1] CRAN (R 3.5.1)
#  fansi         0.4.0   2018-10-05 [1] CRAN (R 3.5.1)
#  ggplot2       3.1.0   2018-10-25 [1] CRAN (R 3.5.1)
#  glue          1.3.1   2019-03-12 [1] CRAN (R 3.5.1)
#  gtable        0.3.0   2019-03-25 [2] CRAN (R 3.5.1)
#  htmltools     0.3.6   2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets   1.3     2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv        1.5.0   2019-03-15 [2] CRAN (R 3.5.1)
#  jsonlite      1.6     2018-12-07 [2] CRAN (R 3.5.1)
#  later         0.8.0   2019-02-11 [2] CRAN (R 3.5.1)
#  lattice       0.20-38 2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval      0.2.2   2019-03-15 [2] CRAN (R 3.5.1)
#  magrittr      1.5     2014-11-22 [1] CRAN (R 3.5.0)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 3.5.1)
#  pillar        1.3.1   2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig     2.0.2   2018-08-16 [1] CRAN (R 3.5.1)
#  plyr          1.8.4   2016-06-08 [2] CRAN (R 3.5.0)
#  png           0.1-7   2013-12-03 [2] CRAN (R 3.5.0)
#  promises      1.0.1   2018-04-13 [2] CRAN (R 3.5.0)
#  purrr         0.3.2   2019-03-15 [2] CRAN (R 3.5.1)
#  R6            2.4.0   2019-02-14 [2] CRAN (R 3.5.1)
#  Rcpp          1.0.0   2018-11-07 [1] CRAN (R 3.5.1)
#  rJava         0.9-10  2018-05-29 [2] CRAN (R 3.5.1)
#  rlang         0.3.1   2019-01-08 [1] CRAN (R 3.5.1)
#  rmote       * 0.3.4   2018-05-02 [1] deltarho (R 3.5.0)
#  scales        1.0.0   2018-08-09 [2] CRAN (R 3.5.1)
#  servr         0.13    2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo * 1.1.1   2018-11-05 [1] CRAN (R 3.5.1)
#  tibble      * 2.0.1   2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect    0.2.5   2018-10-11 [2] CRAN (R 3.5.1)
#  utf8          1.1.4   2018-05-24 [1] CRAN (R 3.5.0)
#  withr         2.1.2   2018-03-15 [2] CRAN (R 3.5.0)
#  xfun          0.5     2019-02-20 [1] CRAN (R 3.5.1)
#  xlsx        * 0.6.1   2018-06-11 [2] CRAN (R 3.5.1)
#  xlsxjars      0.6.1   2014-08-22 [2] CRAN (R 3.5.0)
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
