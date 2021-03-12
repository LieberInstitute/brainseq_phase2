library("here")
library("GenomicRanges")
library("sessioninfo")


files <- dir(
    here(
        "eQTL_full",
        "subset_hancock"
    ),
    "mergedEqtl_output_",
    full.names = TRUE
)

eqtls <- lapply(files, function(f) {
    load(f, verbose = TRUE)
    return(allEqtl)
})

eQTLs_Hancock <- do.call(rbind, eqtls)
dim(eQTLs_Hancock)
# [1] 18507    12
table(eQTLs_Hancock$Region)
# DLPFC       HIPPO Interaction
#  6834        6834        4839


## To fix exon IDs
## From https://github.com/LieberInstitute/brainseq_phase2/blob/master/browser/extract_data.R
load(here("browser", "rda", "exon_name_map.Rdata"), verbose = TRUE)

colnames(eQTLs_Hancock)[colnames(eQTLs_Hancock) == "snps"] <- "snp"
colnames(eQTLs_Hancock)[colnames(eQTLs_Hancock) == "gene"] <- "feature_id"
colnames(eQTLs_Hancock)[colnames(eQTLs_Hancock) == "Region"] <- "eQTL_model"
eQTLs_Hancock$Type <- tolower(eQTLs_Hancock$Type)

m_e <- match(eQTLs_Hancock$feature_id, exon_name_map$libd_bsp2)
table(is.na(m_e), eQTLs_Hancock$Type)
#        exon  gene   jxn    tx
# FALSE 11066     0     0     0
# TRUE      0   408  4738  2295

eQTLs_Hancock$feature_id[!is.na(m_e)] <- exon_name_map$gencode[m_e[!is.na(m_e)]]

fix_csv <- function(df) {
    for (i in seq_len(ncol(df))) {
        if (any(grepl(",", df[, i]))) {
            message(paste(Sys.time(), "fixing column", colnames(df)[i]))
            df[, i] <- gsub(",", ";", df[, i])
        }
    }
    return(df)
}

x <- fix_csv(as.data.frame(eQTLs_Hancock))
write.csv(x, here(
    "eQTL_full",
    "subset_hancock",
    "eQTLs_Hancock.csv"
    )
)


save(eQTLs_Hancock, file = here(
    "eQTL_full",
    "subset_hancock",
    "eQTLs_Hancock.Rdata"
    )
)


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.0.4 RC (2021-02-08 r79975)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-03-12
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package          * version  date       lib source
#  assertthat         0.2.1    2019-03-21 [2] CRAN (R 4.0.3)
#  BiocGenerics     * 0.36.0   2020-10-27 [1] Bioconductor
#  bitops             1.0-6    2013-08-17 [2] CRAN (R 4.0.3)
#  cli                2.2.0    2020-11-20 [1] CRAN (R 4.0.3)
#  colorout           1.2-2    2020-05-09 [1] Github (jalvesaq/colorout@726d681)
#  colorspace         2.0-0    2020-11-11 [2] CRAN (R 4.0.3)
#  crayon             1.4.1    2021-02-08 [1] CRAN (R 4.0.4)
#  digest             0.6.27   2020-10-24 [1] CRAN (R 4.0.3)
#  dplyr              1.0.2    2020-08-18 [1] CRAN (R 4.0.3)
#  ellipsis           0.3.1    2020-05-15 [1] CRAN (R 4.0.3)
#  fansi              0.4.1    2020-01-08 [1] CRAN (R 4.0.0)
#  generics           0.1.0    2020-10-31 [1] CRAN (R 4.0.3)
#  GenomeInfoDb     * 1.26.2   2020-12-08 [1] Bioconductor
#  GenomeInfoDbData   1.2.4    2020-11-30 [2] Bioconductor
#  GenomicRanges    * 1.42.0   2020-10-27 [1] Bioconductor
#  ggplot2            3.3.3    2020-12-30 [1] CRAN (R 4.0.3)
#  glue               1.4.2    2020-08-27 [1] CRAN (R 4.0.3)
#  gtable             0.3.0    2019-03-25 [2] CRAN (R 4.0.3)
#  here             * 1.0.1    2020-12-13 [1] CRAN (R 4.0.3)
#  htmltools          0.5.0    2020-06-16 [1] CRAN (R 4.0.3)
#  htmlwidgets        1.5.3    2020-12-10 [1] CRAN (R 4.0.3)
#  httpuv             1.5.4    2020-06-06 [1] CRAN (R 4.0.3)
#  IRanges          * 2.24.1   2020-12-12 [1] Bioconductor
#  jsonlite           1.7.2    2020-12-09 [2] CRAN (R 4.0.3)
#  later              1.1.0.1  2020-06-05 [1] CRAN (R 4.0.3)
#  lattice            0.20-41  2020-04-02 [3] CRAN (R 4.0.4)
#  lifecycle          0.2.0    2020-03-06 [1] CRAN (R 4.0.0)
#  magrittr           2.0.1    2020-11-17 [1] CRAN (R 4.0.3)
#  munsell            0.5.0    2018-06-12 [2] CRAN (R 4.0.3)
#  pillar             1.4.7    2020-11-20 [1] CRAN (R 4.0.3)
#  pkgconfig          2.0.3    2019-09-22 [1] CRAN (R 4.0.0)
#  png                0.1-7    2013-12-03 [2] CRAN (R 4.0.3)
#  promises           1.1.1    2020-06-09 [1] CRAN (R 4.0.3)
#  purrr              0.3.4    2020-04-17 [1] CRAN (R 4.0.0)
#  R6                 2.5.0    2020-10-28 [2] CRAN (R 4.0.3)
#  Rcpp               1.0.5    2020-07-06 [1] CRAN (R 4.0.3)
#  RCurl              1.98-1.2 2020-04-18 [2] CRAN (R 4.0.3)
#  rlang              0.4.10   2020-12-30 [1] CRAN (R 4.0.3)
#  rmote              0.3.4    2020-05-09 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot          2.0.2    2020-11-15 [2] CRAN (R 4.0.3)
#  S4Vectors        * 0.28.1   2020-12-09 [1] Bioconductor
#  scales             1.1.1    2020-05-11 [2] CRAN (R 4.0.3)
#  servr              0.21     2020-12-14 [1] CRAN (R 4.0.3)
#  sessioninfo      * 1.1.1    2018-11-05 [1] CRAN (R 4.0.3)
#  tibble             3.0.4    2020-10-12 [1] CRAN (R 4.0.3)
#  tidyselect         1.1.0    2020-05-11 [2] CRAN (R 4.0.3)
#  vctrs              0.3.6    2020-12-17 [1] CRAN (R 4.0.3)
#  withr              2.3.0    2020-09-22 [1] CRAN (R 4.0.3)
#  xfun               0.20     2021-01-06 [1] CRAN (R 4.0.3)
#  XVector            0.30.0   2020-10-27 [1] Bioconductor
#  zlibbioc           1.36.0   2020-10-27 [1] Bioconductor
#
# [1] /users/lcollado/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library