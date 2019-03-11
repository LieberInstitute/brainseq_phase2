library('SummarizedExperiment')
library('jaffelab')
library('sessioninfo')

## Load expr data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

## Load gene correlation results
load('rda/gene_pinfo.Rdata', verbose = TRUE)

## Read in housekeeping gene symbols
## List comes from RSeQC http://rseqc.sourceforge.net/
## which we used in
## https://www.pnas.org/content/pnas/suppl/2017/06/19/1617384114.DCSupplemental/pnas.1617384114.sapp.pdf
## _A framework for RNA quality correction in differential expression analysis_ 2017 PNAS paper
## by Jaffe et al
housekeeping <- readLines('/users/ajaffe/Lieber/Projects/RnaQualPaper/OLD/housekeeping.txt')

gr <- rowRanges(rse_gene)
m <- match(rownames(gene_pinfo), names(gr))
stopifnot(!any(is.na(m)))

## Ok, no need to use `m` since they match
stopifnot(identical(m, seq_len(nrow(gene_pinfo))))  s

gene_pinfo <- cbind(gene_pinfo, as.data.frame(gr))

length(housekeeping)
# [1] 553
length(unique(housekeeping))
# [1] 549
gene_pinfo$housekeeping <- tolower(gene_pinfo$Symbol) %in% unique(tolower(housekeeping))

table(gene_pinfo$housekeeping)
# FALSE  TRUE
# 24138   514
table(gene_pinfo$housekeeping) / nrow(gene_pinfo) * 100
#     FALSE      TRUE
# 97.914976  2.085024

house_tab <- with(gene_pinfo, table('High corr' = cleaned.fwer < 0.05, 'Housekeeping' = housekeeping))
addmargins(house_tab)
#          Housekeeping
# High corr FALSE  TRUE   Sum
#     FALSE 23216   506 23722
#     TRUE    922     8   930
#     Sum   24138   514 24652
getOR(house_tab)
# [1] 0.3981035
chisq.test(house_tab)
#
#     Pearson's Chi-squared test with Yates' continuity correction
#
# data:  tab
# X-squared = 6.4919, df = 1, p-value = 0.01084
chisq.test(tab)$p.value
# [1] 0.01083679

## Save for later
save(gene_pinfo, house_tab, file = 'rda/gene_pinfo_housekeeping.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
