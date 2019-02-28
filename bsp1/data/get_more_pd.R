library('SummarizedExperiment')
library('sessioninfo')

load('bsp1_gene.Rdata', verbose = TRUE)
write.table(colData(bsp1_gene)$BrNum, file = 'brnums.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

## Next use LIBD pheno to create the csv file with all the extra info
# scp e:/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/brnums.txt .
# scp /Users/lcollado/Downloads/LIBDpheno_selection_2019-02-26\ 14_30_48.csv e:/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/
pd_all <- read.csv('LIBDpheno_selection_2019-02-26 14_30_48.csv')

vars_of_interest <- c('brnumerical', 'brnum', 'agedeath', 'sex', 'race', 'primarydx', 'pmi', 'clean_cod', 'pmi_confidence_level', 'bmi_calculated', 'autopsy_date', 'date_frozen', 'age_onset_schizo', 'age_onset_mdd')
pd <- pd_all[, vars_of_interest]
summary(pd)
write.csv(pd, file = 'pd_for_hanna.csv', row.names = FALSE, quote = FALSE)
# scp e:/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/pd_for_hanna.csv .

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
