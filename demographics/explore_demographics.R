# Usage:
# Rscript explore_demographics.R  >explore_demographics_log.txt 2>&1

library('SummarizedExperiment')
library('ggplot2')
library('LIBDpheno')
library('devtools')

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)

## Add all the metadata
m <- match(as.integer(gsub('Br', '', colData(rse_gene)$BrNum)), toxicant[['2018-05-31']]$brnumerical)
table(is.na(m))
any(colnames(colData(rse_gene)) %in% colnames(toxicant[['2018-05-31']]))
colData(rse_gene) <- cbind(colData(rse_gene), toxicant[['2018-05-31']][m, ])

## Split by subsets
pds <- list(
    all = colData(rse_gene),
    ctrl = colData(rse_gene)[colData(rse_gene)$Dx == 'Control', ]
)

pds$ctrl_hippo <- pds$ctrl[pds$ctrl$Region == 'HIPPO', ]
pds$ctrl_dlpfc <- pds$ctrl[pds$ctrl$Region == 'DLPFC', ]

pds$ctrl_adult <- pds$ctrl[pds$ctrl$Age >= 18, ]
pds$ctrl_adult_hippo <- pds$ctrl_adult[pds$ctrl_adult$Region == 'HIPPO', ]
pds$ctrl_adult_dlpfc <- pds$ctrl_adult[pds$ctrl_adult$Region == 'DLPFC', ]

pds$ctrl_fetal <- pds$ctrl[pds$ctrl$Age < 0, ]
pds$ctrl_fetal_hippo <- pds$ctrl_fetal[pds$ctrl_fetal$Region == 'HIPPO', ]
pds$ctrl_fetal_dlpfc <- pds$ctrl_fetal[pds$ctrl_fetal$Region == 'DLPFC', ]

load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata', verbose = TRUE)
pds$sczd_hippo <- colData(rse_gene[, keepIndex])
pds$sczd_hippo_case <- pds$sczd_hippo[pds$sczd_hippo$Dx == 'Schizo', ]
pds$sczd_hippo_ctrl <- pds$sczd_hippo[pds$sczd_hippo$Dx == 'Control', ]

load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_DLPFC.Rdata', verbose = TRUE)
pds$sczd_dlpfc <- colData(rse_gene[, keepIndex])
pds$sczd_dlpfc_case <- pds$sczd_dlpfc[pds$sczd_dlpfc$Dx == 'Schizo', ]
pds$sczd_dlpfc_ctrl <- pds$sczd_dlpfc[pds$sczd_dlpfc$Dx == 'Control', ]

#pd <- pds$all

summpd <- function(pd) {
    data.frame(
        n = nrow(pd),
        sexF = mean(pd$Sex == 'F') * 100,
        raceCauc = mean(pd$Race == 'CAUC') * 100,
        age = mean(pd$Age),
        RIN = mean(mean(pd$RIN)),
        PMI = mean(pd$pmi),
        smoke = mean(pd$smoking, na.rm = TRUE) * 100,
        smoke_NA_n = sum(is.na(pd$smoking)),
        death_accident = mean(pd$manner_of_death == 'Accident', na.rm = TRUE) * 100,
        death_homicide = mean(pd$manner_of_death == 'Homicide', na.rm = TRUE) * 100,
        death_natural = mean(pd$manner_of_death == 'Natural', na.rm = TRUE) * 100,
        death_suicide = mean(pd$manner_of_death == 'Suicide', na.rm = TRUE) * 100,
        death_NA_n = sum(is.na(pd$manner_of_death)),
        mapping_perc = mean(mean(pd$overallMapRate)) * 100,
        total_mapped = mean(mean(pd$totalMapped / 1e6)),
        chrM_mapping_perc = mean(mean(pd$mitoRate)) * 100,
        assigned_gene_perc = mean(mean(pd$totalAssignedGene)) * 100,
        rRNA_perc = mean(mean(pd$rRNA_rate)) * 100,
        ageonset = mean(pd$age_onset_schizo, na.rm = TRUE),
        ageonset_NA_n = sum(is.na(pd$age_onset_schizo)),
        antipsychotics = mean(pd$antipsychotics, na.rm = TRUE) * 100,
        antipsychotics_NA_n = sum(is.na(pd$antipsychotics))
    )    
}

#pd1 <- pds$sczd_hippo_case
#pd2 <- pds$sczd_hippo_ctrl

comppd <- function(pd1, pd2) {    
    
    data.frame(
        n = nrow(pd1) + nrow(pd2),
        sexF = t.test(pd1$Sex == 'F', pd2$Sex == 'F')$p.value,
        raceCauc = t.test(pd1$Race == 'CAUC', pd2$Race == 'CAUC')$p.value,
        age = t.test(pd1$Age, pd2$Age)$p.value,
        RIN = t.test(mean(pd1$RIN), mean(pd2$RIN))$p.value,
        PMI = t.test(pd1$pmi, pd2$pmi)$p.value,
        smoke = tryCatch(t.test(pd1$smoking, pd2$smoking)$p.value, error = function(e) return(NA)),
        smoke_NA_n = sum(is.na(pd1$smoking)) + sum(is.na(pd2$smoking)),
        death_accident = tryCatch(t.test(pd1$manner_of_death == 'Accident', pd2$manner_of_death == 'Accident')$p.value, error = function(e) return(NA)),
        death_homicide = tryCatch(t.test(pd1$manner_of_death == 'Homicide', pd2$manner_of_death == 'Homicide')$p.value, error = function(e) return(NA)),
        death_natural = tryCatch(t.test(pd1$manner_of_death == 'Natural', pd2$manner_of_death == 'Natural')$p.value, error = function(e) return(NA)),
        death_suicide = tryCatch(t.test(pd1$manner_of_death == 'Suicide', pd2$manner_of_death == 'Suicide')$p.value, error = function(e) return(NA)),
        death_NA_n = sum(is.na(pd1$manner_of_death)) + sum(is.na(pd2$manner_of_death)),
        mapping_perc = t.test(mean(pd1$overallMapRate), mean(pd2$overallMapRate))$p.value,
        total_mapped = t.test(mean(pd1$totalMapped / 1e6), mean(pd2$totalMapped / 1e6))$p.value,
        chrM_mapping_perc = t.test(mean(pd1$mitoRate), mean(pd2$mitoRate))$p.value,
        assigned_gene_perc = t.test(mean(pd1$totalAssignedGene), mean(pd2$totalAssignedGene))$p.value,
        rRNA_perc = t.test(mean(pd1$rRNA_rate), mean(pd2$rRNA_rate))$p.value,
        ageonset = NA,
        ageonset_NA_n = sum(is.na(pd1$age_onset_schizo)) + sum(is.na(pd2$age_onset_schizo)),
        antipsychotics = NA,
        antipsychotics_NA_n = sum(is.na(pd1$antipsychotics)) + sum(is.na(pd2$antipsychotics))
    )
    
}


summ <- do.call(rbind, lapply(pds, summpd))

pds1 <- pds[grep('sczd.*case|ctrl.*hippo', names(pds))]
pds2 <- pds[grep('sczd.*ctrl|ctrl.*dlpfc', names(pds))]
comps <- mapply(comppd, pds1, pds2)

colnames(comps) <- paste0(names(pds1), '_vs_', names(pds2))
comps <- t(comps)


supptab <- rbind(summ, comps)

options(width = 300)
supptab
t(supptab)

supptab2 <- supptab[c(1:4, 17, 5:7, 18, 8:10, 19, 11:13, 20, 14:16, 21), ]

save(pds, summ, comps, supptab, supptab2, file = 'supptab.Rdata')

write.csv(supptab2, file = 'supptab.csv')

pdf('supptab_plots.pdf', useDingbats = FALSE)
with(pds$sczd_hippo, boxplot(pmi ~ Dx, ylab = 'PMI (in hours)', xlab = 'SCZD status', main = 'HIPPO'))
with(pds$sczd_dlpfc, boxplot(pmi ~ Dx, ylab = 'PMI (in hours)', xlab = 'SCZD status', main = 'DLPFC'))
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
