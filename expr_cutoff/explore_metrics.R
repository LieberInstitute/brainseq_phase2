# Usage:
# module load conda_R/3.4.x
# mkdir -p logs
# Rscript explore_metrics.R  > logs/explore_metrics.txt 2>&1

library('SummarizedExperiment')
library('devtools')
load('rse_gene.Rdata')

dir.create('pdf', showWarnings = FALSE)
pdf('pdf/explore_metrics.pdf')
for(var in c('RIN', 'totalAssignedGene', 'mitoRate', 'overallMapRate', 'ERCCsumLogErr', 'rRNA_rate')) {
    boxplot(mean(colData(rse_gene)[, var]) ~ colData(rse_gene)$Region, ylab = paste('mean', var), xlab = 'Region', col = 'light blue', main = paste('p-value:', signif(t.test(mean(colData(rse_gene)[, var]) ~ colData(rse_gene)$Region)$p.value, 3)))

}
for(var in c('Age', 'snpPC1', 'snpPC2', 'snpPC3', 'snpPC4', 'snpPC5')) {
    boxplot(colData(rse_gene)[, var] ~ colData(rse_gene)$Region, ylab = var, xlab = 'Region', col = 'light blue', main = paste('p-value:', signif(t.test(colData(rse_gene)[, var] ~ colData(rse_gene)$Region)$p.value, 3)))
}

plot(mean(colData(rse_gene)$mitoRate), mean(colData(rse_gene)$rRNA_rate), pch = ifelse(colData(rse_gene)$Region == 'DLPFC', 0, 1))
legend('topright', legend = c('DLPFC', 'HIPPO'), pch = c(0, 1))
plot(mean(colData(rse_gene)$mitoRate), mean(colData(rse_gene)$ERCCsumLogErr), pch = ifelse(colData(rse_gene)$Region == 'DLPFC', 0, 1))
legend('right', legend = c('DLPFC', 'HIPPO'), pch = c(0, 1))
plot(mean(colData(rse_gene)$rRNA_rate), mean(colData(rse_gene)$ERCCsumLogErr), pch = ifelse(colData(rse_gene)$Region == 'DLPFC', 0, 1))
legend('right', legend = c('DLPFC', 'HIPPO'), pch = c(0, 1))

boxplot(mean(colData(rse_gene)$mitoRate) ~ factor(mean(colData(rse_gene)$ERCCsumLogErr) > -200) * colData(rse_gene)$Region, col = 'light blue', ylab = 'mean mitoRate', xlab = 'mean ERCCsumLogErr > -200 by region')
boxplot(mean(colData(rse_gene)$rRNA_rate) ~ factor(mean(colData(rse_gene)$ERCCsumLogErr) > -200) * colData(rse_gene)$Region, col = 'light blue', ylab = 'mean rRNA_rate', xlab = 'mean ERCCsumLogErr > -200 by region')
dev.off()

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
