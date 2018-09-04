###
library('dplyr')
library('SummarizedExperiment')
library('jaffelab')
library('devtools')
library('scales')
library('derfinder')
library('clusterProfiler')

## Load expression data (cleaned and uncleaned) + rse objects for the pheno data
load('rda/expr_and_cleaned.Rdata', verbose = TRUE)
load('rda/rse_and_modQsva.Rdata', verbose = TRUE)
## and correlations
load('rda/indv_corr.Rdata', verbose=TRUE)

## get brainnums
small_pd = as.data.frame(colData(simple_rse$DLPFC$gene)[,c("BrNum", "Race", "Dx")])

## add prs
prs = read.delim("/dcl01/lieber/ajaffe/Brain/Imputation/PRS/prs.pgc2_clumped.postmortem_2017.p0.05.txt",as.is=TRUE)
prs$IID[prs$IID == "Br1060"] = "Br1061"
small_pd$prs_s6 = prs$SCORE[match(small_pd$BrNum, prs$IID)]

## add global corrs
cors = cbind(bind_cols(indv_expr), bind_cols(indv_cleaned))
colnames(cors) = paste0(colnames(cors), "_", rep(c("Obs","Clean"), each=4))
small_pd = cbind(small_pd, cors)

## just cauc
pdCauc = small_pd[small_pd$Race == "CAUC",]

## test
boxplot(prs_s6 ~ Dx, data = pdCauc)
plot(prs_s6 ~ jxnRp10m_Clean, data = pdCauc)

################################################
## association between global corr and PRS? ####
corTestList = apply(pdCauc[,5:ncol(pdCauc)], 2, function(x) cor.test(x, pdCauc$prs_s6))
corDf = data.frame(Pearson = sapply(corTestList, "[[", "estimate"), 
				pvalue = sapply(corTestList, "[[", "p.value"))
rownames(corDf) = ss(rownames(corDf), "\\.")
signif(corDf,3)
               # Pearson pvalue
# geneRpkm_Obs    0.0584  0.535
# exonRpkm_Obs    0.0454  0.630
# jxnRp10m_Obs    0.0202  0.831
# txTpm_Obs       0.0299  0.751
# geneRpkm_Clean  0.1080  0.251
# exonRpkm_Clean  0.1280  0.174
# jxnRp10m_Clean  0.1100  0.242
# txTpm_Clean     0.1150  0.222

########################################################### 
## adjusting for PRS on global corr associations with dx? #
###########################################################

szEffectsCauc = t(apply(pdCauc[,5:ncol(pdCauc)], 2, function(x) {
	summary(lm(x ~ pdCauc$Dx))$coef[2,]
}))
szEffectsCaucAdj = t(apply(pdCauc[,5:ncol(pdCauc)], 2, function(x) {
	summary(lm(x ~ pdCauc$Dx + pdCauc$prs_s6))$coef[2,]
}))
szEffectsCaucInt = t(apply(pdCauc[,5:ncol(pdCauc)], 2, function(x) {
	summary(lm(x ~ pdCauc$Dx *pdCauc$prs_s6))$coef[4,]
}))

##############################
## effect estimates overall  #
##############################
szEffects = t(apply(small_pd[,5:ncol(small_pd)], 2, function(x) {
	summary(lm(x ~ small_pd$Dx))$coef[2,]
}))
raceEffects = t(apply(small_pd[,5:ncol(small_pd)], 2, function(x) {
	summary(lm(x ~ small_pd$Race))$coef[2,]
}))