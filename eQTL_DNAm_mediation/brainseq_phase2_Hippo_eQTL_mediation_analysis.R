#qsub -l bluejay,mf=80G,h_vmem=80G,h_fsize=200G,h_stack=256M -cwd -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH --no-save brainseq_phase2_Hippo_eQTL_mediation_analysis.R

library(SummarizedExperiment)
## Load the map of potential mediations to test
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_DNAm_mediation/PGC2_meQTL_eQTL_Hippo_mergedSet_p01_union.rda')
mergedSetNomSig$se_methylation = mergedSetNomSig$beta_methylation / mergedSetNomSig$statistic_methylation
mergedSetNomSig$se_expression = mergedSetNomSig$beta_expression / mergedSetNomSig$statistic_expression


## Load in the cut methylation data
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_DNAm_mediation/PGC2_Hippo_mediation_bVals.rda')
keepInd = which(pdMeth$Age>13 & pdMeth$Brain.Region =="Hippo")
pdMeth = pdMeth[keepInd, ]
meth_mediation = meth_mediation[,keepInd]
## Load methylation PCs
load('/dcl01/lieber/ajaffe/Steve/Hippo_meQTL/rdas/pcs_methPCs_regions_filtered_over13_Hippo_only.rda')


## Load in the cut expression data
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_DNAm_mediation/PGC2_Hippo_mediation_expression.rda')

## keep adult samples & correct region
keepInd = which(colData(rse_gene_mediation)$Age > 13 & colData(rse_gene_mediation)$Region == "HIPPO")
rse_gene_mediation = rse_gene_mediation[,keepInd]
rse_exon_mediation = rse_exon_mediation[,keepInd]
rse_jxn_mediation = rse_jxn_mediation[,keepInd]
rse_tx_mediation = rse_tx_mediation[,keepInd]

## extract pd and rpkms
pdExprs = colData(rse_gene_mediation)
geneRpkm = assays(rse_gene_mediation)$rpkm
exonRpkm = assays(rse_exon_mediation)$rpkm
jxnRp10m = assays(rse_jxn_mediation)$rp10m
txTpm = assays(rse_tx_mediation)$tpm

## Subset to samples with both methylation and expression data
shared_people = intersect( pdMeth$BrNum, pdExprs$BrNum )
bVals = meth_mediation[, pdMeth[match(shared_people, pdMeth$BrNum),'Chip']  ]
methPCs = methPCs[pdMeth[match(shared_people, pdMeth$BrNum),'Chip'], ]

geneRpkm = geneRpkm[, pdExprs[match(shared_people, pdExprs$BrNum),'RNum'] ]
exonRpkm = exonRpkm[, pdExprs[match(shared_people, pdExprs$BrNum),'RNum'] ]
jxnRp10m = jxnRp10m[, pdExprs[match(shared_people, pdExprs$BrNum),'RNum'] ]
txTpm = txTpm[, pdExprs[match(shared_people, pdExprs$BrNum),'RNum'] ]

exprs4feature = log2(rbind(geneRpkm,exonRpkm,jxnRp10m,txTpm)+1)

## load SNP data
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_DNAm_mediation/PGC2_Hippo_mediation_snps.rda')
snp = snp_mediation[,shared_people]
mds = mds[shared_people,]

## load expression PCs
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/rdas/pcs_hippo_4features_filtered_over13.rda")
genePCs = genePCs[ pdExprs[match(shared_people, pdExprs$BrNum),'RNum'], ]
exonPCs = exonPCs[ pdExprs[match(shared_people, pdExprs$BrNum),'RNum'], ]
jxnPCs = jxnPCs[ pdExprs[match(shared_people, pdExprs$BrNum),'RNum'], ]
txPCs = txPCs[ pdExprs[match(shared_people, pdExprs$BrNum),'RNum'], ]

### 
pd = as.data.frame(pdExprs[ match(shared_people, pdExprs$BrNum), c('Dx','Sex') ])
pd$Dx = factor(pd$Dx, levels=c("Control","Schizo") )
pd$Sex = factor(pd$Sex,levels=c("M","F") )

rownames(pd) = shared_people

mod = cbind(pd[,c('Dx', 'Sex')], as.matrix(mds[,1:5]) )
mod = as.data.frame(model.matrix(~., mod))
####### Test mediation 
getExprsPCs=function(type) switch(type, "Gene"=genePCs, "Exon"=exonPCs, "Jxn"=jxnPCs, "Tx"=txPCs )

#a = for (row_i in 1:nrow(mergedSetNomSig)) 
getMediationStats = function(row_i) {

## Get correct expression PCs
type = mergedSetNomSig[row_i, 'Type' ]
exprsPCs = getExprsPCs(type)

## subset to features of interest
snp_i = mergedSetNomSig[row_i, 'SNP' ]
cpg_i = mergedSetNomSig[row_i, 'cpg' ]
feature_i = mergedSetNomSig[row_i, 'gene']

keepSNP = !is.na(snp[snp_i,])
##### Calculate a: effect of SNP on methylation
mod_a = cbind(mod, methPCs, snp[snp_i,]  )
colnames(mod_a)[ncol(mod_a)] <- c(snp_i)

fastResA=RcppArmadillo::fastLm(X=mod_a[keepSNP,], y=bVals[cpg_i,keepSNP] )
fastResA=data.frame(coef=names(fastResA$coefficients),Estimate=fastResA$coefficients, StdErr=fastResA$stderr)
#fastResA$SNP = snp_i
#fastResA$CpG = cpg_i
fastResA = fastResA[fastResA$coef==snp_i ,]


##### Calculate b: effect of methylation on expression
mod_b = cbind(mod, exprsPCs, snp[snp_i,], bVals[cpg_i,]  )
colnames(mod_b)[(ncol(mod_b)-1):ncol(mod_b)] <- c(snp_i,cpg_i)

fastResB=RcppArmadillo::fastLm(X=mod_b[keepSNP,], y=exprs4feature[feature_i,keepSNP])
fastResB=data.frame(coef=names(fastResB$coefficients),Estimate=fastResB$coefficients, StdErr=fastResB$stderr)
#fastResB$SNP = snp_i
#fastResB$CpG = cpg_i
#fastResB$Feature = feature_i
fastResB = fastResB[fastResB$coef==cpg_i,]
##
res = data.frame(SNP=snp_i,
		CpG=cpg_i,
		Feature=feature_i, 
		a.Estimate_SNP_x_CpG = fastResA$Estimate, 
		a.SE_SNP_x_CpG = fastResA$StdErr,
		b.Estimate_CpG_x_Exprs = fastResB$Estimate,
		b.SE_CpG_x_Exprs = fastResB$StdErr)

if ((row_i %% 500)==0) {print(paste0(as.character( (row_i / nrow(mergedSetNomSig)) *100), "%; ", "row_i=", as.character(row_i)) )}
return(res)
 }

mediationStats = data.table::rbindlist(lapply(1:nrow(mergedSetNomSig), getMediationStats))
# mediationStats = do.call("rbind",lapply(1:nrow(mergedSetNomSig),getMediationStats))
# mediationStats = do.call("rbind",lapply(sample(x=1:nrow(mergedSetNomSig),1000),getMediationStats)) 
save(mediationStats, file='/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_DNAm_mediation/JointModel_Mediation_fullStats.rda')

##
#mediation_spread = data.frame(SNP.Effect = mediationStats$Estimate[mediationStats$SNP==mediationStats$coef], 
#							  SNP.StdErr = mediationStats$StdErr[mediationStats$SNP==mediationStats$coef], 
#							  CpG.Effect = mediationStats$Estimate[mediationStats$CpG==mediationStats$coef] , 
#							  CpG.StdErr = mediationStats$StdErr[mediationStats$CpG==mediationStats$coef]  )

mediation_analysis_results = cbind(mergedSetNomSig, mediationStats)

a=mediation_analysis_results$`a.Estimate_SNP_x_CpG`
sa=mediation_analysis_results$`a.SE_SNP_x_CpG`
b=mediation_analysis_results$`b.Estimate_CpG_x_Exprs`
sb=mediation_analysis_results$`b.SE_CpG_x_Exprs`

## Sobel Test for mediation
mediation_analysis_results$Sobel_Mediation_Z = a*b/sqrt(b^2*sa^2+a^2*sb^2)
mediation_analysis_results$Sobel_Mediation_P = pnorm(-abs(mediation_analysis_results$Sobel_Mediation_Z))*2
mediation_analysis_results$Sobel_Mediation_FDR = p.adjust(mediation_analysis_results$Sobel_Mediation_P, method = "fdr")

## Aroian Test for mediation
mediation_analysis_results$Aroian_Mediation_Z = a*b/sqrt(b^2*sa^2 + a^2*sb^2 + sa^2*sb^2)
mediation_analysis_results$Aroian_Mediation_P = pnorm(-abs(mediation_analysis_results$Aroian_Mediation_Z))*2
mediation_analysis_results$Aroian_Mediation_FDR = p.adjust(mediation_analysis_results$Aroian_Mediation_P, method = "fdr")

table(mediation_analysis_results$Aroian_Mediation_FDR<0.05,mediation_analysis_results$Sobel_Mediation_FDR<0.05 )
pdf('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_DNAm_mediation/plots/eQTL_methylation_mediation_pvalue_distribution.pdf')
hist(mediation_analysis_results$Sobel_Mediation_P)
hist(mediation_analysis_results$Sobel_Mediation_FDR)
hist(mediation_analysis_results$Aroian_Mediation_P)
hist(mediation_analysis_results$Aroian_Mediation_FDR)
dev.off()

save(mediation_analysis_results, file='/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_DNAm_mediation/JointModel_testMediation_Sobel_and_Aroian.rda')

#mod_i = model.matrix(~., mod_i)
#res = summary( lm( exprs4feature[feature_i,] ~ mod_i+0 ) )
#res=as.data.frame(res$coefficients)
#res$SNP = snp_i
#res$CpG = cpg_i
#res$Feature = feature_i
#res$Coefficient = gsub("mod_i","",rownames(res))
#res$Coefficient = gsub("`","",res$Coefficient)
#
#rownames(res) = NULL
#full_res[[i]] = res 
#cat('.')

