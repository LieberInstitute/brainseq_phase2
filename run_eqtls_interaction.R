####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)


######## DLPFC ####
load("count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseTx_merged_n449.rda")
rse_tx_dlpfc = rse_tx
load("count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseJxn_merged_n449.rda")
rse_jxn_dlpfc = rse_jxn
load("count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseExon_merged_n449.rda")
rse_exon_dlpfc = rse_exon
load("count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseGene_merged_n449.rda")
rse_gene_dlpfc = rse_gene

# sum totalMapped IntegerLists (so getRPKM works later)
colData(rse_tx_dlpfc)$totalMapped =
	colData(rse_jxn_dlpfc)$totalMapped =
	colData(rse_exon_dlpfc)$totalMapped = 
	colData(rse_gene_dlpfc)$totalMapped  = sapply(colData(rse_gene_dlpfc)$totalMapped, sum)


######## Hippo ####
load("count_data/hippo_brainseq_phase2_hg38_rseTx_merged_n442.rda")
rse_tx_hippo = rse_tx
load("count_data/hippo_brainseq_phase2_hg38_rseJxn_merged_n442.rda")
rse_jxn_hippo = rse_jxn
load("count_data/hippo_brainseq_phase2_hg38_rseExon_merged_n442.rda")
rse_exon_hippo = rse_exon
load("count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n442.rda")
rse_gene_hippo = rse_gene

## sum totalMapped IntegerLists (so getRPKM works later)
colData(rse_tx_hippo)$totalMapped =
	colData(rse_jxn_hippo)$totalMapped =
	colData(rse_exon_hippo)$totalMapped = 
	colData(rse_gene_hippo)$totalMapped  = sapply(colData(rse_gene_hippo)$totalMapped, sum)

print("....data loaded....")

###################### 
# merge phenotype ####
pd_dlpfc = colData(rse_gene_dlpfc)
pd_hippo = colData(rse_gene_hippo)
cn = intersect(colnames(pd_dlpfc), colnames(pd_hippo))
pd = rbind(pd_dlpfc[,cn], pd_hippo[,cn])


###################### 
# merge expression ###

# genes
temp1 = getRPKM(rse_gene_dlpfc)
temp2 = getRPKM(rse_gene_hippo)
geneRpkm = cbind(temp1,temp2)

# exons
temp1 = getRPKM(rse_exon_dlpfc)
temp2 = getRPKM(rse_exon_hippo)
exonRpkm = cbind(temp1,temp2)

# jxns
temp1 = getRPM(rse_jxn_dlpfc)
temp2 = getRPM(rse_jxn_hippo)
rn = intersect(rownames(rse_jxn_dlpfc), rownames(rse_jxn_hippo))
temp1 = temp1[rn,]
temp2 = temp2[rn,]
jxnRp80m = cbind(temp1,temp2)

rse_jxn_dlpfc = rse_jxn_dlpfc[rn,]
rse_jxn_hippo = rse_jxn_hippo[rn,]
mcols(rse_jxn_dlpfc) = mcols(rse_jxn_hippo) = mcols(rse_jxn_hippo)[,-16]
colData(rse_jxn_dlpfc) = pd_dlpfc[,cn]
colData(rse_jxn_hippo) = pd_hippo[,cn]
rse_jxn = cbind(rse_jxn_dlpfc,rse_jxn_hippo)

# tx
temp1 = assays(rse_tx_dlpfc)$tpm
temp2 = assays(rse_tx_hippo)$tpm
txTpm = cbind(temp1,temp2)

rm(temp1,temp2, rse_gene_hippo,rse_exon_hippo,rse_jxn_hippo,rse_tx_hippo,
				rse_gene_dlpfc,rse_exon_dlpfc,rse_jxn_dlpfc,rse_tx_dlpfc,
				pd_dlpfc, pd_hippo)

print("....data merged....")

###################### 
## load SNP data
load("genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n546.rda")


### make mds and snp dimensions equal to N
###(repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum


## drop SNPs not mapping to hg38
keepIndex = which(!is.na(snpMap$chr_hg38))
snpMap = snpMap[keepIndex,]
snp = snp[keepIndex,]


######################
# statistical model ##
######################
pd$Dx = factor(pd$Dx,
	levels = c("Control", "Schizo"))

mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]) + Region,
	data = pd)
colnames(mod)[4:8] = colnames(mds)[1:5]

######################
# create SNP objects #
######################

theSnps = SlicedData$new(as.matrix(snp))
theSnps$ResliceCombined(sliceSize = 50000)

snpspos = snpMap[,c("SNP","chr_hg38","pos_hg38")]
colnames(snpspos) = c("name","chr","pos")

######################
# create expression ##
######################

## gene
geneIndex = rowMeans(geneRpkm) > 0.5  ## both regions
geneRpkm = geneRpkm[geneIndex,]

## exon
exonIndex = rowMeans(exonRpkm) > 0.5
exonRpkm = exonRpkm[exonIndex,]

## junction
jxnIndex = rowMeans(jxnRp80m) > 0.5 & rowData(rse_jxn)$Class != "Novel"
jxnRp80m = jxnRp80m[jxnIndex,]

## transcript
txIndex = rowMeans(txTpm) > 0.5 
txTpm = txTpm[txIndex,]

rm(geneIndex, exonIndex, jxnIndex, txIndex)
print("....rpkms created....")

#######################
####### do PCA ########
#######################

# pcaGene = prcomp(t(log2(geneRpkm+1)))
# # kGene = num.sv(log2(geneRpkm+1), mod)
# # kGene = ifelse(kGene>25,25,kGene)
# kGene = 25
# genePCs = pcaGene$x[,1:kGene]

# pcaExon = prcomp(t(log2(exonRpkm+1)))
# # kExon = num.sv(log2(exonRpkm+1), mod, vfilter=50000)
# # kExon = ifelse(kExon>25,25,kExon)
# kExon = 25
# exonPCs = pcaExon$x[,1:kExon]

# pcaJxn = prcomp(t(log2(jxnRp80m+1)))
# # kJxn = num.sv(log2(jxnRp80m+1), mod, vfilter=50000)
# # kJxn = ifelse(kJxn>25,25,kJxn)
# kJxn = 25
# jxnPCs = pcaJxn$x[,1:kJxn]

# pcaTx = prcomp(t(log2(txTpm+1)))
# # kTx = num.sv(log2(txTpm+1), mod, vfilter=50000)
# # kTx = ifelse(kTx>25,25,kTx)
# kTx = 25
# txPCs = pcaTx$x[,1:kTx]

# save(genePCs, exonPCs, jxnPCs, txPCs, 
	# file="rdas/pcs_4features_combined_regions.rda")

load("rdas/pcs_4features_combined_regions.rda")

## make covs and move BrainRegion to end
modReg = grep("Region",colnames(mod))
covsGene = SlicedData$new(t(cbind(mod[,-c(1,modReg)], genePCs, mod[,modReg])))
covsExon = SlicedData$new(t(cbind(mod[,-c(1,modReg)], exonPCs, mod[,modReg])))
covsJxn = SlicedData$new(t(cbind(mod[,-c(1,modReg)], jxnPCs, mod[,modReg])))
covsTx = SlicedData$new(t(cbind(mod[,-c(1,modReg)], txPCs, mod[,modReg])))

## fix BrainRegion row label in covs
rownames(covsGene)[nrow(covsGene)] = rownames(covsExon)[nrow(covsExon)] = 
	rownames(covsJxn)[nrow(covsJxn)] = rownames(covsTx)[nrow(covsTx)] = colnames(mod)[modReg]



covsGene = SlicedData$new(t(cbind(mod[,-1],genePCs)))
covsExon = SlicedData$new(t(cbind(mod[,-1],exonPCs)))
covsJxn = SlicedData$new(t(cbind(mod[,-1],jxnPCs)))
covsTx = SlicedData$new(t(cbind(mod[,-1],txPCs)))

rm(genePCs, exonPCs, jxnPCs, txPCs)
print("....pcas created....")

##########################
### feature annotation ###
##########################

###### gene level
posGene = as.data.frame(rowRanges(rse_gene))[,1:3]
posGene$name = rownames(posGene)
posGene = posGene[,c(4,1:3)]

##### exon level 
posExon = as.data.frame(rowRanges(rse_exon))[,1:3]
posExon$name = rownames(posExon)
posExon = posExon[,c(4,1:3)]

##### junction level 
posJxn = as.data.frame(rowRanges(rse_jxn))[,1:3]
posJxn$name = rownames(posJxn)
posJxn = posJxn[,c(4,1:3)]
names(posJxn)[2:4] = c("Chr", "Start","End")

##### transcript level 
posTx = as.data.frame(rowRanges(rse_tx))[,1:3]
posTx$name = rownames(posTx)
posTx = posTx[,c(4,1:3)]
names(posTx)[2:4] = c("Chr", "Start","End")


#############################
### sliced expression data ##
geneSlice = SlicedData$new(log2(geneRpkm+1))
exonSlice = SlicedData$new(log2(exonRpkm+1))
jxnSlice = SlicedData$new(log2(jxnRp80m+1))
txSlice = SlicedData$new(log2(txTpm+1))

geneSlice$ResliceCombined(sliceSize = 5000)
exonSlice$ResliceCombined(sliceSize = 5000)
jxnSlice$ResliceCombined(sliceSize = 5000)
txSlice$ResliceCombined(sliceSize = 5000)


keep = c("theSnps","snpspos","geneSlice","covsGene","posGene","exonSlice","covsExon","posExon",
		"jxnSlice","covsJxn","posJxn","txSlice","covsTx","posTx")
rm(list=ls()[! ls() %in% keep])
print("....beginning eQTL analysis....")


##########################
### Run EQTLs ############
##########################

meGene = Matrix_eQTL_main(snps=theSnps, gene = geneSlice, 
	cvrt = covsGene, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 0.001,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posGene, 
	useModel = modelLINEAR_CROSS,	cisDist=2.5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	

save(meGene, file="eqtl_tables/matrixEqtl_output_interaction_gene.rda")
	
meExon = Matrix_eQTL_main(snps=theSnps, gene = exonSlice, 
	cvrt = covsExon, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 0.001,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posExon, 
	useModel = modelLINEAR_CROSS,	cisDist=2.5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	

save(meExon, file="eqtl_tables/matrixEqtl_output_interaction_exon.rda")
	
meJxn = Matrix_eQTL_main(snps=theSnps, gene = jxnSlice, 
	cvrt = covsJxn, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 0.001,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posJxn, 
	useModel = modelLINEAR_CROSS,	cisDist=2.5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	

save(meJxn,	file="eqtl_tables/matrixEqtl_output_interaction_jxn.rda")
	
meTx = Matrix_eQTL_main(snps=theSnps, gene = txSlice, 
	cvrt = covsTx, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 0.001,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posTx, 
	useModel = modelLINEAR_CROSS,	cisDist=2.5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	

save(meTx, file="eqtl_tables/matrixEqtl_output_interaction_tx.rda")
	
	
save(meGene, meExon, meJxn, meTx,
	file="eqtl_tables/matrixEqtl_output_interaction_4features.rda")

	
######################
###### annotate ######

# load("eqtl_tables/matrixEqtl_output_interaction_4features.rda")
# load("data/zandiHypde_bipolar_rseTx_n511.rda")
# load("data/zandiHypde_bipolar_rseJxn_n511.rda")
# load("data/zandiHypde_bipolar_rseExon_n511.rda")
# load("data/zandiHypde_bipolar_rseGene_n511.rda")

# extract
geneEqtl = meGene$cis$eqtls
geneEqtl$gene = as.character(geneEqtl$gene)
geneEqtl$snps = as.character(geneEqtl$snps)

exonEqtl = meExon$cis$eqtls
exonEqtl$gene = as.character(exonEqtl$gene)
exonEqtl$snps = as.character(exonEqtl$snps)

jxnEqtl = meJxn$cis$eqtls
jxnEqtl$gene = as.character(jxnEqtl$gene)
jxnEqtl$snps = as.character(jxnEqtl$snps)

txEqtl = meTx$cis$eqtls
txEqtl$gene = as.character(txEqtl$gene)
txEqtl$snps = as.character(txEqtl$snps)

################################
# add gene annotation info #####
################################

geneEqtl$Symbol = rowRanges(rse_gene)$Symbol[match(geneEqtl$gene, rownames(rse_gene))]
geneEqtl$EnsemblGeneID = rowRanges(rse_gene)$ensemblID[match(geneEqtl$gene, rownames(rse_gene))]
geneEqtl$Type = "Gene"
geneEqtl$Class = "InGen"
geneEqtl = DataFrame(geneEqtl)
# geneEqtl$gene_type = rowRanges(rse_gene)$gene_type[match(geneEqtl$gene, rownames(rse_gene))]

exonEqtl$Symbol = rowRanges(rse_exon)$Symbol[match(exonEqtl$gene, rownames(rse_exon))]
exonEqtl$EnsemblGeneID = rowRanges(rse_exon)$ensemblID[match(exonEqtl$gene, rownames(rse_exon))]
exonEqtl$Type = "Exon"
exonEqtl$Class = "InGen"
exonEqtl = DataFrame(exonEqtl)
# exonEqtl$gene_type = rowRanges(rse_exon)$gene_type[match(exonEqtl$gene, rownames(rse_exon))]

jxnEqtl$Symbol = rowRanges(rse_jxn)$newGeneSymbol[match(jxnEqtl$gene, rownames(rse_jxn))]
jxnEqtl$EnsemblGeneID = rowRanges(rse_jxn)$newGeneID[match(jxnEqtl$gene, rownames(rse_jxn))]
jxnEqtl$Type = "Jxn"
jxnEqtl$Class = rowRanges(rse_jxn)$Class[match(jxnEqtl$gene, rownames(rse_jxn))]
jxnEqtl = DataFrame(jxnEqtl)
# jxnEqtl$gene_type = rowRanges(rse_jxn)$gene_type[match(jxnEqtl$gene, rownames(rse_jxn))]

txEqtl$Symbol = rowRanges(rse_tx)$gene_name[match(txEqtl$gene, rownames(rse_tx))]
txEqtl$EnsemblGeneID = ss(rowRanges(rse_tx)$gene_id[match(txEqtl$gene, rownames(rse_tx))],"\\.",1)
txEqtl$Type = "Tx"
txEqtl$Class = "InGen"
txEqtl = DataFrame(txEqtl)
# txEqtl$gene_type = rowRanges(rse_tx)$gene_type[match(txEqtl$gene, rownames(rse_tx))]


# merge
allEqtl = rbind(geneEqtl, exonEqtl, jxnEqtl, txEqtl)
allEqtl$gencodeTx = CharacterList(c(as.list(rowRanges(rse_gene)$gencodeTx[match(geneEqtl$gene, 
	rownames(rse_gene))]),
	as.list(rowRanges(rse_exon)$gencodeTx[match(exonEqtl$gene, rownames(rse_exon))]),
	as.list(rowRanges(rse_jxn)$gencodeTx[match(jxnEqtl$gene, rownames(rse_jxn))]),
	as.list(txEqtl$gene)))
save(allEqtl, file="eqtl_tables/mergedEqtl_output_interaction_4features.rda",compress=TRUE)


# # #############
# # # metrics ###
# # sigEqtl = allEqtl[allEqtl$FDR < 0.01,] # fdr significant
# # length(unique(sigEqtl$EnsemblGeneID))

# # sigEqtlList = split(sigEqtl, factor(sigEqtl$Type,
	# # levels=c("Gene","Exon","Jxn", "Tx")))

# # sapply(sigEqtlList, function(x) max(x$pvalue))

# # sapply(sigEqtlList, function(x) length(unique(x$EnsemblGeneID)))
# # sapply(sigEqtlList, function(x) length(unique(x$Symbol)))
# # sapply(sigEqtlList, function(x) length(unique(x$snps)))

# # sapply(sigEqtlList, function(x) quantile(abs(x$beta)))

# # sapply(sigEqtlList, function(x) table(x$Class[!duplicated(x$gene)]))
# # sapply(sigEqtlList, function(x) prop.table(table(x$Class)))

