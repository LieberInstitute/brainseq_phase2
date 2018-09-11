### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library('devtools')

dir.create('rdas', showWarnings = FALSE)
dir.create('eqtl_tables', showWarnings = FALSE)

######################
### load data ####
######################
if(!file.exists('eqtl_tables/matrixEqtl_output_interaction_4features.rda')) {

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/rdas/rse_gtex_gene_subset.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/rdas/rse_gtex_exon_subset.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/rdas/rse_gtex_jxn_subset.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/rdas/rse_gtex_tx_subset.Rdata")

# # fix junction row names
# rownames(rse_gtex_jxn) = paste0(seqnames(rse_gtex_jxn),":",start(rse_gtex_jxn),"-",end(rse_gtex_jxn),"(",strand(rse_gtex_jxn),")")

# # sum totalMapped IntegerLists (so getRPKM works later)
# colData(rse_gtex_tx)$totalMapped =
	# colData(rse_gtex_jxn)$totalMapped =
	# colData(rse_gtex_exon)$totalMapped = 
	# colData(rse_gtex_gene)$totalMapped  = sapply(colData(rse_gtex_gene)$totalMapped, sum)

## keep adult samples - keep both regions
print(dim(rse_gtex_gene))
## Not really filtering by age here
summary(colData(rse_gtex_gene)$age)
keepInd = which(colData(rse_gtex_gene)$age > 13)
print(length(keepInd))
rse_gtex_gene = rse_gtex_gene[,keepInd]
rse_gtex_exon = rse_gtex_exon[,keepInd]
rse_gtex_jxn = rse_gtex_jxn[,keepInd]
rse_gtex_tx = rse_gtex_tx[,keepInd]

## extract pd and rpkms
pd = colData(rse_gtex_gene)
geneRpkm = assays(rse_gtex_gene)$rpkm
exonRpkm = assays(rse_gtex_exon)$rpkm
jxnRp10m = assays(rse_gtex_jxn)$rp10m
txTpm = assays(rse_gtex_tx)$tpm

print("....data loaded....")



######################
### snp data ####
######################

## load SNP data
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/genotypeData_GTEx_hippoPlusDlpfc_simplified.rda', verbose = TRUE)

### make mds and snp dimensions equal to N
m <- match(rownames(pd), pdGtex$sra_accession)
## Same number and location of NAs with
# m <- match(pd$sampid, pdGtex$SAMPID)
# m <- match(rownames(pd), pdGtex$sra_accession)
## This one has a few different matches, but NA locations are identical.
# m <- match(ss(as.character(pd$subjid), '-', 2), pdGtex$SUBJID)
m1 <- which(!is.na(m))
m2 <- m[!is.na(m)]

## Number of samples dropped at this stage and remaining
print(length(m) - length(m1))
print(length(m1))

pd <- pd[m1, ]
rse_gtex_gene <- rse_gtex_gene[, m1]
rse_gtex_exon <- rse_gtex_exon[, m1]
rse_gtex_jxn <- rse_gtex_jxn[, m1]
rse_gtex_tx <- rse_gtex_tx[, m1]
geneRpkm <- geneRpkm[, m1]
exonRpkm <- exonRpkm[, m1]
jxnRp10m <- jxnRp10m[, m1]
txTpm <- txTpm[, m1]

pdGtex <- pdGtex[m2, ]
snpGtex <- snpGtex[, m2]

print(dim(pdGtex))
print(dim(snpGtex))

colnames(snpGtex) = rownames(pd)

######################
# statistical model ##
######################

mod = model.matrix(~ Sex + as.matrix(pdGtex[, paste0('snpPC', 1:5)]) + Region,
	data = pd)
colnames(mod)[3:7] = paste0('snpPC', 1:5)

######################
# create SNP objects #
######################

theSnps = SlicedData$new(as.matrix(snpGtex))
theSnps$ResliceCombined(sliceSize = 50000)

snpspos = snpMapGtex[,c("SNP","chr_hg38","pos_hg38")]
colnames(snpspos) = c("name","chr","pos")

#######################
####### do PCA ########
#######################

pcaGene = prcomp(t(log2(geneRpkm+1)))
kGene = num.sv(log2(geneRpkm+1), mod)
kGene = min(kGene, 25)
genePCs = pcaGene$x[,1:kGene]

pcaExon = prcomp(t(log2(exonRpkm+1)))
kExon = num.sv(log2(exonRpkm+1), mod, vfilter=50000)
kExon = min(kExon, 25)
exonPCs = pcaExon$x[,1:kExon]

pcaJxn = prcomp(t(log2(jxnRp10m+1)))
kJxn = num.sv(log2(jxnRp10m+1), mod, vfilter=50000)
kJxn = min(kJxn, 25)
jxnPCs = pcaJxn$x[,1:kJxn]

pcaTx = prcomp(t(log2(txTpm+1)))
kTx = num.sv(log2(txTpm+1), mod, vfilter=50000)
kTx = min(kTx, 25)
txPCs = pcaTx$x[,1:kTx]

print('k info for the PCAs')
print(c('kGene' = kGene, 'kExon' = kExon, 'kJxn' = kJxn, 'kTx' = kTx))

save(genePCs, exonPCs, jxnPCs, txPCs, 
	file="rdas/pcs_4features_combined_regions_filtered_over13.rda")
# load("rdas/pcs_4features_combined_regions_filtered_over13.rda")

## make covs and move BrainRegion to end
modReg = grep("Region",colnames(mod))
covsGene = SlicedData$new(t(cbind(mod[,-c(1,modReg)], genePCs, mod[,modReg])))
covsExon = SlicedData$new(t(cbind(mod[,-c(1,modReg)], exonPCs, mod[,modReg])))
covsJxn = SlicedData$new(t(cbind(mod[,-c(1,modReg)], jxnPCs, mod[,modReg])))
covsTx = SlicedData$new(t(cbind(mod[,-c(1,modReg)], txPCs, mod[,modReg])))

## fix BrainRegion row label in covs
rownames(covsGene)[nrow(covsGene)] = rownames(covsExon)[nrow(covsExon)] = 
	rownames(covsJxn)[nrow(covsJxn)] = rownames(covsTx)[nrow(covsTx)] = colnames(mod)[modReg]

# covsGene = SlicedData$new(t(cbind(mod[,-1],genePCs)))
# covsExon = SlicedData$new(t(cbind(mod[,-1],exonPCs)))
# covsJxn = SlicedData$new(t(cbind(mod[,-1],jxnPCs)))
# covsTx = SlicedData$new(t(cbind(mod[,-1],txPCs)))

rm(genePCs, exonPCs, jxnPCs, txPCs)
print("....pcas created....")

##########################
### feature annotation ###
##########################

###### gene level
posGene = as.data.frame(rowRanges(rse_gtex_gene))[,1:3]
posGene$name = rownames(posGene)
posGene = posGene[,c(4,1:3)]

##### exon level 
posExon = as.data.frame(rowRanges(rse_gtex_exon))[,1:3]
posExon$name = rownames(posExon)
posExon = posExon[,c(4,1:3)]

##### junction level 
posJxn = as.data.frame(rowRanges(rse_gtex_jxn))[,1:3]
posJxn$name = rownames(posJxn)
posJxn = posJxn[,c(4,1:3)]
names(posJxn)[2:4] = c("Chr", "Start","End")

##### transcript level 
posTx = as.data.frame(rowRanges(rse_gtex_tx))[,1:3]
posTx$name = rownames(posTx)
posTx = posTx[,c(4,1:3)]
names(posTx)[2:4] = c("Chr", "Start","End")


#############################
### sliced expression data ##
geneSlice = SlicedData$new(log2(geneRpkm+1))
exonSlice = SlicedData$new(log2(exonRpkm+1))
jxnSlice = SlicedData$new(log2(jxnRp10m+1))
txSlice = SlicedData$new(log2(txTpm+1))

geneSlice$ResliceCombined(sliceSize = 5000)
exonSlice$ResliceCombined(sliceSize = 5000)
jxnSlice$ResliceCombined(sliceSize = 5000)
txSlice$ResliceCombined(sliceSize = 5000)


#keep = c("theSnps","snpspos","geneSlice","covsGene","posGene","exonSlice","covsExon","posExon",
#		"jxnSlice","covsJxn","posJxn","txSlice","covsTx","posTx")
#rm(list=ls()[! ls() %in% keep])
print("....beginning eQTL analysis....")


##########################
### Run EQTLs ############
##########################

meGene = Matrix_eQTL_main(snps=theSnps, gene = geneSlice, 
	cvrt = covsGene, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posGene, 
	useModel = modelLINEAR_CROSS,	cisDist=2.5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	
save(meGene, file="eqtl_tables/matrixEqtl_output_interaction_gene.rda")
	
meExon = Matrix_eQTL_main(snps=theSnps, gene = exonSlice, 
	cvrt = covsExon, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posExon, 
	useModel = modelLINEAR_CROSS,	cisDist=2.5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	
save(meExon, file="eqtl_tables/matrixEqtl_output_interaction_exon.rda")
	
meJxn = Matrix_eQTL_main(snps=theSnps, gene = jxnSlice, 
	cvrt = covsJxn, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posJxn, 
	useModel = modelLINEAR_CROSS,	cisDist=2.5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	
save(meJxn,	file="eqtl_tables/matrixEqtl_output_interaction_jxn.rda")
	
meTx = Matrix_eQTL_main(snps=theSnps, gene = txSlice, 
	cvrt = covsTx, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posTx, 
	useModel = modelLINEAR_CROSS,	cisDist=2.5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	
save(meTx, file="eqtl_tables/matrixEqtl_output_interaction_tx.rda")
	
	
save(meGene, meExon, meJxn, meTx,
	file="eqtl_tables/matrixEqtl_output_interaction_4features.rda")

	
######################
###### annotate ######

} else {
load("eqtl_tables/matrixEqtl_output_interaction_4features.rda", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/rdas/rse_gtex_gene_subset.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/rdas/rse_gtex_exon_subset.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/rdas/rse_gtex_jxn_subset.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/rdas/rse_gtex_tx_subset.Rdata", verbose = TRUE)
}

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

geneEqtl$Symbol = rowRanges(rse_gtex_gene)$Symbol[match(geneEqtl$gene, rownames(rse_gtex_gene))]
geneEqtl$EnsemblGeneID = rowRanges(rse_gtex_gene)$ensemblID[match(geneEqtl$gene, rownames(rse_gtex_gene))]
geneEqtl$Type = "Gene"
geneEqtl$Class = "InGen"
geneEqtl = DataFrame(geneEqtl)
# geneEqtl$gene_type = rowRanges(rse_gtex_gene)$gene_type[match(geneEqtl$gene, rownames(rse_gtex_gene))]

exonEqtl$Symbol = rowRanges(rse_gtex_exon)$Symbol[match(exonEqtl$gene, rownames(rse_gtex_exon))]
exonEqtl$EnsemblGeneID = rowRanges(rse_gtex_exon)$ensemblID[match(exonEqtl$gene, rownames(rse_gtex_exon))]
exonEqtl$Type = "Exon"
exonEqtl$Class = "InGen"
exonEqtl = DataFrame(exonEqtl)
# exonEqtl$gene_type = rowRanges(rse_gtex_exon)$gene_type[match(exonEqtl$gene, rownames(rse_gtex_exon))]

jxnEqtl$Symbol = rowRanges(rse_gtex_jxn)$newGeneSymbol[match(jxnEqtl$gene, rownames(rse_gtex_jxn))]
jxnEqtl$EnsemblGeneID = rowRanges(rse_gtex_jxn)$newGeneID[match(jxnEqtl$gene, rownames(rse_gtex_jxn))]
jxnEqtl$Type = "Jxn"
jxnEqtl$Class = rowRanges(rse_gtex_jxn)$Class[match(jxnEqtl$gene, rownames(rse_gtex_jxn))]
jxnEqtl = DataFrame(jxnEqtl)
# jxnEqtl$gene_type = rowRanges(rse_gtex_jxn)$gene_type[match(jxnEqtl$gene, rownames(rse_gtex_jxn))]

txEqtl$Symbol = rowRanges(rse_gtex_tx)$gene_name[match(txEqtl$gene, rownames(rse_gtex_tx))]
txEqtl$EnsemblGeneID = ss(rowRanges(rse_gtex_tx)$gene_id[match(txEqtl$gene, rownames(rse_gtex_tx))],"\\.",1)
txEqtl$Type = "Tx"
txEqtl$Class = "InGen"
txEqtl = DataFrame(txEqtl)
# txEqtl$gene_type = rowRanges(rse_gtex_tx)$gene_type[match(txEqtl$gene, rownames(rse_gtex_tx))]


# merge
allEqtl = rbind(geneEqtl, exonEqtl, jxnEqtl, txEqtl)
# allEqtl$gencodeTx = CharacterList(c(as.list(rowRanges(rse_gtex_gene)$gencodeTx[match(geneEqtl$gene,
#     rownames(rse_gtex_gene))]),
#     as.list(rowRanges(rse_gtex_exon)$gencodeTx[match(exonEqtl$gene, rownames(rse_gtex_exon))]),
#     as.list(rowRanges(rse_gtex_jxn)$gencodeTx[match(jxnEqtl$gene, rownames(rse_gtex_jxn))]),
#     as.list(txEqtl$gene)))
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

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

