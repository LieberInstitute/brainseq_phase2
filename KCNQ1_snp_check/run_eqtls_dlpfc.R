####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)

######################
### load data ####
######################

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_exon.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_jxn.Rdata")
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_tx.Rdata")

# # fix junction row names
# rownames(rse_jxn) = paste0(seqnames(rse_jxn),":",start(rse_jxn),"-",end(rse_jxn),"(",strand(rse_jxn),")")

# # sum totalMapped IntegerLists (so getRPKM works later)
# colData(rse_tx)$totalMapped =
	# colData(rse_jxn)$totalMapped =
	# colData(rse_exon)$totalMapped = 
	# colData(rse_gene)$totalMapped  = sapply(colData(rse_gene)$totalMapped, sum)

## keep adult samples & correct region
keepInd = which(colData(rse_gene)$Age > 13 & colData(rse_gene)$Region == "DLPFC")
rse_gene = rse_gene[,keepInd]
rse_exon = rse_exon[,keepInd]
rse_jxn = rse_jxn[,keepInd]
rse_tx = rse_tx[,keepInd]

## extract pd and rpkms
pd = colData(rse_gene)
geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm


######################
### snp data ####
######################

## load SNP data
load("BrainSeq_Phase2_chr11_n551.rda")

### make mds and snp dimensions equal to N
###(repeat rows or columns for BrNum replicates)
# mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
# rownames(mds) = colnames(snp) = pd$RNum
colnames(snp) = pd$RNum


## keep only 1 snp
snpInd = grep("rs34097980", snpMap$SNP)
snpMap = snpMap[snpInd:(snpInd+1),]
snp = snp[snpInd:(snpInd+1),]

table(snp[1,])/ncol(snp)	## MAF = 0.053

# ## drop SNPs not mapping to hg38
# keepIndex = which(!is.na(snpMap$chr_hg38))
# snpMap = snpMap[keepIndex,]
# snp = snp[keepIndex,]

snpMap$name = ss(snpMap$SNP,":")
rownames(snp) = snpMap$SNP
snpMap$chr_hg38 = "chr11"
snpMap$pos_hg38 = c(2854651,2854735)	# from dbSNP

######################
# statistical model ##
######################

pd$Dx = factor(pd$Dx,
	levels = c("Control", "Schizo"))

# mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]),
	# data = pd)
# colnames(mod)[4:8] = colnames(mds)[1:5]

mod = model.matrix(~Dx + Sex, data = pd)

######################
# create SNP objects #
######################

theSnps = SlicedData$new(as.matrix(snp))
theSnps$ResliceCombined(sliceSize = 50000)

snpspos = snpMap[,c("SNP","chr_hg38","pos_hg38")]
colnames(snpspos) = c("name","chr","pos")


#######################
####### do PCA ########
#######################

# pcaGene = prcomp(t(log2(geneRpkm+1)))
# kGene = num.sv(log2(geneRpkm+1), mod)
# genePCs = pcaGene$x[,1:kGene]

# pcaExon = prcomp(t(log2(exonRpkm+1)))
# kExon = num.sv(log2(exonRpkm+1), mod, vfilter=50000)
# exonPCs = pcaExon$x[,1:kExon]

# pcaJxn = prcomp(t(log2(jxnRp10m+1)))
# kJxn = num.sv(log2(jxnRp10m+1), mod, vfilter=50000)
# jxnPCs = pcaJxn$x[,1:kJxn]

# pcaTx = prcomp(t(log2(txTpm+1)))
# kTx = num.sv(log2(txTpm+1), mod, vfilter=50000)
# txPCs = pcaTx$x[,1:kTx]

# save(genePCs, exonPCs, jxnPCs, txPCs, 
	# file="rdas/pcs_dlpfc_4features_filtered_over13.rda")
load("../eQTL_full/eqtl_tables/rdas/pcs_dlpfc_4features_filtered_over13.rda")

covsGene = SlicedData$new(t(cbind(mod[,-1],genePCs)))
covsExon = SlicedData$new(t(cbind(mod[,-1],exonPCs)))
covsJxn = SlicedData$new(t(cbind(mod[,-1],jxnPCs)))
covsTx = SlicedData$new(t(cbind(mod[,-1],txPCs)))

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
jxnSlice = SlicedData$new(log2(jxnRp10m+1))
txSlice = SlicedData$new(log2(txTpm+1))

geneSlice$ResliceCombined(sliceSize = 5000)
exonSlice$ResliceCombined(sliceSize = 5000)
jxnSlice$ResliceCombined(sliceSize = 5000)
txSlice$ResliceCombined(sliceSize = 5000)


##########################
### Run EQTLs ############
##########################
print("Begin eQTL analysis")

meGene = Matrix_eQTL_main(snps=theSnps, gene = geneSlice, 
	cvrt = covsGene, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posGene, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	

meExon = Matrix_eQTL_main(snps=theSnps, gene = exonSlice, 
	cvrt = covsExon, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posExon, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	

meJxn = Matrix_eQTL_main(snps=theSnps, gene = jxnSlice, 
	cvrt = covsJxn, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posJxn, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	
	
meTx = Matrix_eQTL_main(snps=theSnps, gene = txSlice, 
	cvrt = covsTx, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posTx, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)	

save(meGene, meExon, meJxn, meTx,
	file="matrixEqtl_output_dlpfc_rs34097980.rda")

	
######################
###### annotate ######

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
allEqtl = allEqtl[allEqtl$snps=="rs34097980:2875881:T:C",]
allEqtl$FDR = p.adjust(allEqtl$pvalue, "fdr")
save(allEqtl, file="mergedEqtl_output_dlpfc_rs34097980.rda",compress=TRUE)








