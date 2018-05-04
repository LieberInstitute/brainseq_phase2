####





corner = function(

corner = function(mat, n=6){
	if class(mat) == "list" {
		r = min(nrow(mat[[1]]), n)
		c = min(ncol(mat[[1]]), n)
		l = min(length(mat), n)
		lapply(mat[1:l], function(x) x[1:r,1:c] )
	}
	else {
		r = min(nrow(mat), n)
		c = min(ncol(mat), n)
		x[1:r,1:c]
	
	}
}



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

pd = colData(rse_gene)

## load SNP data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")
snpMap$pos_hg19 = paste0(snpMap$CHR, ":", snpMap$POS)

## drop rs10708380:150158001:TG:T (missing info in snpMap (and dbSNP))
snpInd = which(rownames(snpMap) == "rs10708380:150158001:TG:T")
snpMap = snpMap[-snpInd,]
snp = snp[-snpInd,]

#####################
# filter brain region
# make mds and snp dimensions equal to N
# (repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$RNum

## risk loci from PGC paper + rAggr proxy markers
riskLoci = read.csv("rAggr_results_163.csv", stringsAsFactors=FALSE)	# xxx snps
colnames(riskLoci) = gsub("\\.", "_", colnames(riskLoci))
riskLoci$hg19POS = paste0(riskLoci$SNP2_Chr, ":", riskLoci$SNP2_Pos) 

## keep SNPs from list
keepIndex = which(snpMap$pos_hg19 %in% riskLoci$hg19POS)	# keep 10,777 snps
snpMap = snpMap[keepIndex,]
snp = snp[keepIndex,]

snpMap$maf = rowSums(snp, na.rm=TRUE)/(2*rowSums(!is.na(snp))) 


######################
# statistical model ##
######################
pd$PrimaryDx = factor(pd$PrimaryDx,
	levels = c("Control", "Bipolar"))

mod = model.matrix(~PrimaryDx + Sex + as.matrix(mds[,1:5]), data = pd)
colnames(mod)[4:8] = colnames(mds)[1:5]


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

geneRpkm = assays(rse_gene)$rpkm
exonRpkm = assays(rse_exon)$rpkm
jxnRp10m = assays(rse_jxn)$rp10m
txTpm = assays(rse_tx)$tpm

pcaGene = prcomp(t(log2(geneRpkm+1)))
kGene = num.sv(log2(geneRpkm+1), mod)
genePCs = pcaGene$x[,1:kGene]

pcaExon = prcomp(t(log2(exonRpkm+1)))
kExon = num.sv(log2(exonRpkm+1), mod, vfilter=50000)
exonPCs = pcaExon$x[,1:kExon]

pcaJxn = prcomp(t(log2(jxnRp10m+1)))
kJxn = num.sv(log2(jxnRp10m+1), mod, vfilter=50000)
jxnPCs = pcaJxn$x[,1:kJxn]

pcaTx = prcomp(t(log2(txTpm+1)))
kTx = num.sv(log2(txTpm+1), mod, vfilter=50000)
txPCs = pcaTx$x[,1:kTx]

save(genePCs, exonPCs, jxnPCs, txPCs, 
	file="rdas/pcs_4features_amyg.rda")
# load("rdas/pcs_4features_amyg.rda")

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

print("Starting eQTLs")

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
	file="matrixEqtl_output_amyg_raggr_4features.rda")

	
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
save(allEqtl, file="mergedEqtl_output_amyg_raggr_4features.rda",compress=TRUE)









