####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library('sessioninfo')

######################
### load data ####
######################

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata")

## keep adult samples & correct region
keepInd = which(colData(rse_gene)$Age > 13 & colData(rse_gene)$Region == "DLPFC" & colData(rse_gene)$Race == 'CAUC')
rse_gene = rse_gene[,keepInd]
dim(rse_gene)


## extract pd and rpkms
pd = colData(rse_gene)
geneRpkm = assays(rse_gene)$rpkm


######################
### snp data ####
######################

## load SNP data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda")

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

mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]),
	data = pd)
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

pcaGene = prcomp(t(log2(geneRpkm+1)))
kGene = num.sv(log2(geneRpkm+1), mod)
genePCs = pcaGene$x[,1:kGene]

save(genePCs,
	file="rdas/pcs_dlpfc_gene_filtered_over13.rda")

covsGene = SlicedData$new(t(cbind(mod[,-1],genePCs)))


##########################
### feature annotation ###
##########################

###### gene level
posGene = as.data.frame(rowRanges(rse_gene))[,1:3]
posGene$name = rownames(posGene)
posGene = posGene[,c(4,1:3)]

#############################
### sliced expression data ##
geneSlice = SlicedData$new(log2(geneRpkm+1))

geneSlice$ResliceCombined(sliceSize = 5000)

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

save(meGene,
	file="eqtl_tables/matrixEqtl_output_dlpfc_gene.rda")

	
######################
###### annotate ######

# extract
geneEqtl = meGene$cis$eqtls
geneEqtl$gene = as.character(geneEqtl$gene)
geneEqtl$snps = as.character(geneEqtl$snps)

################################
# add gene annotation info #####
################################

geneEqtl$Symbol = rowRanges(rse_gene)$Symbol[match(geneEqtl$gene, rownames(rse_gene))]
geneEqtl$EnsemblGeneID = rowRanges(rse_gene)$ensemblID[match(geneEqtl$gene, rownames(rse_gene))]
geneEqtl$Type = "Gene"
geneEqtl$Class = "InGen"
geneEqtl = DataFrame(geneEqtl)
geneEqtl$gene_type = rowRanges(rse_gene)$gene_type[match(geneEqtl$gene, rownames(rse_gene))]

geneEqtl$gencodeTx = CharacterList(c(as.list(rowRanges(rse_gene)$gencodeTx[match(geneEqtl$gene, 
	rownames(rse_gene))])))

save(geneEqtl, file="eqtl_tables/mergedEqtl_output_dlpfc_gene.rda",compress=TRUE)


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

