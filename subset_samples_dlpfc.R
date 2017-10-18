####

## load packages
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)
library(RColorBrewer)

# make folder
dir.create("qcChecks")

## metrics
metrics = read.csv("preprocessed_data/DLPFC_RiboZero/read_and_alignment_metrics_DLPFC_RiboZero_BrainSeq_Phase2_LIBD.csv",
	as.is=TRUE, row.names=1)
metrics$RNum = ss(metrics$SAMPLE_ID, "_")

pheno = read.csv("preprocessed_data/DLPFC_PhaseII_sample-selection_04_20_2015.csv",as.is=TRUE)
pheno$RNum = paste0("R", pheno$RNum)
pheno$BrNum = paste0("Br", pheno$BrNum)
pheno = pheno[,1:8]

pd = cbind(metrics, pheno[match(metrics$RNum, pheno$RNum),])

# drop samples?
pdf("qcChecks/RIN_vs_overallMap_check_predrop_dlpfc.pdf")
par(mar=c(5,6,2,2),cex.axis=1.8,cex.lab=1.8)
palette(brewer.pal(8,"Set1"))
plot(pd$overallMapRate, pd$totalAssignedGene, 
	ylab = "Exonic Map Rate", cex=1.2,
	xlab = "Overall Mapping Rate", pch =21, bg="grey")
plot(pd$mitoRate, pd$totalAssignedGene, 
	ylab = "Exonic Map Rate", cex=1.2,
	xlab = "chrM Mapping Rate", pch =21, bg="grey")
dev.off()

## drop 38 samples
pd = pd[pd$totalAssignedGene > 0.2 & pd$mitoRate < 0.06,]
table(table(pd$RNum))

##########
## fix some IDs based on genotyping

## load
load("/dcl01/lieber/ajaffe/lab/brain_swap/genotype_match_matrix_allBrains.rda")
snpCorMat = snpCor2[,match(pd$RNum, ss(colnames(snpCor2),"_"))]
colnames(snpCorMat) = pd$RNum

matchInd = as.data.frame(which(snpCorMat > 0.6, arr.ind=TRUE, useNames=TRUE))
matchInd$GenoInfo = rownames(snpCorMat)[matchInd$row]
matchInd$GenoBrNum = ss(matchInd$GenoInfo, "_")
matchInd$RnaInfo = colnames(snpCorMat)[matchInd$col]
matchInd$RnaBrNum = pd$BrNum[match(matchInd$RnaInfo, pd$RNum)]

## add correlations
for(i in 1:nrow(matchInd)) matchInd$genoCor[i] = snpCorMat[matchInd$row[i], matchInd$col[i]]
matchInd = matchInd[,c(7, 3:6)]

# check duplicates
matchInd[matchInd$RnaInfo %in% matchInd$RnaInfo[
	matchInd$RnaBrNum != matchInd$GenoBrNum],]
                    # genoCor        GenoInfo GenoBrNum RnaInfo RnaBrNum
# Br1779_1M         0.8279711       Br1779_1M    Br1779   R3377   Br1794
# Br1794_1M         0.8279711       Br1794_1M    Br1794   R3377   Br1794
# Br1178_1M         0.9591357       Br1178_1M    Br1178   R3927   Br1179
# Br1179_1M         0.9591357       Br1179_1M    Br1179   R3927   Br1179
# Br1563_h650       0.9958781     Br1563_h650    Br1563   R5859   Br1563
# Br2407_1M         0.9958781       Br2407_1M    Br2407   R5859   Br1563
# Br2498_Omni2pt5   0.8595796 Br2498_Omni2pt5    Br2498  R12246   Br2448
# Br2498_Omni2pt5.1 0.8595796 Br2498_Omni2pt5    Br2498  R12246   Br2448
# Br2498_Omni2pt5.2 0.8595796 Br2498_Omni2pt5    Br2498  R12246   Br2448
# Br1051_1M         0.9744902       Br1051_1M    Br1051  R12309   Br1051
# Br11051_Omni5M    0.9241279  Br11051_Omni5M   Br11051  R12309   Br1051
# Br21051_Omni5M    0.9273723  Br21051_Omni5M   Br21051  R12309   Br1051
# Br1060_h650       0.9588512     Br1060_h650    Br1060  R12313   Br1061
# Br1119_Omni5M     0.9015160   Br1119_Omni5M    Br1119  R12324   Br1119
# Br21119_Omni5M    0.9152229  Br21119_Omni5M   Br21119  R12324   Br1119
# Br1119_Omni5M.1   0.9015160   Br1119_Omni5M    Br1119  R12324   Br1119
# Br21119_Omni5M.1  0.9152229  Br21119_Omni5M   Br21119  R12324   Br1119
# Br1142_h650       0.9622655     Br1142_h650    Br1142  R12330   Br1143
# Br1143_1M         0.9666436       Br1143_1M    Br1143  R12330   Br1143
# Br1143_h650       0.9666436     Br1143_h650    Br1143  R12330   Br1143
	
## drop Br1794, Br1143, Br1142, Br1779, Br1178, Br2407 Br2448, Br1563 

pd = pd[ ! pd$BrNum %in% c("Br1794", "Br1143", "Br1142", "Br1779", 
	"Br1178", "Br2407", "Br2448", "Br1563"),]
pd = pd[pd$RNum %in% matchInd$RnaInfo,] # 8 brains

write.table(matchInd, file="qcChecks/genotype_checking_table_DLPFC.txt",
	sep="\t", row.names=FALSE, col.names=FALSE)

########################
### load counts ########
load("preprocessed_data/DLPFC_RiboZero/rpkmCounts_DLPFC_RiboZero_BrainSeq_Phase2_LIBD_n579.rda")
load("preprocessed_data/DLPFC_RiboZero/rawCounts_DLPFC_RiboZero_BrainSeq_Phase2_LIBD_n579.rda")

## filter to same samples
geneCounts = geneCounts[,pd$SAMPLE_ID]
exonCounts = exonCounts[,pd$SAMPLE_ID]

# filter
jIndex = which(jMap$meanExprs > 0.05)
jCounts = jCounts[jIndex,]
jMap = jMap[jIndex,]

## update junction annotation class
tt = rowSums(!is.na(as.data.frame(mcols(jMap)[,
	c("startExon", "endExon")])))
jMap$Class[tt == 1] = "AltStartEnd"
jMap$Class[tt == 2 & !jMap$inGencode] = "ExonSkip"

## add transcripts
txPath = "preprocessed_data/DLPFC_RiboZero/Salmon_tx/"
txFiles = paste0(txPath, pd$SAMPLE_ID, "/quant.sf")
names(txFiles) = rownames(pd)
table(file.exists(txFiles))

txList = lapply(txFiles, read.delim, row.names=1, as.is=TRUE)
tpmMat = sapply(txList, "[[", "TPM")
rownames(tpmMat) = ss(rownames(txList[[1]]),"|", fixed=TRUE)

##################
## compare full counts
pca = prcomp(t(log2(geneCounts+1)))
plot(pca$x) # looks fine

###############
### merge
sIndexes = splitit(pd$RNum)
geneCountMerge = sapply(sIndexes, function(ii) {
	rowSums(t(t(geneCounts[,ii])))
})

exonCountMerge = sapply(sIndexes, function(ii) {
	rowSums(t(t(exonCounts[,ii])))
})

jCountMerge = sapply(sIndexes, function(ii) {
	rowSums(t(t(jCounts[,ii])))
})

## transcripts
combineTxQuant = function (sampIDs, salmonDir) {    
    ##get names of transcripts
    txNames = read.table(file.path(salmonDir, sampIDs[1], "quant.sf"),header = TRUE)$Name
    txNames = as.character(txNames)
    gencodeTx = ss(txNames, "\\|",1)
    
	## read in
	txList = lapply(paste0(salmonDir, sampIDs, "/quant.sf"), read.table, header=TRUE)
	
	## initiate output
    txQuant = data.frame(EffectiveLength=rep(NA,length(gencodeTx)), 
		NumReads=NA, RPK=NA, TPM=NA, row.names=gencodeTx)
    
	## run calculations
    txQuant$EffectiveLength = rowSums(sapply(txList, "[[", "EffectiveLength"))
    txQuant$NumReads = rowSums(sapply(txList, "[[", "NumReads"))
    txQuant$RPK = txQuant$NumReads/txQuant$EffectiveLength
    txQuant$TPM = 1e6 * (txQuant$NumReads/txQuant$EffectiveLength) / sum(txQuant$RPK)

    return(txQuant)
}

## run for multiples
doubleSamples = lapply(sIndexes[lengths(sIndexes)>1], function(ii) pd$SAMPLE_ID[ii])
doubleTxCombine = lapply(doubleSamples, combineTxQuant, salmonDir=txPath)

doubleTx = sapply(doubleTxCombine, "[[", "TPM")
singleTx = tpmMat[,unlist(sIndexes[lengths(sIndexes) == 1])]
tpmMatMerge = cbind(singleTx, doubleTx)
colnames(tpmMatMerge) = ss(colnames(tpmMatMerge), "_")
tpmMatMerge = tpmMatMerge[,names(sIndexes)]

############
## make expression sets
library(SummarizedExperiment)
pdDF = DataFrame(pd)
classIndex = splitit(sapply(pdDF, class))

pdMergeList = sapply(sIndexes, function(ii) {
	x = pdDF[ii,]
	o = pdDF[ii[1],]
	for(j in classIndex$character) o[j] = CharacterList(unique(x[,j]))
	for(j in classIndex$integer) o[j] = IntegerList(unique(x[,j]))
	for(j in classIndex$logical) o[j] = LogicalList(unique(x[,j]))
	for(j in classIndex$numeric) o[j] = NumericList(unique(x[,j]))
	return(o)
})
pdMerge = do.call("rbind", pdMergeList)		
pdMerge$RNum.1 = NULL
colnames(pdMerge)[64] = "Age"

## clean
pdMerge$RNum = sapply(pdMerge$RNum, "[[", 1)
pdMerge$BrNum = sapply(pdMerge$BrNum, "[[", 1)
pdMerge$Age = sapply(pdMerge$Age, "[[", 1)
pdMerge$Sex = sapply(pdMerge$Sex, "[[", 1)
pdMerge$Race = sapply(pdMerge$Race, "[[", 1)
pdMerge$Dx = sapply(pdMerge$Dx, "[[", 1)
pdMerge$Region = sapply(pdMerge$Region, "[[", 1)
rownames(pdMerge) = pdMerge$RNum
identical(rownames(pdMerge), colnames(geneCountMerge))

##########
## gene
geneMapGR = makeGRangesFromDataFrame(geneMap, keep=TRUE)
geneMapGR$gencodeTx = CharacterList(strsplit(geneMapGR$gencodeTx, ";"))
rse_gene = SummarizedExperiment(
	assays = list('counts' = geneCountMerge),
    colData = pdMerge, rowRanges = geneMapGR)

## gene
exonMapGR = makeGRangesFromDataFrame(exonMap, keep=TRUE)
exonMapGR$gencodeTx = CharacterList(strsplit(exonMapGR$gencodeTx, ";"))
rse_exon = SummarizedExperiment(
	assays = list('counts' = exonCountMerge),
    colData = pdMerge, rowRanges = exonMapGR)
	
## junction
rse_jxn = SummarizedExperiment(
	assays = list('counts' = jCountMerge),
    colData = pdMerge, rowRanges = jMap)
	
## transcript
gtf = import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf")
tx = gtf[which(gtf$type == "transcript")]
names(tx) = tx$transcript_id
txMap = tx[rownames(tpmMat)]
rse_tx = SummarizedExperiment(
	assays = list('tpm' = tpmMatMerge),
    colData = pdMerge, rowRanges = txMap)
	
#### function for RPKM
getRPKM = function(rse) {
	require(SummarizedExperiment)
	bg = matrix(rep(colData(rse)$totalMapped), 
		nc = ncol(rse), nr = nrow(rse),	byrow=TRUE)
	wid = matrix(rep(rowData(rse)$Length), 
		nr = nrow(rse), nc = ncol(rse),	byrow=FALSE)
	assays(rse)$counts/(wid/1000)/(bg/1e6)
}

getRPM = function(rse, target = 80e6) {
	require(SummarizedExperiment)
	bg = matrix(rep(colData(rse)$totalMapped/target), 
		nc = ncol(rse), nr = nrow(rse),	byrow=TRUE)
	assays(rse)$counts/bg
}

################
## save ########
save(rse_gene, getRPKM, compress=TRUE,
	file = "count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseGene_merged_n449.rda")
save(rse_exon, getRPKM, compress=TRUE,
	file = "count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseExon_merged_n449.rda")
save(rse_jxn, getRPM, compress=TRUE,
	file = "count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseJxn_merged_n449.rda")
save(rse_tx, compress=TRUE,
	file = "count_data/dlpfc_ribozero_brainseq_phase2_hg38_rseTx_merged_n449.rda")

