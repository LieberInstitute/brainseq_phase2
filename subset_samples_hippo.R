####

## load packages
library(jaffelab)
library(GenomicRanges)
library(rtracklayer)
library(RColorBrewer)

# make folder
dir.create("qcChecks")

## metrics
pd = read.csv("preprocessed_data/Hippo_RiboZero/read_and_alignment_mets_brainseq_hippo_n486.csv",
	as.is=TRUE, row.names=1)
pd$BrNum = paste0("Br", pd$BrNum)
colnames(pd)[2] = "RNum"

# drop samples?
pdf("qcChecks/RIN_vs_overallMap_check_predrop.pdf")
par(mar=c(5,6,2,2),cex.axis=1.8,cex.lab=1.8)
palette(brewer.pal(8,"Set1"))
plot(pd$RIN, pd$overallMapRate, xlab = "RIN", cex=1.2,
	ylab = "Overall Mapping Rate", pch =21, bg="grey")
abline(h=0.7, lty=2)
dev.off()

## drop 4 samples
pd = pd[pd$overallMapRate > 0.7,]

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
              # genoCor    GenoInfo GenoBrNum RnaInfo RnaBrNum
# Br1142_h650 0.9911928 Br1142_h650    Br1142   R4825   Br1143
# Br1143_1M   0.9956069   Br1143_1M    Br1143   R4825   Br1143
# Br1143_h650 0.9956069 Br1143_h650    Br1143   R4825   Br1143
# Br1178_1M   0.9447640   Br1178_1M    Br1178   R4827   Br1179
# Br1179_1M   0.9447640   Br1179_1M    Br1179   R4827   Br1179
# Br1779_1M   0.8465863   Br1779_1M    Br1779   R4915   Br1794
# Br1794_1M   0.8465863   Br1794_1M    Br1794   R4915   Br1794
# Br1563_h650 0.9586196 Br1563_h650    Br1563   R5497   Br1563
# Br2407_1M   0.9586196   Br2407_1M    Br2407   R5497   Br1563
	
## drop Br1142,Br1143, Br1779, Br1178, Br1563, Br2407 genotype 
pd = pd[ ! pd$BrNum %in% c("Br1794", "Br1143", "Br1142", "Br1779", 
	"Br1178", "Br2407", "Br2448", "Br1563"),]

pd = pd[pd$RNum %in% matchInd$RnaInfo,] # one sample

## dropped
                                                # SAMPLE_ID  RNum  BrNum   Age
# R4757_C04PTACXX_CGCTCATT_L00 R4757_C04PTACXX_CGCTCATT_L00 R4757 Br1060 75.52
                             # Sex Race      Dx PhaseII_HIPPO RnumHIPPOPhaseII
# R4757_C04PTACXX_CGCTCATT_L00   F CAUC Control          TRUE             4757

write.table(matchInd, file="qcChecks/genotype_checking_table_hippo.txt",
	sep="\t", row.names=FALSE, col.names=FALSE)

########################
### load counts ########
load("preprocessed_data/Hippo_RiboZero/rawCounts_HIPPO_RiboZero_BrainSeq_Phase2_LIBD_n486.rda")
load("preprocessed_data/Hippo_RiboZero/rpkmCounts_HIPPO_RiboZero_BrainSeq_Phase2_LIBD_n486.rda")

## filter to same samples
geneCounts = geneCounts[,pd$SAMPLE_ID]
exonCounts = exonCounts[,pd$SAMPLE_ID]

# filter
jIndex = which(jMap$meanExprs > 0.05)
jCounts = jCounts[jIndex,]
jMap = jMap[jIndex,]

# convert to matrix
jCounts = as.matrix(as.data.frame(jCounts[,pd$SAMPLE_ID]))

## update junction annotation class
tt = rowSums(!is.na(as.data.frame(mcols(jMap)[,
	c("startExon", "endExon")])))
jMap$Class[tt == 1] = "AltStartEnd"
jMap$Class[tt == 2 & !jMap$inGencode] = "ExonSkip"

## add transcripts
txPath = "preprocessed_data/Hippo_RiboZero/Salmon_tx/"
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
combineTxQuant = function (sampID1, sampID2, salmonDir)
{    
    ##get names of transcripts
    txNames = read.table(file.path(salmonDir, sampID1, "quant.sf"),header = TRUE)$Name
    txNames = as.character(txNames)
    gencodeTx = ss(txNames, "\\|",1)
    
    txQuant1 = read.table(file.path(salmonDir, sampID1, "quant.sf"),header = TRUE)
    txQuant2 = read.table(file.path(salmonDir, sampID2, "quant.sf"),header = TRUE)
    
    txQuant = data.frame(EffectiveLength=rep(NA,length(gencodeTx)), NumReads=NA, RPK=NA, TPM=NA, row.names=gencodeTx)
    
    txQuant$EffectiveLength = txQuant1$EffectiveLength + txQuant2$EffectiveLength
    txQuant$NumReads = txQuant1$NumReads + txQuant2$NumReads
    txQuant$RPK = txQuant$NumReads/txQuant$EffectiveLength
    txQuant$TPM = 1e6 * (txQuant$NumReads/txQuant$EffectiveLength) / sum(txQuant$RPK)

    return(txQuant)
}
doubleSamples = lapply(sIndexes[lengths(sIndexes)>1], function(ii) pd$SAMPLE_ID[ii])
doubleSamples = do.call("rbind",doubleSamples)

doubleTxCombine = apply(doubleSamples, 1, function(x) {
	cat(".")
	combineTxQuant(x[1], x[2], txPath)
})
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
pdMerge$RNum = sapply(pdMerge$RNum, "[[", 1)
pdMerge$BrNum = sapply(pdMerge$BrNum, "[[", 1)
pdMerge$Age = sapply(pdMerge$Age, "[[", 1)
pdMerge$Sex = sapply(pdMerge$Sex, "[[", 1)
pdMerge$Race = sapply(pdMerge$Race, "[[", 1)
pdMerge$Dx = sapply(pdMerge$Dx, "[[", 1)
pdMerge$PhaseII_HIPPO = sapply(pdMerge$PhaseII_HIPPO, "[[", 1)
pdMerge$RnumHIPPOPhaseII = sapply(pdMerge$RnumHIPPOPhaseII, "[[", 1)
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
	file = "count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n442.rda")
save(rse_exon, getRPKM, compress=TRUE,
	file = "count_data/hippo_brainseq_phase2_hg38_rseExon_merged_n442.rda")
save(rse_jxn, getRPM, compress=TRUE,
	file = "count_data/hippo_brainseq_phase2_hg38_rseJxn_merged_n442.rda")
save(rse_tx, compress=TRUE,
	file = "count_data/hippo_brainseq_phase2_hg38_rseTx_merged_n442.rda")

