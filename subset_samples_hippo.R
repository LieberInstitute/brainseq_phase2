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
	

## drop DNAs: Br1779_1M, Br1178_1M, Br1142_h650, Br11836_Omni5M, 
	##			Br2407_1M, Br1357_1M, Br2173_1M

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

## subset 
txTpm = txTpm[,pd$SAMPLE_ID]

# txList = lapply(txFiles, read.delim, row.names=1, as.is=TRUE)
# tpmMat = sapply(txList, "[[", "TPM")
# rownames(tpmMat) = ss(rownames(txList[[1]]),"|", fixed=TRUE)

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
doubleTxCombine = mclapply(doubleSamples, combineTxQuant, salmonDir=txPath,mc.cores=4)
doubleTx = sapply(doubleTxCombine, "[[", "TPM")

singleTx = txTpm[,unlist(sIndexes[lengths(sIndexes) == 1])]
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
txMap = tx[rownames(txTpm)]
rse_tx = SummarizedExperiment(
	assays = list('tpm' = tpmMatMerge),
    colData = pdMerge, rowRanges = txMap)
	
################
## save ########
save(rse_gene, compress=TRUE,
	file = "count_data/hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda")
save(rse_exon, compress=TRUE,
	file = "count_data/hippo_brainseq_phase2_hg38_rseExon_merged_n447.rda")
save(rse_jxn, compress=TRUE,
	file = "count_data/hippo_brainseq_phase2_hg38_rseJxn_merged_n447.rda")
save(rse_tx, compress=TRUE,
	file = "count_data/hippo_brainseq_phase2_hg38_rseTx_merged_n447.rda")

