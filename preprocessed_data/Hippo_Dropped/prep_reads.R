##
library(readxl)
library(jaffelab)

## pheno
pd = read_excel("Hippo_Consortium2_omitted_list.xlsx")
colnames(pd)[1] = "BrNum"
pd$BrNum = paste0("Br", pd$BrNum)
pd$RNum = paste0("R", pd$RNum)

## get reads
fastqPath = "/dcl01/lieber/ajaffe/Nina/BrainSeq_PhaseII/Hippo_RiboZero/Hippo_omitted_samples"
fastqFiles = list.files(fastqPath, full.names=TRUE, pattern = "fastq.gz$", recursive=TRUE)
names(fastqFiles) = ss(list.files(fastqPath, pattern = "fastq.gz$", recursive=TRUE), "/", 2)

## rename as RNum_Flow_Lane
names(fastqFiles) = paste0(ss(names(fastqFiles), "_"), "_",
	ss(names(fastqFiles), "_",2), "_", ss(names(fastqFiles), "_",4))
	
fastqList = split(fastqFiles, names(fastqFiles))
fqMat = do.call("rbind", fastqList)

dat = data.frame(leftRead = fqMat[,1],
	leftMd5 = 0, rightRead = fqMat[,2],
	SampleID = ss(rownames(fqMat),"_L0"), 
	stringsAsFactors=FALSE)
write.table(dat, file = "samples.manifest",
	sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)