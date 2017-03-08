###

library(jaffelab)
library(readxl)

## DLPFC
load("/users/ajaffe/Lieber/Projects/RNAseq/Degradation/rdas/annotated_phenotype_data_degradation_twoLibraries.rda")
pdDlpfc = pd[,1:15] # drop old alignmnet states
pdDlpfc$Region = "DLPFC"
rm(pd)

pdDlpfc$SampleID = paste0("DLPFC_", pdDlpfc$SampleID)
rownames(pdDlpfc) = pdDlpfc$SampleID

###############
## hippo

## Read in data
pdHippo = read_excel("/users/ajaffe/Lieber/Projects/RnaQualPaper/Degradation/Hippo/RNA-Seq Aliquots Hippo degradation 20161109.xlsx",
	skip=1)

# change column names
colnames(pdHippo)[c(1:2,7,14)] = c("Position","ID","DegradationTime", "Dx")
pdHippo = pdHippo[!is.na(pdHippo$ID),]
colnames(pdHippo)[3] = gsub("#", "Num", colnames(pdHippo)[3])

pdHippo$SampleID = paste0("HIPPO_", pdHippo$BrNum, "_",
	pdHippo$Position, "_ribo")
rownames(pdHippo) = pdHippo$SampleID
pdHippo$Flowcell = "AHC7KYBBXX"
pdHippo$LibraryProtocol = "RiboZero"

### get reads
thePath = "/dcl01/lieber/ajaffe/Nina/161116_J00144_0062_AHC7KYBBXX/degradation/"
reads = list.files(thePath, pattern = "*fastq.gz$")
names(reads) = rownames(pdHippo)[match(ss(reads, "_"), pdHippo$ID )]

readList = split(reads, names(reads))
readMat = do.call("rbind", readList)

mm = match(rownames(pdHippo), rownames(readMat))
pdHippo = pdHippo[!is.na(mm),]
readMat = readMat[mm[!is.na(mm)],]
pdHippo$leftRead = paste0(thePath, readMat[,1])
pdHippo$rightRead = paste0(thePath, readMat[,2])

## keep consistent columns
n = intersect(names(pdHippo), names(pdDlpfc))
pdDlpfc = pdDlpfc[,n]
pdHippo = pdHippo[,n]
pd = rbind(pdDlpfc, pdHippo)
pd$Dataset = paste0(pd$Region, "_", pd$LibraryProtocol)
save(pd, file="overall_degradation_pheno.rda")

#########################
## write output for pipeline
left = strsplit(pd$leftRead, ",")
right = strsplit(pd$rightRead, ",")
names(left) = names(right) = rownames(pd)
left = unlist(left)
right = unlist(right)

sample.manifest = data.frame(leftRead = left, 
	md5 = 0, rightRead = right, md5=0,
	lab = names(left),stringsAsFactors=FALSE)
sample.manifest$lab = gsub("ribo[0-9]", "ribo", sample.manifest$lab)
sample.manifest$dataset = paste0(ss(sample.manifest$lab,"_"),
					"_",ss(sample.manifest$lab,"_",4))
sampleList = split(sample.manifest,sample.manifest$dataset)

## write out
write.table(sampleList$DLPFC_poly[,1:5], 
	file = "DLPFC_polyA/samples.manifest", sep="\t",
	col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(sampleList$DLPFC_ribo[,1:5], 
	file = "DLPFC_RiboZero/samples.manifest", sep="\t",
	col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(sampleList$HIPPO_ribo[,1:5], 
	file = "Hippo_RiboZero/samples.manifest", sep="\t",
	col.names=FALSE, row.names=FALSE, quote=FALSE)