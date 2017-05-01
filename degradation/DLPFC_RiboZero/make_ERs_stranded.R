###

library(rtracklayer)
library(derfinder)
library(jaffelab)

## chromosome info
chrInfo <- read.table('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg19.chrom.sizes.gencode',
    header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
chrInfo <- subset(chrInfo, chr %in% paste0('chr', c(1:22, 'X', 'Y')))

## load stranded mean bigwigs
meanBWs = paste0("Coverage/mean.",c("Forward","Reverse"), ".bw")
names(meanBWs) = c("Forward","Reverse")
meanList = lapply(meanBWs, import)
strand(meanList$Forward) = Rle("+")
strand(meanList$Reverse) = Rle("-")

## make ERs
meanListFilter = GRangesList(lapply(meanList, function(x) x[abs(x$score) >= 5]))
reduceList = endoapply(meanListFilter, reduce,min.gapwidth=2)

## at least 12 BP and on main chromosomes
erList = endoapply(reduceList, function(x) 
	x[width(x) >= 12 & seqnames(x) %in% chrInfo$chr])

## write out
export(erList$Forward, con = "bed/DLPFC_RiboZero_ERs_cut5_Forward.bed")
export(erList$Reverse, con = "bed/DLPFC_RiboZero_ERs_cut5_Reverse.bed")

#######################################
####### Stop and Run:			  #####
####### `qsub get_ER_coverage.sh` #####
#######################################

## phenotype
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/overall_degradation_pheno.rda")
pd = pd[pd$Dataset == "DLPFC_RiboZero",]

## read back in
## forward
bedForward = import("bed/DLPFC_RiboZero_ERs_cut5_Forward.bed")
fnForward = paste0("Coverage/cut5/",rownames(pd), ".Forward.cut5.txt")
names(fnForward) = rownames(pd)
erCovForward = sapply(fnForward, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f10", x)), 
		as.is=TRUE)$sum
})

## reverse
bedReverse = import("bed/DLPFC_RiboZero_ERs_cut5_Reverse.bed")
fnReverse = paste0("Coverage/cut5/",rownames(pd), ".Reverse.cut5.txt")
names(fnReverse) = rownames(pd)
erCovReverse = sapply(fnReverse, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f10", x)), 
		as.is=TRUE)$sum*-1 # flip sign
})

## join
regions = c(bedForward, bedReverse)
regionMat = rbind(degCovForward, degCovReverse)/100 # read length
rownames(regionMat) = names(regions) = paste0(seqnames(regions), ":", start(regions), 
								"-", end(regions), "(", strand(regions), ")")
regions$score = rowMeans(regionMat)

###############
## annotate ###
###############

library(biomaRt)
library(GenomicRanges)

## load genomic state
load('/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/gs/gs_gencode_v25_hg38.Rdata')
gs = gs_gencode_v25_hg38$fullGenome

all <- GRanges(names(seqlengths(gs)), IRanges(start = 1, end = seqlengths(gs)), strand = '*')
new_intergenic <- lapply(c('+', '-'), function(str) {
    onestrand <- gs[strand(gs) == str | strand(gs) == '*']
    onestrand_simple <- onestrand
    mcols(onestrand_simple) <- NULL
    pieces <- disjoin(c(all, onestrand_simple), ignore.strand = TRUE)
    ov <- countOverlaps(pieces, onestrand_simple, ignore.strand = TRUE)
    new_intergenic <- pieces[ov == 0]
    strand(new_intergenic) <- str
    new_intergenic$theRegion <- paste0('intergenic', str)
    new_intergenic$tx_id <- IntegerList(NA)
    new_intergenic$tx_name <- CharacterList(NA)
    new_intergenic$gene <- IntegerList(NA)
    return(new_intergenic)
})
gs_stranded <- c(gs, unlist(GRangesList(new_intergenic)))

ensemblAnno = annotateRegions(regions,gs_stranded,ignore.strand=FALSE,minoverlap=1)
countTable = ensemblAnno$countTable ## STRAND ISSUES, NEED TO FIX
countTable[,2] = rowSums(countTable[,2:4])
countTable = countTable[,c(1,2,5)]

## gene annotation
geneMap = read.delim("/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Counts/gene/DLPFC_Br2074_D4_ribo_Gencode.v25.hg38_Genes.counts",
	header=TRUE, as.is=TRUE,skip=1)[,1:6]
rownames(geneMap) = geneMap$Geneid
geneMap$Chr = ss(geneMap$Chr, ";")
geneMap$Start = as.numeric(ss(geneMap$Start, ";"))
tmp = strsplit(geneMap$End, ";")
geneMap$End = as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand = ss(geneMap$Strand, ";")
rownames(geneMap) = geneMap$Geneid
geneMap$gencodeID = geneMap$Geneid
geneMap$ensemblID = ss(geneMap$Geneid, "\\.")
geneMap$Geneid = NULL

## add symbol
ensembl = useMart("ENSEMBL_MART_ENSEMBL",  
		dataset="hsapiens_gene_ensembl", host="jul2016.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
		values=geneMap$ensemblID, mart=ensembl)
geneMap$Symbol = sym$hgnc_symbol[match(geneMap$ensemblID, sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene[match(geneMap$ensemblID, sym$ensembl_gene_id)]
geneMapGR = makeGRangesFromDataFrame(geneMap, keep=TRUE)

## overlaps
dA = distanceToNearest(regions, geneMapGR)
regions$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
regions$nearestID = names(geneMapGR)[subjectHits(dA)]
regions$distToGene = mcols(dA)$distance
mcols(regions) = cbind(mcols(regions), countTable)

## add additional annotation
regions$annoClass = NA
regions$annoClass[regions$exon > 0 & 
	regions$intron == 0 &
	regions$intergenic == 0] = "strictExonic"
regions$annoClass[regions$exon == 0 & 
	regions$intron > 0 &
	regions$intergenic == 0] = "strictIntronic"
regions$annoClass[regions$exon == 0 & 
	regions$intron == 0 &
	regions$intergenic > 0] = "strictIntergenic"
regions$annoClass[regions$exon > 0 & 
	regions$intron > 0 &
	regions$intergenic == 0] = "exonIntron"
regions$annoClass[regions$exon > 0 & 
	regions$intergenic > 0] = "extendUTR"

save(pd, regions, regionMat, file = "expressedRegions_DLPFC_RiboZero_degradation_cut5_hg38_n20.rda")

##########
## compare to derfinder code
# bw =  paste0("/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Coverage/",
	# rep(rownames(pd), each=2), ".", rep(c("Forward","Reverse"), times=12), ".bw")
bwForward =  paste0("/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Coverage/",
	rownames(pd), ".Forward.bw")
names(bwForward) = rownames(pd)
fullCovForward = fullCoverage(bwForward, chrs = paste0("chr", 1:22), mc.cores=4)

bwReverse =  paste0("/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Coverage/",
	rownames(pd), ".Reverse.bw")
names(bwReverse) = rownames(pd)
fullCovReverse = fullCoverage(bwReverse, chrs = paste0("chr", 1:22), mc.cores=4)

save(fullCovForward, fullCovReverse, pd, 
	file = "fullCoverage_bothStrands_DLPFC_RiboZero.rda")

regionMatForward = regionMatrix(fullCovForward, cutoff = 5, L = 100, 
	returnBP=FALSE, mc.cores=6)
regionMatReverse = regionMatrix(fullCovReverse, cutoff = 5, L = 100, 
	returnBP=FALSE, mc.cores=6)
save(regionMatForward, regionMatReverse, pd, 
	file = "regionMats_viaDerfinder_bothStrands_DLPFC_RiboZero.rda")
