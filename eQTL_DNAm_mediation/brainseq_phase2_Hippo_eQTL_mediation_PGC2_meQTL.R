#qsub -l bluejay,mf=150G,h_vmem=150G,h_fsize=200G,h_stack=256M -cwd -b y -M stephensemick@gmail.com -o log -e log R CMD BATCH --no-save 03a_meQTL_Hippo.R

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(minfi)
library(GenomicRanges)

######################
### load data ####
######################
load('/dcl01/lieber/ajaffe/Steve/Hippo_meQTL/rdas/cleanSamples_n694_processed_data_postfiltered.rda')

## keep adult samples - keep both regions
keepInd = which(pd$Age > 13 & pd$`Brain.Region`=="Hippo")
pd = pd[keepInd, ]
meth = as.matrix(bVals[,keepInd])

##### Subset to probes available
load('/dcl01/lieber/ajaffe/Steve/meth450k_annotation_hg38/hg38_out/rdas/hg38_goldset_annotation.rda') #load hg38 position annotation
mp = rownames(meth)[rownames(meth)%in%goldset$Name]
table(rownames(meth)%in%goldset$Name)
meth=meth[mp,] #drop probes that don't map to hg38
goldset = goldset[match(mp,goldset$Name),]

print("....data loaded....")

######################
### snp data #########
######################

## load SNP data
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551.rda')
### make mds and snp dimensions equal to N
###(repeat rows or columns for BrNum replicates)
mds = mds[pd$BrNum,]
snp = snp[,pd$BrNum]
rownames(mds) = colnames(snp) = pd$Chip

## drop SNPs not mapping to hg38
keepIndex = which(!is.na(snpMap$chr_hg38))
snpMap = snpMap[keepIndex,]
snp = snp[keepIndex,]

## subset to PGC2 SNPs of interest 
load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_DNAm_mediation/PGC2_theSnps_full.rda')
keep_PGC_SNP_index = which(paste0('chr',snpMap$CHR,':',snpMap$POS) %in% theSnps$chrpos)
snpMap = snpMap[keep_PGC_SNP_index,]
snp = snp[keep_PGC_SNP_index,]

length(keep_PGC_SNP_index)
#[1] 89469

######################
# statistical model ##
######################
pd$Dx = factor(pd$Dx, levels = c("Control", "Schizo"))
pd$Sex <- as.factor(pd$Sex)

mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]), data = pd)
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

#pcaMeth = prcomp(t(meth))
#kMeth = num.sv(meth, mod)
#kMeth = min(kMeth, 25)
#methPCs = pcaMeth$x[,1:kMeth]
#
#save(methPCs, file="rdas/pcs_methPCs_regions_filtered_over13_Hippo_only.rda")
load("/dcl01/lieber/ajaffe/Steve/Hippo_meQTL/rdas/pcs_methPCs_regions_filtered_over13_Hippo_only.rda")

print("....pcas created....")

######################
# create covariate objects #
######################
covs = SlicedData$new(t(cbind(mod[,-1], methPCs ) ) )

print("....covariate object created....")

#######################
# create meth objects #
#######################
		
#### make meth objects
meth = SlicedData$new(meth)
meth$ResliceCombined(sliceSize = 5000)

### meth position
posMeth = data.frame(probeid = goldset$Name,
	chr = goldset$chr_hg38, start = goldset$predictedPos, end=goldset$predictedPos)

#######################
# Run matrix eQTL #
#######################

print("....beginning eQTL analysis....")
	
meMeth = Matrix_eQTL_main(snps=theSnps, 
						 gene = meth, 
						 cvrt = covs, 
						 output_file_name.cis =  "HippoMeth.ctxt" ,
						 pvOutputThreshold.cis = 1, 
						 pvOutputThreshold=0,
						 snpspos = snpspos, 
						 genepos = posMeth, 
						 useModel = modelLINEAR,
						 cisDist=2e4,
						 pvalue.hist = 100,
						 min.pv.by.genesnp = TRUE)
save(meMeth, file="PGC2_meqtl/PGC2_cis_raw_output_hippo_methylation_adultsOnly.rda")

#######
## annotate
#load('/dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_DNAm_mediation/PGC2_meqtl/PGC2_cis_raw_output_hippo_methylation_adultsOnly.rda')
### filter
meqtl = meMeth$cis$eqtl
meqtl$snps = as.character(meqtl$snps)
colnames(meqtl)[1:2] = c("SNP","cpg")
meqtl$cpg = as.character(meqtl$cpg)
#
### annotate
m = match(meqtl$SNP, snpMap$SNP)
meqtl$snpChr = snpMap$chr_hg38[m]
meqtl$snpPos = snpMap$pos_hg38[m]
meqtl$snpRsNum = snpMap$name[m]
meqtl$snpCounted = snpMap$newCount[m]
meqtl$snpAlt = snpMap$newRef[m]
meqtl$snpType=snpMap$Type[m]
meqtl$inSampleMAF = rowSums(snp[m,],na.rm=TRUE)/ (2*rowSums(!is.na(snp[m,])))
#	
#### see if SNP distrupts a CpG
library(BSgenome.Hsapiens.UCSC.hg38)
meqtl$snpChr = gsub("chr23","chrX", meqtl$snpChr)
gr = GRanges(meqtl$snpChr, IRanges(meqtl$snpPos-1, meqtl$snpPos+1))
trio = getSeq(Hsapiens, gr)
meqtl$disruptCpG = vcountPattern("CG", trio)
#
cpgIndex=match(meqtl$cpg, goldset$Name)
meqtl$methChr = as.character(goldset[cpgIndex,'chr_hg38'])
meqtl$methPos = as.numeric(goldset[cpgIndex,'predictedPos'])
meqtl$methProbeRs = as.character(goldset[cpgIndex,'Probe_rs'])
meqtl$methCpgRs = as.character(goldset[cpgIndex,'CpG_rs'])
meqtl$methSbeRs = as.character(goldset[cpgIndex,'SBE_rs'])
meqtl$methRelationToIsland = as.character(goldset[cpgIndex,'Relation_to_Island'])
#
meqtl$distMethToSnp = meqtl$methPos - meqtl$snpPos
meqtl$meanMeth = rowMeans(meth)[meqtl$cpg]
#
save(meqtl, file='PGC2_meqtl/PGC2_cis_annotated_output_hippo_methylation_adultsOnly.rda')

#### take best per meqtl
#sig = meqtl[!duplicated(meqtl$cpg),]
#
#tt = table(meqtl$cpg)
#sig$numSnps = tt[sig$cpg]
#dis = sapply(split(meqtl$disruptCpG, meqtl$cpg),sum)
#sig$disruptCpG = dis[sig$cpg]
#
#snprange = t(sapply(split(meqtl$snpPos, meqtl$cpg), range))
#snprange= snprange[sig$cpg,]
#sig$startSnpLD = snprange[,1]
#sig$endSnpLD = snprange[,2]
#sig$snpLDLength = snprange[,2]-snprange[,1]
#
#### add gene info
#library(bumphunter)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#theTranscripts = annotateTranscripts(
#	TxDb.Hsapiens.UCSC.hg38.knownGene,codingOnly=TRUE)
#
#an = annotateNearest(map, theTranscripts)	
#map$nearestGene = as.character(theTranscripts$Gene)[an$subjectHits]
#map$nearestGeneDist = an$dist

#sig$nearestGene = map$nearestGene[match(sig$cpg, names(map))]	
#sig$nearestGeneDist = map$nearestGeneDist[match(sig$cpg, names(map))]	
#meqtl$nearestGene = map$nearestGene[match(meqtl$cpg, names(map))]	
#meqtl$nearestGeneDist = map$nearestGeneDist[match(meqtl$cpg, names(map))]	
#
##### add fetal effect
#fIndex=which(pd$Age < 0)
#pdFetal= pd[fIndex,]
#pcaFetal = prcomp(t(p[,fIndex]))
#modFetal = model.matrix(~ pdFetal$snpPC1+
#	pdFetal$snpPC2 + pdFetal$snpPC3)
#nsvFetal = num.sv(p[,fIndex], modFetal)
#
#snpFetal = as.matrix(snp[match(sig$snps,rownames(snp)),fIndex])
#pFetal = as.matrix(p[sig$cpg,fIndex])
#outFetal = matrix(nrow = nrow(sig), nc = 3)
#for(j in 1:nrow(sig)) {
#	if(j %% 5000 == 0) cat(".")
#	outFetal[j,] = summary(lm(pFetal[j,]~snpFetal[j,] + 
#		cbind(modFetal[,-1], pcaFetal$x[,1:nsvFetal])))$coef[2,c(1,3,4)]
#}
#colnames(outFetal) = paste0("fetal_", c("slope","tstat","pval"))
#sig = cbind(sig, outFetal)
#sig$fetal_MAF = rowSums(snpFetal, na.rm=TRUE)/(2*rowSums(!is.na(snpFetal)))
#
#### add adult effect
#aIndex=which(pd$Age >13)
#pdAdult= pd[aIndex,]
#pcaAdult = prcomp(t(p[,aIndex]))
#modAdult = model.matrix(~ pdAdult$snpPC1+
#	pdAdult$snpPC2 + pdAdult$snpPC3)
#nsvAdult = num.sv(p[,aIndex], modAdult)
#
#snpAdult = as.matrix(snp[match(sig$snps,rownames(snp)),aIndex])
#pAdult = as.matrix(p[sig$cpg,aIndex])
#outAdult = matrix(nrow = nrow(sig), nc = 3)
#for(j in 1:nrow(sig)) {
#	if(j %% 5000 == 0) cat(".")
#	outAdult[j,] = summary(lm(pAdult[j,]~snpAdult[j,] + 
#		cbind(modAdult[,-1], pcaAdult$x[,1:nsvAdult])))$coef[2,c(1,3,4)]
#}
#colnames(outAdult) = paste0("adult_", c("slope","tstat","pval"))
#sig = cbind(sig, outAdult)
#sig$adult_MAF = rowSums(snpAdult, na.rm=TRUE)/(2*rowSums(!is.na(snpAdult)))
#
#
### race MAF
#rIndexes=split(seq(along=pd2$Race), pd2$Race)
#mafRace = sapply(rIndexes, function(ii){
#	ssnp = as.matrix(snp2[sig$snps,ii])
#	rowSums(ssnp, na.rm=TRUE)/(2*rowSums(!is.na(ssnp)))
#})
#colnames(mafRace) = paste0(colnames(mafRace),"_inSampleMAF")
#sig = cbind(sig, mafRace)
#	
#save(sig, meqtl, compress=TRUE)
#