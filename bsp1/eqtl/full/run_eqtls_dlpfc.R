####

### libraries
library('SummarizedExperiment')
library('jaffelab')
library('MatrixEQTL')
library('sva')
library('sessioninfo')

######################
### load data ####
######################

load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/bsp1_gene.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/bsp1_exon.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/bsp1_jxn.Rdata", verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/bsp1/data/bsp1_tx.Rdata", verbose = TRUE)

# # fix junction row names
# rownames(rse_jxn) = paste0(seqnames(rse_jxn),":",start(rse_jxn),"-",end(rse_jxn),"(",strand(rse_jxn),")")

# # sum totalMapped IntegerLists (so getRPKM works later)
# colData(rse_tx)$totalMapped =
	# colData(rse_jxn)$totalMapped =
	# colData(rse_exon)$totalMapped = 
	# colData(rse_gene)$totalMapped  = sapply(colData(rse_gene)$totalMapped, sum)

## keep adult samples & correct region
keepInd = which(colData(bsp1_gene)$Age > 13 & colData(bsp1_gene)$Region == "DLPFC")
rse_gene = bsp1_gene[,keepInd]
rse_exon = bsp1_exon[,keepInd]
rse_jxn = bsp1_jxn[,keepInd]
rse_tx = bsp1_tx[,keepInd]

## Drop unused objects and print some info
rm(bsp1_gene, bsp1_exon, bsp1_jxn, bsp1_tx)
table(rse_gene$Dx)


if(!file.exists('eqtl_tables/matrixEqtl_output_dlpfc_4features.rda')) {
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
    load('/dcl01/ajaffe/data/lab/brainseq_phase1/genotype_data/brainseq_phase1_Genotypes_n732.rda', verbose = TRUE)

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
    	levels = c("Control", "Schizo", "MDD", "Bipolar"))

    mod = model.matrix(~Dx + Sex + as.matrix(mds[,1:5]),
    	data = pd)
    colnames(mod)[grep('mds', colnames(mod))] = colnames(mds)[1:5]


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
    	file="eqtl_tables/rdas/pcs_dlpfc_4features_filtered_over13.rda")
    # load("eqtl_tables/rdas/pcs_dlpfc_4features_filtered_over13.rda")

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
    	file="eqtl_tables/matrixEqtl_output_dlpfc_4features.rda")
} else {
    rm(keepInd)
    message(paste(Sys.time(), 'loading eqtl_tables/matrixEqtl_output_dlpfc_4features.rda'))
    load("eqtl_tables/matrixEqtl_output_dlpfc_4features.rda", verbose = TRUE)
}

######################
###### annotate ######

if(!file.exists('eqtl_tables/mergedEqtl_output_dlpfc_4features_in_progress.rda')) {
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

    rm(meGene, meExon, meJxn, meTx)
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

    ## Process each one individually to reduce the memory blow up
    rm(rse_tx)
    geneEqtl$gencodeTx <- CharacterList(
        as.list(rowRanges(rse_gene)$gencodeTx[match(geneEqtl$gene, rownames(rse_gene))])
    )
    rm(rse_gene)
    exonEqtl$gencodeTx <- CharacterList(
        as.list(rowRanges(rse_exon)$gencodeTx[match(exonEqtl$gene, rownames(rse_exon))])
    )
    rm(rse_exon)
    jxnEqtl$gencodeTx <- CharacterList(
        as.list(rowRanges(rse_jxn)$gencodeTx[match(jxnEqtl$gene, rownames(rse_jxn))])
    )
    rm(rse_jxn)
    txEqtl$gencodeTx <- CharacterList(
        as.list(txEqtl$gene)
    )

    save(geneEqtl, exonEqtl, jxnEqtl, txEqtl, file = 'eqtl_tables/mergedEqtl_output_dlpfc_4features_in_progress.rda')
} else {
    load('eqtl_tables/mergedEqtl_output_dlpfc_4features_in_progress.rda', verbose = TRUE)
}

# merge
allEqtl = rbind(geneEqtl, exonEqtl, jxnEqtl, txEqtl)
save(allEqtl, file="eqtl_tables/mergedEqtl_output_dlpfc_4features.rda",compress=TRUE)


# #############
# # metrics ###
# sigEqtl = allEqtl[allEqtl$FDR < 0.01,] # fdr significant
# length(unique(sigEqtl$EnsemblGeneID))

# sigEqtlList = split(sigEqtl, factor(sigEqtl$Type,
	# levels=c("Gene","Exon","Jxn", "Tx")))

# sapply(sigEqtlList, function(x) max(x$pvalue))

# sapply(sigEqtlList, function(x) length(unique(x$EnsemblGeneID)))
# sapply(sigEqtlList, function(x) length(unique(x$Symbol)))
# sapply(sigEqtlList, function(x) length(unique(x$snps)))

# sapply(sigEqtlList, function(x) quantile(abs(x$beta)))

# sapply(sigEqtlList, function(x) table(x$Class[!duplicated(x$gene)]))
# sapply(sigEqtlList, function(x) prop.table(table(x$Class)))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
