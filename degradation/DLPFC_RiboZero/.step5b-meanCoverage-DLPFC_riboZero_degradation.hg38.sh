#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=120G,h_fsize=100G
#$ -N step5b-meanCoverage-DLPFC_riboZero_degradation.hg38
#$ -o ./logs/meanCoverage-DLPFC_riboZero_degradation.txt
#$ -e ./logs/meanCoverage-DLPFC_riboZero_degradation.txt
#$ -hold_jid pipeline_setup,step5-coverage-DLPFC_riboZero_degradation.hg38
#$ -m a
echo "**** Job starts ****"
date

if [ ! -f "inferred_strandness_pattern.txt" ]
then
    echo "Missing the file inferred_strandness_pattern.txt"
    exit 1
fi

## Load required software
module load wiggletools/default
module load ucsctools

## Read the strandness information
STRANDRULE=$(cat inferred_strandness_pattern.txt)

if [ ${STRANDRULE} == "none" ]
then
    ## Locate normalized BigWig files and concatenate them in a space separated list
    BIGWIGS=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/samples.manifest | awk '{print "/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Coverage/"$NF".bw"}' | paste -sd " ")
    
    ## Create mean of normalized bigwigs
    wiggletools write /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Coverage/mean.wig mean ${BIGWIGS}
    wigToBigWig /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Coverage/mean.wig /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Coverage/mean.bw
else
    for strand in Forward Reverse
    do
        echo "Processing strand ${strand}"
        ## Locate normalized BigWig files and concatenate them in a space separated list
        BIGWIGS=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/samples.manifest | awk -v strand="${strand}" '{print "/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Coverage/"$NF"."strand".bw"}' | paste -sd " ")
    
        ## Create mean of normalized bigwigs
        wiggletools write /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Coverage/mean.${strand}.wig mean ${BIGWIGS}
        wigToBigWig /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Coverage/mean.${strand}.wig /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Coverage/mean.${strand}.bw
    done
fi

## Remove temp files
rm /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Coverage/mean*.wig

echo "**** Job ends ****"
date
