#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=200G,h_vmem=240G,h_fsize=100G
#$ -N step5b-meanCoverage-HIPPO_RiboZero_BrainSeq_Phase2.LIBD
#$ -o ./logs/meanCoverage-HIPPO_RiboZero_BrainSeq_Phase2.txt
#$ -e ./logs/meanCoverage-HIPPO_RiboZero_BrainSeq_Phase2.txt
#$ -hold_jid pipeline_setup,step5-coverage-HIPPO_RiboZero_BrainSeq_Phase2.LIBD
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

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
    BIGWIGS=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/samples.manifest | awk '{print "/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/"$NF".bw"}' | paste -sd " ")
    
    ## Create mean of normalized bigwigs
    wiggletools write /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/mean.wig mean ${BIGWIGS}
    wigToBigWig /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/mean.wig /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/mean.bw
else
    for strand in Forward Reverse
    do
        echo "Processing strand ${strand}"
        ## Locate normalized BigWig files and concatenate them in a space separated list
        BIGWIGS=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/samples.manifest | awk -v strand="${strand}" '{print "/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/"$NF"."strand".bw"}' | paste -sd " ")
    
        ## Create mean of normalized bigwigs
        wiggletools write /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/mean.${strand}.wig mean ${BIGWIGS}
        wigToBigWig /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/mean.${strand}.wig /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/mean.${strand}.bw
    done
fi

## Remove temp files
rm /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/mean*.wig

echo "**** Job ends ****"
date
