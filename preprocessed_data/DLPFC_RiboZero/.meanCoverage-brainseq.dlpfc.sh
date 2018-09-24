#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=200G,h_vmem=240G,h_fsize=100G
#$ -N meanCoverage-brainseq.dlpfc
#$ -o ./logs/meanCoverage-brainseq.o.txt
#$ -e ./logs/meanCoverage-brainseq.e.txt
#$ -hold_jid pipeline_setup,coverage-brainseq.dlpfc
#$ -m a
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\na2b1cda9a67211320865c002ee1a20429cce32b1"

## Locate normalized BigWig files and concatenate them in a space separated list
BIGWIGS=$(while read line; do ID=$(basename ${line}); echo "/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Coverage/${ID}.bw"; done < /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/SAMPLE_IDs.txt | paste -sd " ")

## Create mean of normalized bigwigs
module load wiggletools/default
module load ucsctools
wiggletools write /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Coverage/mean.wig mean ${BIGWIGS}
wigToBigWig /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Coverage/mean.wig /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Coverage/mean.bw

## Remove temp files
rm /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Coverage/mean.wig

echo "**** Job ends ****"
date
