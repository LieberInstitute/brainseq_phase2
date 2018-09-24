#!/bin/bash
#$ -cwd
#$ -l mem_free=200G,h_vmem=240G,h_fsize=100G
#$ -N meanCoverage-brainseq.hippo
#$ -o ./logs/meanCoverage-brainseq.o.txt
#$ -e ./logs/meanCoverage-brainseq.e.txt
#$ -hold_jid pipeline_setup,coverage-brainseq.hippo
#$ -m a
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\n3de0c0d9f6fd9a855ca1cad7cba3b746f3ef6f4f"

## Locate normalized BigWig files and concatenate them in a space separated list
BIGWIGS=$(while read line; do ID=$(basename ${line}); echo "/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/${ID}.bw"; done < /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/SAMPLE_IDs.txt | paste -sd " ")

## Create mean of normalized bigwigs
module load wiggletools/default
module load ucsctools
wiggletools write /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/mean.wig mean ${BIGWIGS}
wigToBigWig /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/mean.wig /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/mean.bw

## Remove temp files
rm /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/mean.wig

echo "**** Job ends ****"
date
