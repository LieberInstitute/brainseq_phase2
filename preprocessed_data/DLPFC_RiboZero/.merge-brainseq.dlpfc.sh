#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=6G,h_vmem=10G,h_fsize=150G
#$ -N merge-brainseq.dlpfc
#$ -pe local 8
#$ -o ./logs/merge-brainseq.o.txt
#$ -e ./logs/merge-brainseq.e.txt
#$ -hold_jid pipeline_setup
#$ -m a

echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\na2b1cda9a67211320865c002ee1a20429cce32b1"

Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step00-merge.R -s /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/SAMPLE_IDs.txt -o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/brainseq/dlpfc/merged_fastq -p TRUE -e fastq.gz -c 8

echo "**** Job ends ****"
date
