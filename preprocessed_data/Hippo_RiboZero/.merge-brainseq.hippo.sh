#!/bin/bash
#$ -cwd
#$ -l mem_free=6G,h_vmem=10G,h_fsize=150G
#$ -N merge-brainseq.hippo
#$ -pe local 8
#$ -o ./logs/merge-brainseq.o.txt
#$ -e ./logs/merge-brainseq.e.txt
#$ -hold_jid pipeline_setup
#$ -m a

echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\n3de0c0d9f6fd9a855ca1cad7cba3b746f3ef6f4f"

Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step00-merge.R -s /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/SAMPLE_IDs.txt -o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/brainseq/hippo/merged_fastq -p TRUE -e fastq.gz -c 8

echo "**** Job ends ****"
date
