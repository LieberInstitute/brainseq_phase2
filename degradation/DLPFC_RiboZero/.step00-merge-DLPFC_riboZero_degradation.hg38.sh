#!/bin/bash
#$ -cwd
#$ -l mem_free=3G,h_vmem=5G,h_fsize=150G
#$ -N step00-merge-DLPFC_riboZero_degradation.hg38
#$ -pe local 8
#$ -o ./logs/merge-DLPFC_riboZero_degradation.txt
#$ -e ./logs/merge-DLPFC_riboZero_degradation.txt
#$ -hold_jid pipeline_setup
#$ -m a

echo "**** Job starts ****"
date

Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step00-merge.R -s /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/.samples_unmerged.manifest -o /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/merged_fastq -c 8

echo "**** Job ends ****"
date
