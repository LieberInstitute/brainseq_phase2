#!/bin/bash
#$ -cwd
#$ -l mem_free=3G,h_vmem=5G,h_fsize=150G
#$ -N step00-merge-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples
#$ -pe local 8
#$ -o ./logs/merge-HIPPO_riboZero_BrainSeq_Phase2.txt
#$ -e ./logs/merge-HIPPO_riboZero_BrainSeq_Phase2.txt
#$ -hold_jid pipeline_setup
#$ -m a

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step00-merge.R -s /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/.samples_unmerged.manifest -o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/merged_fastq -c 8

echo "**** Job ends ****"
date
