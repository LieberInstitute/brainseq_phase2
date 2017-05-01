#!/bin/bash
#$ -cwd
#$ -N step3b-infer-strandness-DLPFC_polyA_degradation.hg38
#$ -o ./logs/infer-strandness-DLPFC_polyA_degradation.txt
#$ -e ./logs/infer-strandness-DLPFC_polyA_degradation.txt
#$ -hold_jid pipeline_setup,step3-hisat2-DLPFC_polyA_degradation.hg38
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Process the infer experiment info
Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step3b_infer_strandness.R

echo "**** Job ends ****"
date
