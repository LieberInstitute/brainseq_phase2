#!/bin/bash
#$ -cwd
#$ -l mem_free=250G,h_vmem=250G,h_fsize=200G
#$ -N dlpfc_subset
#$ -o logs/dlpfc_subset.txt
#$ -e logs/dlpfc_subset.txt
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/results_subset_dlpfc.R

echo "**** Job ends ****"
date
