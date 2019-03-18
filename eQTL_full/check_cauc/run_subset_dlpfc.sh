#!/bin/bash
#$ -cwd
#$ -l mem_free=150G,h_vmem=150G,h_fsize=200G
#$ -N dlpfc_subset_cauc
#$ -o logs/dlpfc_subset_cauc.txt
#$ -e logs/dlpfc_subset_cauc.txt
#$ -m e
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript results_subset_dlpfc.R

echo "**** Job ends ****"
date
