#!/bin/bash
#$ -cwd
#$ -l mem_free=440G,h_vmem=440G,h_fsize=200G
#$ -N dlpfc_subset
#$ -o logs/dlpfc_subset.txt
#$ -e logs/dlpfc_subset.txt
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
