#!/bin/bash
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=200G
#$ -N hippo_subset_cauc
#$ -o logs/hippo_subset_cauc.txt
#$ -e logs/hippo_subset_cauc.txt
#$ -m e
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript results_subset_hippo.R

echo "**** Job ends ****"
date
