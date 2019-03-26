#!/bin/bash 
#$ -cwd
# Specify log file names
#$ -o logs/dlpfc_check_sczd_cell_prop.txt
#$ -e logs/dlpfc_check_sczd_cell_prop.txt
#$ -m e
#$ -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
#$ -N dlpfc_check_sczd_cell_prop

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript casectrl_DLPFC_allFeatures.R

echo "**** Job ends ****"
date
