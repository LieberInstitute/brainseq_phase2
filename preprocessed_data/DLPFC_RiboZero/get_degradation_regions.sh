#!/bin/bash
#$ -cwd
#$ -N degradation_DLPFC
#$ -o ./logs/degradation_DLPFC.txt
#$ -e ./logs/degradation_DLPFC.txt
#$ -m e
#$ -pe local 8
#$ -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=100G

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript get_degradation_regions.R

echo "**** Job ends ****"
date
