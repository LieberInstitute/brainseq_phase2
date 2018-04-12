#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
#$ -N explore_limma_dev
#$ -o ./logs/explore_limma_dev.txt
#$ -e ./logs/explore_limma_dev.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript explore_limma_dev.R

echo "**** Job ends ****"
date
