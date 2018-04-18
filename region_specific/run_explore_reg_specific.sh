#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
#$ -N explore_reg_specific
#$ -o ./logs/explore_reg_specific.txt
#$ -e ./logs/explore_reg_specific.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript explore_reg_specific.R

echo "**** Job ends ****"
date
