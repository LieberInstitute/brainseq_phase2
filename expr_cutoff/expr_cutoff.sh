#!/bin/bash
#$ -cwd
#$ -N expr_cutoff
#$ -o ./logs/expr_cutoff.txt
#$ -e ./logs/expr_cutoff.txt
#$ -m e
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript expr_cutoff.R

echo "**** Job ends ****"
date
