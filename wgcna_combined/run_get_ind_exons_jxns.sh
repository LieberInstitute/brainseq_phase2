#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=15G,h_vmem=15G,h_fsize=100G -pe local 8
#$ -N get_ind_exons_jxns
#$ -o ./logs/get_ind_exons_jxns.txt
#$ -e ./logs/get_ind_exons_jxns.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript get_ind_exons_jxns.R

echo "**** Job ends ****"
date
