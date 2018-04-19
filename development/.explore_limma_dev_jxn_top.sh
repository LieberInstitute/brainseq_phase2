#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=60G,h_vmem=60G,h_fsize=100G
#$ -N explore_limma_dev_jxn_top
#$ -o ./logs/explore_limma_dev_jxn_top.txt
#$ -e ./logs/explore_limma_dev_jxn_top.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

module load conda_R/3.4.x
Rscript explore_limma_dev_top.R -t jxn

echo "**** Job ends ****"
date
