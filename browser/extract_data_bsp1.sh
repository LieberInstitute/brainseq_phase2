#!/bin/bash
#$ -cwd
#$ -l mem_free=330G,h_vmem=330G,h_fsize=200G
#$ -N extract_data_bsp1
#$ -o logs/extract_data_bsp1.txt
#$ -e logs/extract_data_bsp1.txt
#$ -m e
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"


Rscript extract_data_bsp1.R

echo "**** Job ends ****"
date
