#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=40G,h_vmem=40G,h_fsize=100G
#$ -N span_limma_reg_specific_jxn_adult
#$ -o ./logs/span_limma_reg_specific_jxn_adult.txt
#$ -e ./logs/span_limma_reg_specific_jxn_adult.txt
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
Rscript span_limma_reg_specific.R -t jxn -a adult

echo "**** Job ends ****"
date
