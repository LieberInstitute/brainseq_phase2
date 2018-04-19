#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -N explore_reg_specific_jxn_adult_top
#$ -o ./logs/explore_reg_specific_jxn_adult_top.txt
#$ -e ./logs/explore_reg_specific_jxn_adult_top.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript explore_reg_specific_top.R -t jxn -a adult

echo "**** Job ends ****"
date
