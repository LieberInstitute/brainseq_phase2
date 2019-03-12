#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N limma_reg_specific_gene_adult_subsample
#$ -o ./subsample/logs/limma_reg_specific_gene_adult_subsample.$TASK_ID.txt
#$ -e ./subsample/logs/limma_reg_specific_gene_adult_subsample.$TASK_ID.txt
#$ -t 1-100
#$ -m a

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

# module load conda_R/3.4.x
Rscript limma_reg_specific_subsample.R -t gene -a adult -i ${SGE_TASK_ID}

echo "**** Job ends ****"
date
