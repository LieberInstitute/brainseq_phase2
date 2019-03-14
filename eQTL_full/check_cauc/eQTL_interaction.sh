#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=60G,h_vmem=60G,h_fsize=200G
#$ -N inter_eQTL_bsII_CAUC
#$ -o eqtl_tables/logs/eQTL_interaction_gene.txt
#$ -e eqtl_tables/logs/eQTL_interaction_gene.txt
#$ -m e

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

echo "**** Job starts ****"
date

Rscript run_eqtls_interaction.R

echo "**** Job ends ****"
date
