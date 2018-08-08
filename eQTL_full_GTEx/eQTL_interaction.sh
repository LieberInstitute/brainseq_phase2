#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=200G,h_vmem=200G,h_fsize=200G
#$ -N inter_eQTL_bsII
#$ -o eqtl_tables/logs/eQTL_interaction_4_features2.txt
#$ -e eqtl_tables/logs/eQTL_interaction_4_features2.txt
#$ -m e
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript run_eqtls_interaction.R

echo "**** Job ends ****"
date
