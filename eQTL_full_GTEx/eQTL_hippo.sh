#!/bin/bash
#$ -cwd
#$ -l mem_free=175G,h_vmem=175G,h_fsize=200G
#$ -N bsII_hippo_eQTL
#$ -o ./eqtl_tables/logs/eQTL_hippo_4_features.txt
#$ -e ./eqtl_tables/logs/eQTL_hippo_4_features.txt
#$ -m e
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript run_eqtls_hippo.R

echo "**** Job ends ****"
date
