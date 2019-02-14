#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=160G,h_vmem=160G,h_fsize=200G
#$ -N bsI_raggr_dlpfc_eQTL
#$ -o ./eqtl_tables/logs/eQTL_dlpfc_4_features.txt
#$ -e ./eqtl_tables/logs/eQTL_dlpfc_4_features.txt
#$ -m e
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

## Pre-reqs:
# mkdir -p eqtl_tables
# mkdir -p eqtl_tables/logs
# mkdir -p eqtl_tables/rdas

Rscript raggr_run_eqtls_dlpfc.R

echo "**** Job ends ****"
date
