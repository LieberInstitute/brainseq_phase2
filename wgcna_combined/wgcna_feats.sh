#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=200G
#$ -N wgcna_dlpfc_feats
#$ -o logs/wgcna_dlpfc_feats.txt
#$ -e logs/wgcna_dlpfc_feats.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript wgcna_exon_jxn.R

echo "**** Job ends ****"
date

