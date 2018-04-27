#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=60G,h_vmem=60G,h_fsize=200G
#$ -N wgcna_feats
#$ -o logs/wgcna_feats.txt
#$ -e logs/wgcna_feats.txt
#$ -pe local 4
#$ -m e
#$ -hold_jid get_ind_exons_jxns

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

