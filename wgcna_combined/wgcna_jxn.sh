#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_fsize=300G
#$ -N wgcna_jxn
#$ -o logs/wgcna_jxn.txt
#$ -e logs/wgcna_jxn.txt
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

module load conda_R/3.4.x
Rscript wgcna_jxn.R

echo "**** Job ends ****"
date
