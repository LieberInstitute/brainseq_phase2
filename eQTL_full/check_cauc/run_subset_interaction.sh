#!/bin/bash
#$ -cwd
#$ -l mem_free=110G,h_vmem=110G,h_fsize=200G
#$ -N interaction_subset_cauc
#$ -o logs/interaction_subset_cauc.txt
#$ -e logs/interaction_subset_cauc.txt
#$ -m e
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"


Rscript results_subset_interaction.R

echo "**** Job ends ****"
date
