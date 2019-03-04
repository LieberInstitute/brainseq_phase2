#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -N chmod
#$ -o ./logs/chmod.txt
#$ -e ./logs/cmod.txt
#$ -m e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

chmod -R 770 /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas

echo "**** Job ends ****"
date
