#!/bin/bash
#$ -cwd
#$ -l mem_free=45G,h_vmem=45G,h_fsize=100G
#$ -N limma_casecontrolint_tx
#$ -o ./logs/limma_casecontrolint_tx.txt
#$ -e ./logs/limma_casecontrolint_tx.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

#module load conda_R/3.4.x
Rscript limma_casecontrolint.R -t tx

echo "**** Job ends ****"
date
