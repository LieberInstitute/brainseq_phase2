#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N ligases_bsp2
#$ -o logs/ligases_bsp2.$TASK_ID.txt
#$ -e logs/ligases_bsp2.$TASK_ID.txt
#$ -m a
#$ -t 1-316
#$ -tc 20

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
SYMBOL=$(awk '{ print $0}' ligases_found.txt | awk "NR==${SGE_TASK_ID}")
Rscript age_custom_plot.R -t gene -s ${SYMBOL} -o ligases_pdf

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.0
## available from http://research.libd.org/sgejobs/

## Manual ones Jim requested:
# Rscript age_custom_plot.R -t gene -s CRBN -o ligases_pdf
# Rscript age_custom_plot.R -t gene -s VHL -o ligases_pdf
