#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=6G,h_fsize=100G
#$ -N span_limma_reg_specific_jxn_adult
#$ -o ./logs/span_limma_reg_specific_jxn_adult.txt
#$ -e ./logs/span_limma_reg_specific_jxn_adult.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: lcollado"
echo "Job id: "
echo "Job name: "
echo "Hostname: compute-066"
echo "Task id: "

module load conda_R/3.4.x
Rscript span_limma_reg_specific.R -t jxn -a adult

echo "**** Job ends ****"
date
