#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=40G,h_vmem=40G,h_fsize=100G
#$ -N limma_reg_specific_jxn_adult
#$ -o ./logs/limma_reg_specific_jxn_adult.txt
#$ -e ./logs/limma_reg_specific_jxn_adult.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: lcollado"
echo "Job id: "
echo "Job name: "
echo "Hostname: jhpce01"
echo "Task id: "

module load conda_R/3.4.x
Rscript limma_reg_specific.R -t jxn -a adult

echo "**** Job ends ****"
date
