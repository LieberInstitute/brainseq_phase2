#!/bin/bash
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -N limma_reg_specific_exon_adult
#$ -o ./logs/limma_reg_specific_exon_adult.txt
#$ -e ./logs/limma_reg_specific_exon_adult.txt
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
Rscript limma_reg_specific.R -t exon -a adult

echo "**** Job ends ****"
date
