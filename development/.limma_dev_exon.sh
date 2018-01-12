#!/bin/bash
#$ -cwd
#$ -l mem_free=45G,h_vmem=45G,h_fsize=100G
#$ -N limma_dev_exon
#$ -o ./logs/limma_dev_exon.txt
#$ -e ./logs/limma_dev_exon.txt
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
Rscript limma_dev.R -t exon

echo "**** Job ends ****"
date
