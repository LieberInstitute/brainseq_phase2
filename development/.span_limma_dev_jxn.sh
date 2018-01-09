#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=45G,h_vmem=45G,h_fsize=100G
#$ -N span_limma_dev_jxn
#$ -o ./logs/span_limma_dev_jxn.txt
#$ -e ./logs/span_limma_dev_jxn.txt
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
Rscript span_limma_dev.R -t jxn

echo "**** Job ends ****"
date
