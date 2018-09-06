#!/bin/bash
#$ -cwd
#$ -l mem_free=250G,h_vmem=250G,h_fsize=200G
#$ -N hippo_subset
#$ -o logs/hippo_subset.txt
#$ -e logs/hippo_subset.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full_GTEx/results_subset_hippo.R

echo "**** Job ends ****"
date
