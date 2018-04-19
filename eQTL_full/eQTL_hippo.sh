#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=200G
#$ -N bsII_hippo_eQTL
#$ -o ./eqtl_tables/logs/eQTL_hippo_4_features.txt
#$ -e ./eqtl_tables/logs/eQTL_hippo_4_features.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/run_eqtls_hippo.R

echo "**** Job ends ****"
date
