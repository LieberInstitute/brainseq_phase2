#!/bin/bash
#$ -cwd
#$ -l mem_free=160G,h_vmem=160G,h_fsize=200G
#$ -N bsII_dlpfc_eQTL
#$ -o ./eqtl_tables/logs/eQTL_dlpfc_4_features.txt
#$ -e ./eqtl_tables/logs/eQTL_dlpfc_4_features.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/run_eqtls_dlpfc.R

echo "**** Job ends ****"
date
