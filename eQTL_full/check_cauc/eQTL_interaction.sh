#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=200G,h_vmem=200G,h_fsize=200G
#$ -N inter_eQTL_bsII
#$ -o eqtl_tables/logs/eQTL_interaction_4_features2.txt
#$ -e eqtl_tables/logs/eQTL_interaction_4_features2.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/run_eqtls_interaction.R

echo "**** Job ends ****"
date
