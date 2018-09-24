#!/bin/bash
#$ -cwd
#$ -pe local 8
#$ -l bluejay,mem_free=20G,h_vmem=24G,h_fsize=200G
#$ -N Rcounts-brainseq.dlpfc
#$ -o ./logs/Rcounts-brainseq.txt
#$ -e ./logs/Rcounts-brainseq.txt
#$ -hold_jid pipeline_setup,featCounts-brainseq.dlpfc
#$ -m a
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\na2b1cda9a67211320865c002ee1a20429cce32b1"

Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/.create_count_objects-human.R

echo "**** Job ends ****"
date
