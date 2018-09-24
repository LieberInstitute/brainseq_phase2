#!/bin/bash
#$ -cwd
#$ -pe local 8
#$ -l bluejay,mem_free=20G,h_vmem=24G,h_fsize=200G
#$ -N Rcounts-brainseq.hippo
#$ -o ./logs/Rcounts-brainseq.txt
#$ -e ./logs/Rcounts-brainseq.txt
#$ -hold_jid pipeline_setup,featCounts-brainseq.hippo
#$ -m a
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\n3de0c0d9f6fd9a855ca1cad7cba3b746f3ef6f4f"

Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/.create_count_objects-human.R

echo "**** Job ends ****"
date
