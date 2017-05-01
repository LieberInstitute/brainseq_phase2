#!/bin/bash
#$ -cwd
#$ -pe local 8
#$ -l mem_free=5G,h_vmem=6G,h_fsize=200G
#$ -N step7-Rcounts-DLPFC_riboZero_degradation.hg38
#$ -o ./logs/Rcounts-DLPFC_riboZero_degradation.txt
#$ -e ./logs/Rcounts-DLPFC_riboZero_degradation.txt
#$ -hold_jid pipeline_setup,step4-featCounts-DLPFC_riboZero_degradation.hg38,step6-txQuant-DLPFC_riboZero_degradation.hg38
#$ -m a
echo "**** Job starts ****"
date


Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/.create_count_objects-human.R -o hg38 -m /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero -e DLPFC_riboZero_degradation -p hg38 -l TRUE -c FALSE -t 8

echo "**** Job ends ****"
date
