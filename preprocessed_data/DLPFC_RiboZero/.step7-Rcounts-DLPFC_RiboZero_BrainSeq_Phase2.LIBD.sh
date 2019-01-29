#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l bluejay,mem_free=28G,h_vmem=30G,h_fsize=200G
#$ -N step7-Rcounts-DLPFC_RiboZero_BrainSeq_Phase2.LIBD
#$ -o ./logs/Rcounts-DLPFC_RiboZero_BrainSeq_Phase2.txt
#$ -e ./logs/Rcounts-DLPFC_RiboZero_BrainSeq_Phase2.txt
#$ -hold_jid pipeline_setup,step4-featCounts-DLPFC_RiboZero_BrainSeq_Phase2.LIBD,step6-txQuant-DLPFC_RiboZero_BrainSeq_Phase2.LIBD
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/.create_count_objects-human.R -o hg38 -m /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero -e DLPFC_RiboZero_BrainSeq_Phase2 -p LIBD -l TRUE -c TRUE -t 5 -s reverse

echo "**** Job ends ****"
date