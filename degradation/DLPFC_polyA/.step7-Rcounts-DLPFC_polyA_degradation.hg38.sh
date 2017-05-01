#!/bin/bash
#$ -cwd
#$ -pe local 8
#$ -l mem_free=5G,h_vmem=6G,h_fsize=200G
#$ -N step7-Rcounts-DLPFC_polyA_degradation.hg38
#$ -o ./logs/Rcounts-DLPFC_polyA_degradation.txt
#$ -e ./logs/Rcounts-DLPFC_polyA_degradation.txt
#$ -hold_jid pipeline_setup,step4-featCounts-DLPFC_polyA_degradation.hg38,step6-txQuant-DLPFC_polyA_degradation.hg38
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"


Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_polyA/.create_count_objects-human.R -o hg38 -m /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_polyA -e DLPFC_polyA_degradation -p hg38 -l TRUE -c FALSE -t 8 -s FALSE

echo "**** Job ends ****"
date
