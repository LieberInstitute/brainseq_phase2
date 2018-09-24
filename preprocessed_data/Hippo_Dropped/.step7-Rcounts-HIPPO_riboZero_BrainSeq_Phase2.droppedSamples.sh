#!/bin/bash
#$ -cwd
#$ -pe local 8
#$ -l mem_free=5G,h_vmem=6G,h_fsize=200G
#$ -N step7-Rcounts-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples
#$ -o ./logs/Rcounts-HIPPO_riboZero_BrainSeq_Phase2.txt
#$ -e ./logs/Rcounts-HIPPO_riboZero_BrainSeq_Phase2.txt
#$ -hold_jid pipeline_setup,step4-featCounts-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples,step6-txQuant-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"


Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/.create_count_objects-human.R -o hg38 -m /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped -e HIPPO_riboZero_BrainSeq_Phase2 -p droppedSamples -l TRUE -c TRUE -t 8 -s reverse

echo "**** Job ends ****"
date
