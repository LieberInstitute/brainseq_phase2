#!/bin/bash
#$ -cwd
#$ -N step4-featCounts-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples_clean
#$ -o ./logs/featCounts-HIPPO_riboZero_BrainSeq_Phase2_clean.txt
#$ -e ./logs/featCounts-HIPPO_riboZero_BrainSeq_Phase2_clean.txt
#$ -hold_jid pipeline_setup,step4-featCounts-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"


## Delete temporary files after they have been used
rm -rf /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/Counts/junction/tmpdir

echo "**** Job ends ****"
date
