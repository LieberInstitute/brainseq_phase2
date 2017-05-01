#!/bin/bash
#$ -cwd
#$ -N step4-featCounts-DLPFC_polyA_degradation.hg38_clean
#$ -o ./logs/featCounts-DLPFC_polyA_degradation_clean.txt
#$ -e ./logs/featCounts-DLPFC_polyA_degradation_clean.txt
#$ -hold_jid pipeline_setup,step4-featCounts-DLPFC_polyA_degradation.hg38
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
rm -rf /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_polyA/Counts/junction/tmpdir

echo "**** Job ends ****"
date
