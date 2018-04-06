#!/bin/bash
#$ -cwd
#$ -pe local 5
#$ -l mem_free=18G,h_vmem=20G,h_fsize=200G
#$ -N step7-Rcounts-gtex.hippo
#$ -o ./logs/Rcounts-gtex.txt
#$ -e ./logs/Rcounts-gtex.txt
#$ -hold_jid pipeline_setup,step4-featCounts-gtex.hippo,step6-txQuant-gtex.hippo
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

## Force R 3.3.x in JHPCE (to avoid some issues with conda_R)
module unload conda_R
module load R/3.3.x

Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/.step7-create_count_objects-human.R -o hg38 -m /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex -e gtex -p hippo -l TRUE -c FALSE -t 5 -s reverse

echo "**** Job ends ****"
date
