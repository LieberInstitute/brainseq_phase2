#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=14G,h_fsize=200G
#$ -N fix_step3
#$ -o ./logs/fix_step3.$TASK_ID.txt
#$ -e ./logs/fix_step3.$TASK_ID.txt
#$ -t 1-579
#$ -tc 15
#$ -m a

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

## Locate file and ids
FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ FALSE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")

SORTEDBAM=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/${ID}_accepted_hits.sorted

## Run infer experiment
echo "**** Inferring strandedness with infer_experiment.py ****"
date
module load python/2.7.9
~/.local/bin/infer_experiment.py -i ${SORTEDBAM}.bam -r /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/RSeQC/hg38.bed 1> /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/infer_strandness/${ID}.txt 2>&1

echo "**** Job ends ****"
date
