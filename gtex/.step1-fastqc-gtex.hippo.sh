#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=7G,h_fsize=100G
#$ -N step1-fastqc-gtex.hippo
#$ -o ./logs/fastqc-gtex.$TASK_ID.txt
#$ -e ./logs/fastqc-gtex.$TASK_ID.txt
#$ -t 1-94
#$ -tc 100
#$ -hold_jid pipeline_setup,step00-merge-gtex.hippo
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
echo "Sample id: $(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

## Locate file and ids
FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")

mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/FastQC/Untrimmed/${ID}

if [ TRUE == "TRUE" ]
then 
    /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/FastQC_v0.11.5/fastqc ${FILE1} ${FILE2} --outdir=/dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/FastQC/Untrimmed/${ID} --extract
else
    /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/FastQC_v0.11.5/fastqc ${FILE1} --outdir=/dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/FastQC/Untrimmed/${ID} --extract
fi

echo "**** Job ends ****"
date
