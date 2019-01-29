#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=14G,h_fsize=100G
#$ -N fastqc-brainseq.dlpfc
#$ -o ./logs/fastqc-brainseq.o.$TASK_ID.txt
#$ -e ./logs/fastqc-brainseq.e.$TASK_ID.txt
#$ -t 1-579
#$ -tc 100
#$ -hold_jid pipeline_setup,merge-brainseq.dlpfc
#$ -m a
echo "**** Job starts ****"
date


echo -e "**** Pipeline version: GitHub sha ****\na2b1cda9a67211320865c002ee1a20429cce32b1"

FILEID=$(awk "NR==${SGE_TASK_ID}" /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/SAMPLE_IDs.txt )
ID=$(basename "${FILEID}")
mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/FastQC/Untrimmed/${ID}

if [ TRUE == "TRUE" ]
then 
    /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/FastQC_v0.11.5/fastqc ${FILEID}_R1_001.fastq.gz ${FILEID}_R2_001.fastq.gz --outdir=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/FastQC/Untrimmed/${ID} --extract
else
    /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/FastQC_v0.11.5/fastqc ${FILEID}.fastq.gz --outdir=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/FastQC/Untrimmed/${ID} --extract
fi

echo "**** Job ends ****"
date