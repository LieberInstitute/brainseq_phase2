#!/bin/bash
#$ -cwd
#$ -l mem_free=6G,h_vmem=10G,h_fsize=100G
#$ -N ercc-brainseq.hippo
#$ -pe local 6
#$ -o ./logs/ercc-brainseq.o.$TASK_ID.txt
#$ -e ./logs/ercc-brainseq.e.$TASK_ID.txt
#$ -t 1-20
#$ -tc 20
#$ -hold_jid pipeline_setup,merge-brainseq.hippo
#$ -m a
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\n3de0c0d9f6fd9a855ca1cad7cba3b746f3ef6f4f"

FILEID=$(awk "NR==${SGE_TASK_ID}" /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/SAMPLE_IDs.txt )
ID=$(basename "${FILEID}")
mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Ercc/${ID}

if [ TRUE == "TRUE" ] ; then 
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/kallisto quant -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx -o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Ercc/${ID} -t 6 --rf-stranded ${FILEID}_R1_001.fastq.gz ${FILEID}_R2_001.fastq.gz
else
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/kallisto quant -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx -o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Ercc/${ID} -t 6 --single ${FILEID}.fastq.gz

fi

echo "**** Job ends ****"
date
