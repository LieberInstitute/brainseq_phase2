#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=25G,h_vmem=38G
#$ -N salmon-brainseq.dlpfc
#$ -o ./logs/salmon-brainseq.o.$TASK_ID.txt
#$ -e ./logs/salmon-brainseq.e.$TASK_ID.txt
#$ -t 1-579
#$ -tc 20

echo "**** Job starts ****"
date

FILEID=$(awk "NR==${SGE_TASK_ID}" /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/SAMPLE_IDs.txt )
ID=$(basename "${FILEID}")

mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Salmon_tx/${ID}

/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Salmon-0.7.2_linux_x86_64/bin/salmon \
quant -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/transcripts/salmon_index_gencode.v25.transcripts \
-l ISR -1 ${FILEID}_R1_001.fastq.gz -2 ${FILEID}_R2_001.fastq.gz \
-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Salmon_tx/${ID}


echo "**** Job ends ****"
date
