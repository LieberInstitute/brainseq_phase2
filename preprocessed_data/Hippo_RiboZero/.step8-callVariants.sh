#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=8G,h_fsize=100G
#$ -N varCall-brainseq.hippo
#$ -pe local 1
#$ -o ./logs/varCall-brainseq.o.$TASK_ID.txt
#$ -e ./logs/varCall-brainseq.e.$TASK_ID.txt
#$ -t 1-483
#$ -tc 50
#$ -m a
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\n3de0c0d9f6fd9a855ca1cad7cba3b746f3ef6f4f"

mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Genotypes/

FILEID=$(awk "NR==${SGE_TASK_ID}" /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/SAMPLE_IDs.txt )
ID=$(basename "${FILEID}")
BAM=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/HISAT2_out/${ID}_accepted_hits.sorted.bam

SNPBED=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed
SNPTMP=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Genotypes/${ID}_calledVariants.tmp
SNPOUT=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Genotypes/${ID}_calledVariants.txt
GENOME=/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/GRCh38.primary_assembly.genome.fa
module load samtools
samtools mpileup -l $SNPBED -AB -q0 -Q0 -d1000000 -t DP -f $GENOME $BAM | /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/pileVar.pl > $SNPOUT
rm $SNPTMP