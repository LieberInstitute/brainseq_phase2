#!/bin/bash
#$ -cwd
#$ -l mem_free=60G,h_vmem=80G,h_fsize=100G
#$ -N coverage-brainseq.hippo
#$ -o ./logs/coverage-brainseq.o.$TASK_ID.txt
#$ -e ./logs/coverage-brainseq.e.$TASK_ID.txt
#$ -t 1-452
#$ -tc 40
#$ -hold_jid pipeline_setup,featCounts-brainseq.hippo
#$ -m a
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\n3de0c0d9f6fd9a855ca1cad7cba3b746f3ef6f4f"

FILEID=$(awk "NR==${SGE_TASK_ID}" /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/SAMPLE_IDs.txt )
ID=$(basename "${FILEID}")

## Normalizing bigwigs to 40 million 100 bp reads
module load python/2.7.9
module load ucsctools
python ~/.local/bin/bam2wig.py -s /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -i /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/HISAT2_out/${ID}_accepted_hits.sorted.bam -t 4000000000 -o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/${ID}

## Remove temp files
rm /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Coverage/${ID}.wig

echo "**** Job ends ****"
date
