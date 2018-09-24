#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=18G,h_vmem=20G,h_fsize=100G
#$ -N step5-coverage-DLPFC_RiboZero_BrainSeq_Phase2.LIBD
#$ -o ./logs/coverage-DLPFC_RiboZero_BrainSeq_Phase2.$TASK_ID.txt
#$ -e ./logs/coverage-DLPFC_RiboZero_BrainSeq_Phase2.$TASK_ID.txt
#$ -t 1-579
#$ -tc 40
#$ -hold_jid pipeline_setup,step3-hisat2-DLPFC_RiboZero_BrainSeq_Phase2.LIBD,step3b-infer-strandness-DLPFC_RiboZero_BrainSeq_Phase2.LIBD
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
echo "Sample id: $(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

if [ ! -f "inferred_strandness_pattern.txt" ]
then
    echo "Missing the file inferred_strandness_pattern.txt"
    exit 1
fi

FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [[ TRUE == "TRUE" ]]
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
STRANDRULE=$(cat inferred_strandness_pattern.txt)

## Can only use -d when the data is stranded
if [ ${STRANDRULE} == "none" ]
then
    STRANDOPTION=""
else
    STRANDOPTION="-d \"${STRANDRULE}\""
fi

## Normalizing bigwigs to 40 million 100 bp reads
module load python/2.7.9
module load ucsctools
#python ~/.local/bin/bam2wig.py -s /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -i /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/${ID}_accepted_hits.sorted.bam -t 4000000000 -o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Coverage/${ID} "${STRANDOPTION}"
python ~/.local/bin/bam2wig.py -s /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -i /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/${ID}_accepted_hits.sorted.bam -t 4000000000 -o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Coverage/${ID} -d "1+-,1-+,2++,2--"


## Remove temp files
rm /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Coverage/${ID}*.wig

echo "**** Job ends ****"
date
