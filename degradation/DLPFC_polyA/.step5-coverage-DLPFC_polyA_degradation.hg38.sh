#!/bin/bash
#$ -cwd
#$ -l mem_free=30G,h_vmem=40G,h_fsize=100G
#$ -N step5-coverage-DLPFC_polyA_degradation.hg38
#$ -o ./logs/coverage-DLPFC_polyA_degradation.$TASK_ID.txt
#$ -e ./logs/coverage-DLPFC_polyA_degradation.$TASK_ID.txt
#$ -t 1-20
#$ -tc 40
#$ -hold_jid pipeline_setup,step3-hisat2-DLPFC_polyA_degradation.hg38,step3b-infer-strandness-DLPFC_polyA_degradation.hg38
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

if [ ! -f "inferred_strandness_pattern.txt" ]
then
    echo "Missing the file inferred_strandness_pattern.txt"
    exit 1
fi

FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_polyA/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [[ TRUE == "TRUE" ]]
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_polyA/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_polyA/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
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
python ~/.local/bin/bam2wig.py -s /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -i /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_polyA/HISAT2_out/${ID}_accepted_hits.sorted.bam -t 4000000000 -o /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_polyA/Coverage/${ID} "${STRANDOPTION}"

## Remove temp files
rm /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_polyA/Coverage/${ID}*.wig

echo "**** Job ends ****"
date
