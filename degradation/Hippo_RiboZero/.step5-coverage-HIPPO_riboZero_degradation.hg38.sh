#!/bin/bash
#$ -cwd
#$ -l mem_free=30G,h_vmem=40G,h_fsize=100G
#$ -N step5-coverage-HIPPO_riboZero_degradation.hg38
#$ -o ./logs/coverage-HIPPO_riboZero_degradation.$TASK_ID.txt
#$ -e ./logs/coverage-HIPPO_riboZero_degradation.$TASK_ID.txt
#$ -t 1-12
#$ -tc 40
#$ -hold_jid pipeline_setup,step3-hisat2-HIPPO_riboZero_degradation.hg38,step3b-infer-strandness-HIPPO_riboZero_degradation.hg38
#$ -m a
echo "**** Job starts ****"
date

if [ ! -f "inferred_strandness_pattern.txt" ]
then
    echo "Missing the file inferred_strandness_pattern.txt"
    exit 1
fi

FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [[ TRUE == "TRUE" ]]
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
STRANDRULE=$(cat inferred_strandness_pattern.txt)

## Normalizing bigwigs to 40 million 100 bp reads
module load python/2.7.9
module load ucsctools
python ~/.local/bin/bam2wig.py -s /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -i /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/HISAT2_out/${ID}_accepted_hits.sorted.bam -t 4000000000 -o /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Coverage/${ID} -d "${STRANDRULE}"

## Remove temp files
rm /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Coverage/${ID}*.wig

echo "**** Job ends ****"
date
