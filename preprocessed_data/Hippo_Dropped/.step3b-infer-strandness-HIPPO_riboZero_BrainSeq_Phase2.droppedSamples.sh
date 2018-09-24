#!/bin/bash
#$ -cwd
#$ -N step3b-infer-strandness-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples
#$ -o ./logs/infer-strandness-HIPPO_riboZero_BrainSeq_Phase2.txt
#$ -e ./logs/infer-strandness-HIPPO_riboZero_BrainSeq_Phase2.txt
#$ -hold_jid pipeline_setup,step3-hisat2-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Process the infer experiment info
Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step3b_infer_strandness.R -o "HISAT2_out/infer_strandess" -p "inferred_strandness_pattern.txt"

echo "**** Job ends ****"
date
