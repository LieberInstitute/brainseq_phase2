#!/bin/bash
#$ -cwd
#$ -N step3b-infer-strandness-HIPPO_riboZero_degradation.hg38
#$ -o ./logs/infer-strandness-HIPPO_riboZero_degradation.txt
#$ -e ./logs/infer-strandness-HIPPO_riboZero_degradation.txt
#$ -hold_jid pipeline_setup,step3-hisat2-HIPPO_riboZero_degradation.hg38
#$ -m a
echo "**** Job starts ****"
date

## Process the infer experiment info
Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step3b_infer_strandness.R

echo "**** Job ends ****"
date
