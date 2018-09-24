#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=10G
#$ -N step8b-mergeVariantCalls-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples
#$ -o ./logs/mergeVariantCalls-HIPPO_riboZero_BrainSeq_Phase2.txt
#$ -e ./logs/mergeVariantCalls-HIPPO_riboZero_BrainSeq_Phase2.txt
#$ -hold_jid pipeline_setup,step8-callVariants-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

VCFS=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/samples.manifest | awk '{print "/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/Genotypes/"$NF".vcf.gz"}' | paste -sd " ")

module load vcftools
module load htslib
vcf-merge ${VCFS} | bgzip -c > /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/Genotypes/mergedVariants.vcf.gz

echo "**** Job ends ****"
date
