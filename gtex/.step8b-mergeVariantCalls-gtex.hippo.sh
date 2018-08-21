#!/bin/bash
#$ -cwd
#$ -l mem_free=8G,h_vmem=10G
#$ -N step8b-mergeVariantCalls-gtex.hippo
#$ -o ./logs/mergeVariantCalls-gtex.txt
#$ -e ./logs/mergeVariantCalls-gtex.txt
#$ -hold_jid pipeline_setup,step8-callVariants-gtex.hippo
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

VCFS=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/samples.manifest | awk '{print "/dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/Genotypes/"$NF".vcf.gz"}' | paste -sd " ")

module load vcftools
module load htslib
vcf-merge ${VCFS} | bgzip -c > /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/Genotypes/mergedVariants.vcf.gz

echo "**** Job ends ****"
date
