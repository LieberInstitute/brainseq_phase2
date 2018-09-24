#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=8G,h_fsize=100G
#$ -N step8-callVariants-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples
#$ -o ./logs/callVariants-HIPPO_riboZero_BrainSeq_Phase2.$TASK_ID.txt
#$ -e ./logs/callVariants-HIPPO_riboZero_BrainSeq_Phase2.$TASK_ID.txt
#$ -t 1-54
#$ -tc 50
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
echo "****"
echo "Sample id: $(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/Genotypes/

ID=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
BAM=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/HISAT2_out/${ID}_accepted_hits.sorted.bam

SNPTMP=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/Genotypes/${ID}_tmp.vcf
SNPOUTGZ=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/Genotypes/${ID}.vcf.gz

module load bcftools
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools mpileup -l /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed -AB -q0 -Q13 -d1000000 -uf /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/GRCh38.primary_assembly.genome.fa ${BAM} -o ${SNPTMP}
bcftools call -mv -Oz ${SNPTMP} > ${SNPOUTGZ}

module load htslib
tabix -p vcf ${SNPOUTGZ}

rm ${SNPTMP}

echo "**** Job ends ****"
date
