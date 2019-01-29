#!/bin/bash
#$ -cwd
#$ -l mem_free=80G,h_vmem=90G,h_fsize=100G
#$ -N step6-txQuant-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples
#$ -pe local 1
#$ -o ./logs/txQuant-HIPPO_riboZero_BrainSeq_Phase2.$TASK_ID.txt
#$ -e ./logs/txQuant-HIPPO_riboZero_BrainSeq_Phase2.$TASK_ID.txt
#$ -t 1-54
#$ -tc 15
#$ -hold_jid pipeline_setup,step4-featCounts-HIPPO_riboZero_BrainSeq_Phase2.droppedSamples
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")


mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/Salmon_tx/${ID}

if [ TRUE == "TRUE" ] ; then 
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Salmon-0.7.2_linux_x86_64/bin/salmon quant 	-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/transcripts/salmon_index_gencode.v25.transcripts -p 1 -l ISR 	-1 ${FILE1} -2 ${FILE2} 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/Salmon_tx/${ID}
else
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Salmon-0.7.2_linux_x86_64/bin/salmon quant 	-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/transcripts/salmon_index_gencode.v25.transcripts -p 1 -l ISR 	-r ${FILE1} 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_Dropped/Salmon_tx/${ID}
fi


echo "**** Job ends ****"
date