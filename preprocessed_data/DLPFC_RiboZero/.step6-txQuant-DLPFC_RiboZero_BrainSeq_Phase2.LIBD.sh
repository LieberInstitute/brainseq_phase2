#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=60G,h_fsize=100G
#$ -N step6-txQuant-DLPFC_RiboZero_BrainSeq_Phase2.LIBD
#$ -pe local 2
#$ -o ./logs/txQuant-DLPFC_RiboZero_BrainSeq_Phase2.$TASK_ID.txt
#$ -e ./logs/txQuant-DLPFC_RiboZero_BrainSeq_Phase2.$TASK_ID.txt
#$ -t 1-0
#$ -tc 15
#$ -hold_jid pipeline_setup,step4-featCounts-DLPFC_RiboZero_BrainSeq_Phase2.LIBD
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

FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ FALSE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")


mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Salmon_tx/${ID}

if [ FALSE == "TRUE" ] ; then 
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Salmon-0.7.2_linux_x86_64/bin/salmon quant 	-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/transcripts/salmon_index_gencode.v25.transcripts -p 1 -l SR 	-1 ${FILE1} -2 ${FILE2} 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Salmon_tx/${ID}
else
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Salmon-0.7.2_linux_x86_64/bin/salmon quant 	-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/transcripts/salmon_index_gencode.v25.transcripts -p 1 -l SR 	-r ${FILE1} 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Salmon_tx/${ID}
fi


echo "**** Job ends ****"
date
