#!/bin/bash
#$ -cwd
#$ -l mem_free=40G,h_vmem=45G,h_fsize=100G
#$ -N step6-txQuant-gtex.dlpfc
#$ -pe local 1
#$ -o ./logs/txQuant-gtex.$TASK_ID.txt
#$ -e ./logs/txQuant-gtex.$TASK_ID.txt
#$ -t 1-108
#$ -tc 15
#$ -hold_jid pipeline_setup,step4-featCounts-gtex.dlpfc
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
echo "Sample id: $(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_dlpfc/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_dlpfc/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_dlpfc/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_dlpfc/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")


mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_dlpfc/Salmon_tx/${ID}

if [ TRUE == "TRUE" ] ; then 
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Salmon-0.8.2_linux_x86_64/bin/salmon quant 	-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/transcripts/salmon_0.8.2_index_gencode.v25.transcripts -p 1 -l IU 	-1 ${FILE1} -2 ${FILE2} 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_dlpfc/Salmon_tx/${ID}
else
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Salmon-0.8.2_linux_x86_64/bin/salmon quant 	-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/transcripts/salmon_0.8.2_index_gencode.v25.transcripts -p 1 -l IU 	-r ${FILE1} 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex_dlpfc/Salmon_tx/${ID}
fi


echo "**** Job ends ****"
date
