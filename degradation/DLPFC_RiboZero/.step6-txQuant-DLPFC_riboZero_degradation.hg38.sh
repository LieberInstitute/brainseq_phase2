#!/bin/bash
#$ -cwd
#$ -l mem_free=80G,h_vmem=90G,h_fsize=100G
#$ -N step6-txQuant-DLPFC_riboZero_degradation.hg38
#$ -pe local 1
#$ -o ./logs/txQuant-DLPFC_riboZero_degradation.$TASK_ID.txt
#$ -e ./logs/txQuant-DLPFC_riboZero_degradation.$TASK_ID.txt
#$ -t 1-20
#$ -tc 15
#$ -hold_jid pipeline_setup,step4-featCounts-DLPFC_riboZero_degradation.hg38
#$ -m a
echo "**** Job starts ****"
date

FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")


mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Salmon_tx/${ID}

if [ TRUE == "TRUE" ] ; then 
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Salmon-0.7.2_linux_x86_64/bin/salmon quant 	-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/transcripts/salmon_index_gencode.v25.transcripts -p 8 -l ISR 	-1 ${FILE1} -2 ${FILE2} 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Salmon_tx/${ID}
else
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Salmon-0.7.2_linux_x86_64/bin/salmon quant 	-i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/transcripts/salmon_index_gencode.v25.transcripts -p 8 -l U 	-r ${FILE1} 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/DLPFC_RiboZero/Salmon_tx/${ID}
fi


echo "**** Job ends ****"
date
