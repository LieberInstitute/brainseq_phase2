#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=20G,h_vmem=30G,h_fsize=100G
#$ -N trim-brainseq.dlpfc
#$ -pe local 8
#$ -o ./logs/trim-brainseq.o.$TASK_ID.txt
#$ -e ./logs/trim-brainseq.e.$TASK_ID.txt
#$ -t 1-579
#$ -tc 8
#$ -hold_jid pipeline_setup,fastqc-brainseq.dlpfc
#$ -m a
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\na2b1cda9a67211320865c002ee1a20429cce32b1"

FILEID=$(awk "NR==${SGE_TASK_ID}" /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/SAMPLE_IDs.txt )
ID=$(basename "${FILEID}")

if [ TRUE == "TRUE" ] ; then 
	REPORT1=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/FastQC/Untrimmed/${ID}/${ID}_R1_001_fastqc/summary.txt
	REPORT2=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/FastQC/Untrimmed/${ID}/${ID}_R2_001_fastqc/summary.txt
	RESULT1=$(grep "Adapter Content" $REPORT1 | cut -c1-4)
	RESULT2=$(grep "Adapter Content" $REPORT2 | cut -c1-4)

	if [[ $RESULT1 == "FAIL" || $RESULT2 == "FAIL" ]] ; then
		## trim, rerun fastQC
		echo "End 1 adapters: $RESULT1"
		echo "End 2 adapters: $RESULT2"
		echo "Trimming will occur."
		
		mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq
		FP=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq/${ID}_trimmed_forward_paired.fq.gz
		FU=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq/${ID}_trimmed_forward_unpaired.fq.gz
		RP=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq/${ID}_trimmed_reverse_paired.fq.gz
		RU=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq/${ID}_trimmed_reverse_unpaired.fq.gz
		
		## trim adapters
		java -jar /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 		${FILEID}_R1_001.fastq.gz ${FILEID}_R2_001.fastq.gz $FP $FU $RP $RU 		ILLUMINACLIP:/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10:1 		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75

		## rerun fastqc
		mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/FastQC/Trimmed/${ID}
		/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/FastQC_v0.11.5/fastqc 		$FP $FU $RP $RU 		--outdir=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/FastQC/Trimmed/${ID} --extract
	else
		echo "No trimming required!"
	fi

else
	## reads are single-end
	REPORT1=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/FastQC/Untrimmed/${ID}/${ID}_fastqc/summary.txt
	RESULT1=$(grep "Adapter Content" $REPORT1 | cut -c1-4)

	if [ $RESULT1 == "FAIL" ] ; then
		## trim, rerun fastQC
		echo "Adapters: $RESULT1"
		echo "Trimming will occur."
		
		mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq
		OUT=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq/${ID}_trimmed.fastq.gz
		
		## trim adapters
		java -jar /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 		${FILEID}.fastq.gz $OUT 		ILLUMINACLIP:/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/Trimmomatic-0.36/adapters/TruSeq2-SE.fa:2:30:10:1 		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		## rerun fastqc
		mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/FastQC/Trimmed/${ID}
		/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/FastQC_v0.11.5/fastqc $OUT 		--outdir=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/FastQC/Trimmed/${ID} --extract
	else
		echo "No trimming required!"
	fi
fi

	
echo "**** Job ends ****"
date
