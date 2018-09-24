#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=14G,h_fsize=100G
#$ -N hisat2-brainseq.dlpfc
#$ -pe local 8
#$ -o ./logs/hisat2-brainseq.o.$TASK_ID.txt
#$ -e ./logs/hisat2-brainseq.e.$TASK_ID.txt
#$ -t 1-579
#$ -tc 15
#$ -hold_jid pipeline_setup,trim-brainseq.dlpfc
#$ -m a
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\na2b1cda9a67211320865c002ee1a20429cce32b1"

FILEID=$(awk "NR==${SGE_TASK_ID}" /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/SAMPLE_IDs.txt )
ID=$(basename "${FILEID}")

if [ -e /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq/${ID}_trimmed_forward_paired.fq.gz ] ; then
	## Trimmed, paired-end
	echo "HISAT2 alignment run on trimmed paired-end reads"
	FP=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq/${ID}_trimmed_forward_paired.fq.gz
	FU=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq/${ID}_trimmed_forward_unpaired.fq.gz
	RP=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq/${ID}_trimmed_reverse_paired.fq.gz
	RU=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq/${ID}_trimmed_reverse_unpaired.fq.gz
	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2 -p 8 	-x /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/hisat2_GRCh38primary -1 $FP -2 $RP -U ${FU},${RU} 	-S /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/${ID}_hisat_out.sam --rna-strandness RF --phred33 	2>/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/align_summaries/${ID}_summary.txt
	
elif  [ -e /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq/${ID}_trimmed.fq.gz ] ; then
	## Trimmed, single-end
	echo "HISAT2 alignment run on trimmed single-end reads"
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2 -p 8 	-x /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/hisat2_GRCh38primary -U /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/trimmed_fq/${ID}_trimmed.fq.gz 	-S /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/${ID}_hisat_out.sam --phred33 	2>/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/align_summaries/${ID}_summary.txt

elif [ TRUE == "TRUE" ] ; then
	## Untrimmed, pair-end
	echo "HISAT2 alignment run on original untrimmed paired-end reads"
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2 -p 8 	-x /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/hisat2_GRCh38primary -1 ${FILEID}_R1_001.fastq.gz -2 ${FILEID}_R2_001.fastq.gz 	-S /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/${ID}_hisat_out.sam --rna-strandness RF --phred33 	2>/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/align_summaries/${ID}_summary.txt

else
	## Untrimmed, single-end
	echo "HISAT2 alignment run on original untrimmed single-end reads"
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2 -p 8 	-x /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/hisat2_GRCh38primary -U ${FILEID}.fastq.gz 	-S /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/${ID}_hisat_out.sam --phred33 	2>/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/align_summaries/${ID}_summary.txt
fi


###sam to bam
SAM=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/${ID}_hisat_out.sam
BAMACC=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/${ID}_accepted_hits.bam
BAMS=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/${ID}_accepted_hits.sorted

#filter unmapped segments
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools view -bh -F 4 $SAM > $BAMACC
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools sort -@ 8 $BAMACC $BAMS
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools index ${BAMS}.bam

rm $SAM

echo "**** Job ends ****"
date
