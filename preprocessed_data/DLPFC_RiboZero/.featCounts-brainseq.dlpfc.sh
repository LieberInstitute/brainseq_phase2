#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=20G,h_vmem=24G,h_fsize=100G
#$ -N featCounts-brainseq.dlpfc
#$ -pe local 8
#$ -o ./logs/featCounts-brainseq.o.$TASK_ID.txt
#$ -e ./logs/featCounts-brainseq.e.$TASK_ID.txt
#$ -t 1-579
#$ -tc 10
#$ -hold_jid pipeline_setup,hisat2-brainseq.dlpfc
#$ -m a
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\na2b1cda9a67211320865c002ee1a20429cce32b1"

FILEID=$(awk "NR==${SGE_TASK_ID}" /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/SAMPLE_IDs.txt )
ID=$(basename "${FILEID}")
BAM=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/HISAT2_out/${ID}_accepted_hits.sorted.bam

if [ TRUE == "TRUE" ] ; then 
	# genes	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts 	-s 2 -p -T 8 -a /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Counts/gene/${ID}_Gencode.v25.hg38_Genes.counts $BAM
	# exons	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts 	-s 2 -p -O -f -T 8 -a /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Counts/exon/${ID}_Gencode.v25.hg38_Exons.counts $BAM
else
	# genes	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts 	-T 8 -a /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Counts/gene/${ID}_Gencode.v25.hg38_Genes.counts $BAM
	# exons	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts 	-O -f -T 8 -a /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Counts/exon/${ID}_Gencode.v25.hg38_Exons.counts $BAM
fi
	
# junctions	
OUTJXN=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Counts/junction/${ID}_junctions_primaryOnly_regtools.bed
OUTCOUNT=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Counts/junction/${ID}_junctions_primaryOnly_regtools.count
TMPDIR=/dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Counts/junction/tmpdir
TMPBAM=${TMPDIR}/${ID}.bam
#filter only primary alignments
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools view -@ 8 -bh -F 0x100 $BAM > ${TMPBAM}
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools index ${TMPBAM}

## Load python 2.7.9 since the default one cannot run:
# python
# import site
module load python/2.7.9
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/regtools/build/regtools junctions extract -i 9 -o ${OUTJXN} ${TMPBAM}
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/bed_to_juncs_withCount < ${OUTJXN} > ${OUTCOUNT}


echo "**** Job ends ****"
date
