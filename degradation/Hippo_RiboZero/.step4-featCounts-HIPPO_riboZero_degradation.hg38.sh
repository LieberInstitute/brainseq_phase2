#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=12G,h_fsize=100G
#$ -N step4-featCounts-HIPPO_riboZero_degradation.hg38
#$ -pe local 8
#$ -o ./logs/featCounts-HIPPO_riboZero_degradation.$TASK_ID.txt
#$ -e ./logs/featCounts-HIPPO_riboZero_degradation.$TASK_ID.txt
#$ -t 1-12
#$ -tc 10
#$ -hold_jid pipeline_setup,step3-hisat2-HIPPO_riboZero_degradation.hg38
#$ -m a
echo "**** Job starts ****"
date


# Directories
mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Counts/gene
mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Counts/exon
mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Counts/junction/tmpdir

FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
BAM=/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/HISAT2_out/${ID}_accepted_hits.sorted.bam

if [ TRUE == "TRUE" ] ; then 
	# genes	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts 	-s 2 -p -T 8 -a /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Counts/gene/${ID}_Gencode.v25.hg38_Genes.counts $BAM
	# exons	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts 	-s 2 -p -O -f -T 8 -a /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Counts/exon/${ID}_Gencode.v25.hg38_Exons.counts $BAM
else
	# genes	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts 	-T 8 -a /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Counts/gene/${ID}_Gencode.v25.hg38_Genes.counts $BAM
	# exons	
	/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts 	-O -f -T 8 -a /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf 	-o /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Counts/exon/${ID}_Gencode.v25.hg38_Exons.counts $BAM
fi
	
# junctions	
OUTJXN=/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Counts/junction/${ID}_junctions_primaryOnly_regtools.bed
OUTCOUNT=/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Counts/junction/${ID}_junctions_primaryOnly_regtools.count
TMPDIR=/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Counts/junction/tmpdir
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
