#!/bin/bash
#$ -cwd
#$ -pe local 10
#$ -l mem_free=15G,h_vmem=15G,h_fsize=200G
#$ -N step9-findERs-gtex.hippo
#$ -o ./logs/findERs-gtex.txt
#$ -e ./logs/findERs-gtex.txt
#$ -hold_jid pipeline_setup,step5b-meanCoverage-gtex.hippo
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

for meanFile in Coverage/mean*.bw
do
    echo "************************************"
    date
    echo "Initializing script for ${meanFile}"
    echo "************************************"
    Rscript /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/.step9-find_expressed_regions.R -m ${meanFile} -o /dcl01/lieber/ajaffe/lab/brainseq_phase2/gtex/ERs -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode -c 10
done

echo "**** Job ends ****"
date
