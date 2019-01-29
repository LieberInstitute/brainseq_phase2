#!/bin/bash

## Usage:
# sh apply_weights.sh

mkdir -p logs

for region in HIPPO DLPFC
do
    
    # for feature in gene exon jxn tx
    for feature in gene
    do
        
        # set of summary stats
        for summstats in pgc2 psycm
        do
        
            SHORT="apply_weights_${region}_${feature}_${summstats}"

            # Construct shell file
            echo "Creating script for chromosome ${region} at the ${feature} level for ${summstats}"
            cat > .${SHORT}.sh <<EOF

#!/bin/bash
#$ -cwd
#$ -l mem_free=30G,h_vmem=30G,h_fsize=100G
#$ -N ${SHORT}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"

## Load dependencies
module load fusion_twas/github

## List current modules
module list

## Choose the correct GWAS summary statistics file
if [ "${summstats}" == "psycm" ]
then
    summstatsfile="/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/psycm/clozuk_pgc2.meta.reformatted.sumstats_hg38_ourname"
elif [ "${summstats}" == "pgc2" ]
then
    summstatsfile="/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/pgc_scz2_sumstats/PGC2.SCZ.sumstats_hg38_ourname"
else
    echo "Unexpected ${summstats} input"
fi

## Apply weights for the given region/feature pair and the given GWAS summary statistics
cd ${region}/${feature}
mkdir -p ${summstats}
cd ${summstats}


for chr in {1..22}
do
    echo "*************************"
    echo ""
    echo "processing chromosome \${chr}"
    date
    echo ""

## Create summarized analysis
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.assoc_test.R \
    --sumstats \${summstatsfile} \
    --weights /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/${region}/${feature}/${region}_${feature}.pos \
    --weights_dir /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/${region}/${feature} \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr ${chr} \
    --out ${summstats}.\${chr}.dat
    
    echo ""
    echo "making plots for chromosome \${chr}"
    date
    echo ""

## companion plotting step
Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.post_process.R \
    --sumstats /\${summstatsfile} \
    --input ${summstats}.\${chr}.dat \
    --out ${summstats}.\${chr}.analysis \
    --ref_ld_chr /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/reference_hg38/LDREF_hg38/1000G.EUR. \
    --chr ${chr} \
    --plot --locus_win 100000 --verbose 2 --plot_individual --plot_eqtl --plot_corr \
    --glist_path "/jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/glist-hg38"

done

echo "**** Job ends ****"
date
EOF
            call="qsub .${SHORT}.sh"
            echo $call
            $call
        done
    done
done