#!/bin/bash

## Usage:
# sh run_limma_reg_specific.sh

mkdir -p logs
mkdir -p rda
mkdir -p pdf

for featuretype in gene exon jxn tx
do
    for agegroup in adult fetal
    do

SHORT="limma_reg_specific_${featuretype}_${agegroup}"

# Construct shell file
echo "Creating script for feature type ${featuretype} and age group ${agegroup}"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=40G,h_vmem=40G,h_fsize=100G
#$ -N ${SHORT}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

module load conda_R/3.4.x
Rscript limma_reg_specific.R -t ${featuretype} -a ${agegroup}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
    done
done
