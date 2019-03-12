#!/bin/bash

## Usage:
# sh run_limma_reg_specific_subsample.sh

mkdir -p subsample
mkdir -p subsample/logs
mkdir -p subsample/rda
mkdir -p subsample/pdf

for featuretype in gene #exon jxn tx
do
    for agegroup in adult #fetal
    do

SHORT="limma_reg_specific_${featuretype}_${agegroup}_subsample"

# Construct shell file
echo "Creating script for feature type ${featuretype} and age group ${agegroup}"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N ${SHORT}
#$ -o ./subsample/logs/${SHORT}.\$TASK_ID.txt
#$ -e ./subsample/logs/${SHORT}.\$TASK_ID.txt
#$ -t 1-100
#$ -m a

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"

# module load conda_R/3.4.x
Rscript limma_reg_specific_subsample.R -t ${featuretype} -a ${agegroup} -i \${SGE_TASK_ID}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
    done
done
