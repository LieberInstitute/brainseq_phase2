#!/bin/bash

## Usage:
# sh run_explore_reg_specific_top.sh

mkdir -p logs
mkdir -p rda
mkdir -p pdf

for featuretype in gene exon jxn tx
do
    for agegroup in adult fetal
    do

SHORT="explore_reg_specific_${featuretype}_${agegroup}_top"

# Construct shell file
echo "Creating script for feature type ${featuretype} and age group ${agegroup}"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=25G,h_vmem=25G,h_fsize=100G
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

Rscript explore_reg_specific_top.R -t ${featuretype} -a ${agegroup}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
    done
done
