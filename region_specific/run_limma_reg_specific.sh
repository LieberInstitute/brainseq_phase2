#!/bin/bash

## Usage:
# sh run_limma_reg_specific.sh

mkdir -p logs

for featuretype in gene exon jxn tx
do
    for agegroup in adult fetal
    do

SHORT="limma_reg_specific_${featuretype}_${agegroup}"

# Construct shell file
echo "Creating script for feature type ${featuretype} and age group ${age}"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=80G,h_vmem=80G,h_fsize=100G
#$ -N ${SHORT}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

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
