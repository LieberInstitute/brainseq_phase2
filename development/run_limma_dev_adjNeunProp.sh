#!/bin/bash

## Usage:
# sh run_limma_dev.sh

mkdir -p logs
mkdir -p rda
mkdir -p pdf

for featuretype in gene exon jxn tx
do

SHORT="limma_dev_${featuretype}"

# Construct shell file
echo "Creating script for feature type ${featuretype}"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l mem_free=45G,h_vmem=45G,h_fsize=100G
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

module load conda_R/3.4.x
Rscript limma_dev.R -t ${featuretype}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
done
