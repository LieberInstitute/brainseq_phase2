#!/bin/bash

## Usage:
# sh run_limma_dev.sh

mkdir -p logs

for featuretype in gene exon jxn
do

SHORT="limma_dev_${featuretype}"

# Construct shell file
echo "Creating script for feature type ${featuretype}"

cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=150G,h_vmem=150G,h_fsize=100G
#$ -N ${SHORT}
#$ -o ./logs/${SHORT}.txt
#$ -e ./logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

module load conda_R/3.4.x
Rscript limma_dev.R -t ${featuretype}

echo "**** Job ends ****"
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
done
