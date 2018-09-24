#!/bin/bash
#$ -cwd
#$ -N featCounts-brainseq.dlpfc_clean
#$ -o ./logs/featCounts-brainseq_clean.o.txt
#$ -e ./logs/featCounts-brainseq_clean.e.txt
#$ -hold_jid pipeline_setup,featCounts-brainseq.dlpfc
#$ -m a
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\na2b1cda9a67211320865c002ee1a20429cce32b1"

## Delete temporary files after they have been used
rm -rf /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/DLPFC_RiboZero/Counts/junction/tmpdir

echo "**** Job ends ****"
date
