#!/bin/bash
#$ -cwd
#$ -N featCounts-brainseq.hippo_clean
#$ -o ./logs/featCounts-brainseq_clean.o.txt
#$ -e ./logs/featCounts-brainseq_clean.e.txt
#$ -hold_jid pipeline_setup,featCounts-brainseq.hippo
#$ -m a
echo "**** Job starts ****"
date

echo -e "**** Pipeline version: GitHub sha ****\n3de0c0d9f6fd9a855ca1cad7cba3b746f3ef6f4f"

## Delete temporary files after they have been used
rm -rf /dcl01/lieber/ajaffe/lab/brainseq_phase2/preprocessed_data/Hippo_RiboZero/Counts/junction/tmpdir

echo "**** Job ends ****"
date
