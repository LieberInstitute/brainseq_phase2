#!/bin/bash
#$ -cwd
#$ -N step4-featCounts-HIPPO_riboZero_degradation.hg38_clean
#$ -o ./logs/featCounts-HIPPO_riboZero_degradation_clean.txt
#$ -e ./logs/featCounts-HIPPO_riboZero_degradation_clean.txt
#$ -hold_jid pipeline_setup,step4-featCounts-HIPPO_riboZero_degradation.hg38
#$ -m a
echo "**** Job starts ****"
date


## Delete temporary files after they have been used
rm -rf /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Counts/junction/tmpdir

echo "**** Job ends ****"
date
