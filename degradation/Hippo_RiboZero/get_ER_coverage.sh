#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=7G,h_fsize=100G
#$ -N step8-bwtool-er-HIPPO_riboZero_degradation.hg38
#$ -o ./logs/bwtool-er-HIPPO_riboZero_degradation.$TASK_ID.txt
#$ -e ./logs/bwtool-er-HIPPO_riboZero_degradation.$TASK_ID.txt
#$ -t 1-12
#$ -tc 100
echo "**** Job starts ****"
date

## get ID
ID=$(cat /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")

## set cut label
CUT=5
mkdir -p /dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Coverage/cut${CUT}

## forward
bedForward=/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/bed/HIPPO_RiboZero_ERs_cut${CUT}_Forward.bed
bwForward=/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Coverage/$ID.Forward.bw
outForward=/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Coverage/cut${CUT}/$ID.Forward.cut${CUT}.txt
bwtool summary $bedForward $bwForward $outForward -header -fill=0 -with-sum

## reverse
bedReverse=/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/bed/HIPPO_RiboZero_ERs_cut${CUT}_Reverse.bed
bwReverse=/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Coverage/$ID.Reverse.bw
outReverse=/dcl01/lieber/ajaffe/lab/brainseq_phase2/degradation/Hippo_RiboZero/Coverage/cut${CUT}/$ID.Reverse.cut${CUT}.txt
bwtool summary $bedReverse $bwReverse $outReverse -header -fill=0 -with-sum
