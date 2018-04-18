#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=200G
#$ -N wgcna_dlpfc_feats
#$ -o logs/wgcna_dlpfc_feats.txt
#$ -e logs/wgcna_dlpfc_feats.txt
#$ -m a
echo "**** Job starts ****"
date

Rscript wgcna_exon_jxn.R

echo "**** Job ends ****"
date

