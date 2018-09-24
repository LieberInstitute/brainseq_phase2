BrainSeq Phase 2 analysis
=========================

This repository contains analysis code for the BrainSeq Phase II project from the BrainSeq Consortium carried out by researchers at the [Lieber Institute for Brain Development](https://www.libd.org/).

If you wish to visualize the eQTL results described in this project, please use the [LIBD eQTL browser](http://eqtl.brainseq.org/).

## Citation

If you use anything in this repository please cite the following publication:

Collado-Torres L, Burke EE, Peterson A, Shin JH, Straub RE, Rajpurohit A, Semick SA, Ulrich WS, BrainSeq Consortium, Valencia C, Tao R, Deep-Soboslay A, Hyde TM, Kleinman JE, Weinberger DR, Jaffe AE. Regional heterogeneity in gene expression, regulation and coherence in hippocampus  and dorsolateral prefrontal cortex across development and in schizophrenia. JOURNAL + DOI. 2018.

## Files

| directory | contents |
| --------- | -------- |
| [`brainspan`](brainspan/) | Code for processing the BrainSpan data. |
| [`browser`](browser/) | Code for creating the files for the eQTL browser. Contains a detailed README file. |
| [`caseControl`](caseControl/) | Initial code for the SCZD case-control analysis. Final code is at the [`qsva_brain`](https://github.com/LieberInstitute/qsva_brain) repository. |
| [`caseControl_HIPPO_checks`](caseControl_HIPPO_checks/) | Code for checking the SCZD case-control HIPPO results. Final code at [`qsva_brain`](https://github.com/LieberInstitute/qsva_brain) repo. |
| [`casecontrolint`](casecontrolint/) | Code for the brain region and SCZD diagnosis status interaction DE analysis. |
| [`cellComp`](cellComp/) | Un-used cell composition code. |
| [`correlation`](correlation/) | DLPFC and HIPPO expression correlation analyses. |
| [`degradation`](degradation/) | Un-used degradation files. Final degradation code at [`qsva_brain`](https://github.com/LieberInstitute/qsva_brain) repo. |
| [`degradation_strand_minus`](degradation_strand_minus/) | Un-used degradation files. Final degradation code at [`qsva_brain`](https://github.com/LieberInstitute/qsva_brain) repo. |
| [`degradation_strand_positive`](degradation_strand_positive/) | Un-used degradation files. Final degradation code at [`qsva_brain`](https://github.com/LieberInstitute/qsva_brain) repo. |
| [`demographics`](demographics/) | Code for exploring demographic variables such as RIN. |
| [`development`](development/) | Code for the DE analysis across age using a linear spline model. |
| [`eQTL_DNAm_mediation`](eQTL_DNAm_mediation/) | Un-used DNA mediation files. |
| [`eQTL_GWAS_riskSNPs`](eQTL_GWAS_riskSNPs/) | eQTL analysis using PGC2 GWAS risk SNPs and neighboring SNPs. |
| [`eQTL_full`](eQTL_full/) | Genome wide eQTL analyses. |
| [`eQTL_full_GTEx`](eQTL_full_GTEx/) | Replication eQTL analyses using GTEx data. |
| [`expr_cutoff`](expr_cutoff/) | Code for filtering the features with low expression values and creating the RSE objects used throughout the project. |
| [`gtex`](gtex/) | Code for processing the HIPPO GTEx data and preparing the genotype data for the eQTL analysis. |
| [`gtex_both`](gtex_both/) | Code for merging the DLPFC and HIPPO GTEx data. |
| [`gtex_dlpfc`](gtex_dlpfc/) | Code for processing the HIPPO GTEx data. |
| [`region_specific`](region_specific/) | Code for the DE analyses between HIPPO and DLPFC for prenatal and adult controls. |
| [`supp_tabs`](supp_tabs/) | Code for creating the supplementary tables. |
| [`wgcna`](wgcna/) | Un-used WGCNA code for DLPFC and HIPPO. |
| [`wgcna_combined`](wgcna_combined/) | Un-used WGCNA code for a combined analysis using DLPFC and HIPPO. |
| [`misc`](misc/) | Files from early explorations including checking for some sample swaps and quality checks. |

## LIBD internal:

JHPCE location: `/dcl01/lieber/ajaffe/lab/brainseq_phase2`
