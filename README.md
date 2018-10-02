BrainSeq Phase II analysis
==========================

This repository contains analysis code for the BrainSeq Phase II project from the BrainSeq Consortium carried out by researchers at the [Lieber Institute for Brain Development](https://www.libd.org/).

If you wish to visualize the eQTL results described in this project, please use the [LIBD eQTL browser](http://eqtl.brainseq.org/).

## License

<img src="https://licensebuttons.net/l/by-nc/3.0/88x31.png" alt width="88" height="31" scale="0">
Attribution-NonCommercial: CC BY-NC

This license lets others remix, tweak, and build upon our work non-commercially as long as they acknowledge our work.

[View License Deed](https://creativecommons.org/licenses/by-nc/4.0) | [View Legal Code](https://creativecommons.org/licenses/by-nc/4.0/legalcode)

## Citation

If you use anything in this repository please cite the following publication:

Collado-Torres L, Burke EE, Peterson A, Shin JH, Straub RE, Rajpurohit A, Semick SA, Ulrich WS, BrainSeq Consortium, Valencia C, Tao R, Deep-Soboslay A, Hyde TM, Kleinman JE, Weinberger DR, Jaffe AE. Regional heterogeneity in gene expression, regulation and coherence in hippocampus  and dorsolateral prefrontal cortex across development and in schizophrenia. bioRxiv, DOI 10.1101/426213. 2018.

Pre-print URL: https://www.biorxiv.org/content/early/2018/09/26/426213.

## Files

| directory | contents |
| --------- | -------- |
| [`brainspan`](brainspan/) | Code for processing the BrainSpan data. |
| [`browser`](browser/) | Code for creating the files for the eQTL browser. Contains a detailed README file. |
| [`caseControl`](caseControl/) | Initial (un-used) exploratory code for the SCZD case-control analysis. Final code is at the [`qsva_brain`](https://github.com/LieberInstitute/qsva_brain) repository. |
| [`caseControl_HIPPO_checks`](caseControl_HIPPO_checks/) | Code for checking the SCZD case-control HIPPO results. Final code at [`qsva_brain`](https://github.com/LieberInstitute/qsva_brain) repo. |
| [`casecontrolint`](casecontrolint/) | Code for the brain region and SCZD diagnosis status interaction DE analysis. |
| [`correlation`](correlation/) | DLPFC and HIPPO expression correlation analyses. |
| [`demographics`](demographics/) | Code for exploring demographic variables such as RIN. |
| [`development`](development/) | Code for the DE analysis across age using a linear spline model. |
| [`eQTL_GWAS_riskSNPs`](eQTL_GWAS_riskSNPs/) | eQTL analysis using PGC2 GWAS risk SNPs and neighboring SNPs. |
| [`eQTL_full`](eQTL_full/) | Genome wide eQTL analyses. |
| [`eQTL_full_GTEx`](eQTL_full_GTEx/) | Replication eQTL analyses using GTEx data. |
| [`expr_cutoff`](expr_cutoff/) | Code for filtering the features with low expression values and creating the RSE objects used throughout the project. |
| [`gtex`](gtex/) | Code for processing the HIPPO GTEx data and preparing the genotype data for the eQTL analysis. |
| [`gtex_both`](gtex_both/) | Code for merging the DLPFC and HIPPO GTEx data. |
| [`gtex_dlpfc`](gtex_dlpfc/) | Code for processing the HIPPO GTEx data. |
| [`preprocessed_data`](preprocessed_data/) | Code for processing the RNA-seq reads. Uses the LIBD RNA-seq pipeline developed by EE Burke, L Collado-Torres, and AE Jaffe. |
| [`region_specific`](region_specific/) | Code for the DE analyses between HIPPO and DLPFC for prenatal and adult controls. |
| [`supp_tabs`](supp_tabs/) | Code for creating the supplementary tables. |
| [`misc`](misc/) | Files from early explorations including checking for some sample swaps and quality checks. |

## LIBD internal:

JHPCE location: `/dcl01/lieber/ajaffe/lab/brainseq_phase2`
