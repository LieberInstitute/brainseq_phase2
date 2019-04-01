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
| [`bsp1`](bsp1/) | eQTL replication with BrainSeq Phase I DLPFC polyA+ data. See this [README](https://github.com/LieberInstitute/brainseq_phase2/tree/master/bsp1/eqtl/full) for the results. |
| [`caseControl`](caseControl/) | Initial (un-used) exploratory code for the SCZD case-control analysis. Final code is at the [`qsva_brain`](https://github.com/LieberInstitute/qsva_brain) repository. |
| [`caseControl_HIPPO_checks`](caseControl_HIPPO_checks/) | Code for checking the SCZD case-control HIPPO results. Final code at [`qsva_brain`](https://github.com/LieberInstitute/qsva_brain) repo. |
| [`casecontrolint`](casecontrolint/) | Code for the brain region and SCZD diagnosis status interaction DE analysis. |
| [`cellComp`](cellComp/) | RNA cell fraction deconvolution. Contains a detailed README file. |
| [`check_expr`](check_expr/) | Check genes not expressed at other feature levels. Contains a detailed README file. |
| [`check_noQsva`](check_noQsva/) | Check gene-level SCZD vs control DEG results without adjusting for qSVs. Contains a detailed README file. |
| [`check_protein_coding`](check_protein_coding/) | Check for protein coding and non-coding enrichment/depletion. Contains a detailed README file. |
| [`check_sex/casecontrol`](check_sex/casecontrol/) | Check SCZD case control results by sex. Contains a detailed README file. |
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
| [`supp_tabs`](supp_tabs/) | Code for creating some supplementary tables. |
| [`twas`](twas/) | Perform TWAS analysis using the FUSION TWAS software by Gusev et al., Nature Genetics, 2016. Contains a detailed README file. |
| [`misc`](misc/) | Files from early explorations including checking for some sample swaps and quality checks. |
| [`KCNQ1_snp_check`](KCNQ1_snp_check/) | Random check. |


## Public Files

| Description | JHPCE path | URL |
| --- | --- | --- |
| Table S15. Genome wide significant eQTL snp-feature pairs | `/dcl01/lieber/ajaffe/lab/brainseq_phase2/supp_tabs/SupplementaryTableXX_eQTL.tar.gz` | [AWS](https://s3.us-east-2.amazonaws.com/libd-brainseq2/SupplementaryTable15_eQTL.tar.gz) |
| Unfiltered gene [RangedSummarizedExperiment](http://bioconductor.org/packages/SummarizedExperiment) object | `/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_gene_unfiltered.Rdata` | [AWS](https://s3.us-east-2.amazonaws.com/libd-brainseq2/rse_gene_unfiltered.Rdata) |
| Unfiltered exon [RangedSummarizedExperiment](http://bioconductor.org/packages/SummarizedExperiment) object | `/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_exon_unfiltered.Rdata` | [AWS](https://s3.us-east-2.amazonaws.com/libd-brainseq2/rse_exon_unfiltered.Rdata) |
| Unfiltered exon-exon junction [RangedSummarizedExperiment](http://bioconductor.org/packages/SummarizedExperiment) object | `/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_jxn_unfiltered.Rdata` | [AWS](https://s3.us-east-2.amazonaws.com/libd-brainseq2/rse_jxn_unfiltered.Rdata) |
| Unfiltered transcript [RangedSummarizedExperiment](http://bioconductor.org/packages/SummarizedExperiment) object | `/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/unfiltered/rse_tx_unfiltered.Rdata` | [AWS](https://s3.us-east-2.amazonaws.com/libd-brainseq2/rse_tx_unfiltered.Rdata) |
| DLPFC vs HIPPO DEG objects (adult and prenatal, includes BrainSpan replication and cell RNA fraction sensitivity results)| `/dcl01/lieber/ajaffe/lab/brainseq_phase2/region_specific/rda/RegionSpecificDEGobjects.tar.gz` | [AWS](https://s3.us-east-2.amazonaws.com/libd-brainseq2/RegionSpecificDEGobjects.tar.gz) |
| Development DEG objects (includes BrainSpan replication and cell RNA fraction sensitivity results) | `/dcl01/lieber/ajaffe/lab/brainseq_phase2/development/rda/DevelopmentDEGobjects.tar.gz` | [AWS](https://s3.us-east-2.amazonaws.com/libd-brainseq2/DevelopmentDEGobjects.tar.gz) |
| SCZD vs non-psychiatric control DEG objects (includes qSVs as well as results for the interaction and no-qSVA gene-level sensitivity analyses). See [BrainSeq Phase I SCZD DE features](https://s3.us-east-2.amazonaws.com/jaffe-nat-neuro-2018/expressed_de_features.rda) for more. | `/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/SCZDvsControlDEGobjects.tar.gz` | [AWS](https://s3.us-east-2.amazonaws.com/libd-brainseq2/) |
| Demographic table including cell RNA fraction estimates | `/dcl01/lieber/ajaffe/lab/brainseq_phase2/cellComp/methprop_pd.Rdata` | [AWS](https://s3.us-east-2.amazonaws.com/libd-brainseq2/methprop_pd.Rdata) |
| TWAS DLPFC weights | `/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/DLPFC/DLPFC_weights.tar.gz` | [AWS](https://s3.us-east-2.amazonaws.com/libd-brainseq2/DLPFC_weights.tar.gz) |
| TWAS HIPPO weights | `/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/HIPPO/HIPPO_weights.tar.gz` | [AWS](https://s3.us-east-2.amazonaws.com/libd-brainseq2/HIPPO_weights.tar.gz) |
| TWAS results R objects | `/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/rda/TWAS_results.tar.gz` | [AWS](https://s3.us-east-2.amazonaws.com/libd-brainseq2/TWAS_results.tar.gz) |

Check the corresponding scripts or use GitHub's search feature to find where each of the R objects were created. If you have questions about the files, please open an [issue](https://github.com/LieberInstitute/brainseq_phase2/issues).

## LIBD internal:

JHPCE location: `/dcl01/lieber/ajaffe/lab/brainseq_phase2`
