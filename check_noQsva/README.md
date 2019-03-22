SCZD vs non-psychiatric controls without qSVA
=============================================


## With vs without qSVA

The following plots compare the SCZD vs control t-statistics for each region with and without qSVA adjustment at the gene expression level. Points in pink are significant (FDR<5%) in only one of the analyses, while those in blue are significant (FDR<5%) in both.

| HIPPO  | DLPFC  |
|---|---|
| ![HIPPO no qSVA vs HIPPO with qSVA](pdf/scatter_models_Page_6.png)  | ![DLPFC no qSVA vs DLPFC with qSVA](pdf/scatter_models_Page_7.png)  |

## HIPPO vs DLPFC

The following plots show the comparisons between HIPPO and DLPFC either adjusting for quality surrogate variables (qSVs) or without at the gene expression level.

| without qSVA  | with qSVA  |
|---|---|
| ![HIPPO no qSVA vs DLPFC no qSVA](pdf/scatter_models_Page_1.png)  | ![HIPPO with qSVA vs DLPFC no qSVA](pdf/original/scatter_models_Page_1.png)  |


## Gene set enrichment

The following plots show the results from the gene set enrichment analysis for biological processes (BP), molecular functions (MF), cellular components (CC) and KEGG pathways with and without adjusting for qSVs.

### DLPFC

| without qSVA  | with qSVA  |
|---|---|
| ![Enriched BP without qSVA](pdf/gse_dlpfc_Page_1.png)  | ![Enriched BP with qSVA](pdf/original/gse_dlpfc_Page_1.png)  |
| ![Enriched MF without qSVA](pdf/gse_dlpfc_Page_2.png)  |  ![Enriched MF with qSVA](pdf/original/gse_dlpfc_Page_2.png) |
| ![Enriched CC without qSVA](pdf/gse_dlpfc_Page_3.png)  |  ![Enriched CC with qSVA](pdf/original/gse_dlpfc_Page_3.png) |
| ![Enriched KEGG pathways without qSVA](pdf/gse_dlpfc_Page_4.png) |  ![Enriched KEGG pathways with qSVA](pdf/original/gse_dlpfc_Page_4.png) |


### HIPPO

| without qSVA  | with qSVA  |
|---|---|
| ![Enriched BP without qSVA](pdf/gse_hippo_Page_1.png)  | ![Enriched BP with qSVA](pdf/original/gse_hippo_Page_1.png)  |
| ![Enriched MF without qSVA](pdf/gse_hippo_Page_2.png)  |  ![Enriched MF with qSVA](pdf/original/gse_hippo_Page_2.png) |
| ![Enriched CC without qSVA](pdf/gse_hippo_Page_3.png)  |  ![Enriched CC with qSVA](pdf/original/gse_hippo_Page_3.png) |
| ![Enriched KEGG pathways without qSVA](pdf/gse_hippo_Page_4.png) |  ![Enriched KEGG pathways with qSVA](pdf/original/gse_hippo_Page_4.png) |

## Venn genes DLFPC FDR<10% and HIPPO FDR<20%

Venn diagrams of the genes at FDR<10% for DLPFC and FDR<20% for HIPPO, with and without qSVA adjustment. First, for genes with higher expression in controls than SCZD, then the reverse.

| without qSVA  | with qSVA  |
|---|---|
| ![Control > SCZD without qSVA](pdf/venn_de_genes_by_sign_Page_02.png)  | ![Control > SCZD with qSVA](pdf/original/venn_de_genes_by_sign_Page_02.png)  |
| ![SCZD > Control without qSVA](pdf/venn_de_genes_by_sign_Page_03.png)  |  ![SCZD > Control with qSVA](pdf/original/venn_de_genes_by_sign_Page_03.png) |


## Enriched GOs

### DLFPC FDR<10% and HIPPO FDR<20%

The following plots show the top enriched gene ontology terms with the DLPFC FDR<10% and the HIPPO FDR<20% features.

| without qSVA  | with qSVA  |
|---|---|
| ![Enriched BP without qSVA](pdf/go_de_genes_Page_1.png)  | ![Enriched BP with qSVA](pdf/original/go_de_genes_Page_1.png)  |
| ![Enriched MF without qSVA](pdf/go_de_genes_Page_2.png)  |  ![Enriched MF with qSVA](pdf/original/go_de_genes_Page_2.png) |
| ![Enriched CC without qSVA](pdf/go_de_genes_Page_3.png)  |  ![Enriched CC with qSVA](pdf/original/go_de_genes_Page_3.png) |
| ![Enriched KEGG pathways without qSVA](pdf/go_de_genes_Page_4.png) |  ![Enriched KEGG pathways with qSVA](pdf/original/go_de_genes_Page_4.png) |


### Top 200 features

The following plots show the top enriched gene ontology terms with the top 200 features per brain region. This is more comparable across the sets with and without qSVA since there are more genes with FDR<10% in DLPFC and FDR<20% in HIPPO when you don't adjust for qSVs.

| without qSVA  | with qSVA  |
|---|---|
| ![Enriched BP without qSVA](pdf/go_de_genes_top200_Page_1.png)  | ![Enriched BP with qSVA](pdf/original/go_de_genes_top200_Page_1.png)  |
| ![Enriched MF without qSVA](pdf/go_de_genes_top200_Page_2.png)  |  ![Enriched MF with qSVA](pdf/original/go_de_genes_top200_Page_2.png) |
| ![Enriched CC without qSVA](pdf/go_de_genes_top200_Page_3.png)  |  ![Enriched CC with qSVA](pdf/original/go_de_genes_top200_Page_3.png) |
| ![Enriched KEGG pathways without qSVA](pdf/go_de_genes_top200_Page_4.png) |  ![Enriched KEGG pathways with qSVA](pdf/original/go_de_genes_top200_Page_4.png) |


## BSP1 and CMC

These plots show how each brain region compares against BrainSeq Phase I DLPFC polyA+ data (BSP1) or the CommonMind Consortium (CMC) data, with and without adjusting for qSVs. These plots are at the gene expression level.

| without qSVA  | with qSVA  |
|---|---|
| ![HIPPO no qSVA vs BSP1](pdf/scatter_models_Page_2.png)  | ![HIPPO with qSVA vs BSP1](pdf/original/scatter_models_Page_2.png)  |
| ![HIPPO no qSVA vs CMC](pdf/scatter_models_Page_3.png)  |  ![HIPPO with qSVA vs CMC](pdf/original/scatter_models_Page_3.png) |
| ![DLPFC no qSVA vs BSP1](pdf/scatter_models_Page_4.png)  |  ![DLPFC with qSVA vs BSP1](pdf/original/scatter_models_Page_4.png) |
| ![DLPFC no qSVA vs CMC](pdf/scatter_models_Page_5.png)  |  ![DLPFC with qSVA vs CMC](pdf/original/scatter_models_Page_5.png) |



