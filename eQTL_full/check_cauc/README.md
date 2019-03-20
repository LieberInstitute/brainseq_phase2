eQTL CAUC-only sensitivity analysis
===================================

We carried out our eQTL analyses using a mixed population of subjects, most of which are either from AA or CAUC races, while adjusting for quantitative genomic ancestry components. This analysis design is in line with our previous work (Jaffe et al 2018). The statistical power gains by combined analysis, coupled with the overall lack of eQTL effect estimates by race, motivated these combined analyses here.

We however ran a sensitivity analysis by subsetting each of our eQTL models (DLPFC, HIPPO and interaction between the two brain regions) to just the self-reported CAUC samples (43.8 to 46.1% of the original sample size, see table below) at the gene expression feature level. We still adjusted for quantitative genomic ancestry components to account for differences between self-reported and genomic ancestry, as well as to further control for samples being genotyped across multiple microarray platforms.

```{r}
sample_sizes
#     n  AA CAUC HISP AS CAUC_percent       model
# 1 792 417  356   10  9     44.94949 interaction
# 2 397 204  174   10  9     43.82872       DLPFC
# 3 395 213  182    0  0     46.07595       HIPPO
```

Each of the three eQTL models showed high directionality agreement between the genotype regression coefficient from the significant mixed ethnicity eQTL SNP-gene pairs and those from the CAUC-only analysis: interaction: 94.5%, DLPFC: 94.1%, HIPPO: 94.5%. While allelic directionally was highly concordant between ethnicities among those eQTLs identified in combined analysis, many of these eQTLs are much less significant within CAUCs only. For example, while between 82-87% of mixed-ethnicity eQTLs were concordant and at least marginally significant in CAUC at p < 0.05 (interaction 87.2%, DLPFC: 82.3%, HIPPO: 83.1%), many fewer eQTLs remain significant at more stringent p-values, like p < 0.001 (interaction: 57.9%, DLPFC: 57.5%, HIPPO: 59.1%). Given the high directional consistency (~95%), we believe this decrease in identified eQTLs in CAUC-only is likely due to the decrease in power from using a smaller number of individuals. This is further supported by the vast majority of the SNP-gene pairs in CAUCs at FDR < 0.01 were themselves genome-wide significant (FDR < 0.01) in our original mixed ancestry analyses (interaction 78.4%, DLPFC: 85.2%, HIPPO: 84.8%, see table below) with identical directionality for those observed in both analyses (interaction 100%, DLPFC: 99.99%, HIPPO: 99.99%).

```{r}
cauc_sig_in_brainseq_tab
#             N_FALSE N_TRUE Percent_FALSE Percent_TRUE
# interaction    2393   8669      21.63262     78.36738
# dlpfc        124627 720099      14.75354     85.24646
# hippo         87323 485973      15.23175     84.76825
```
