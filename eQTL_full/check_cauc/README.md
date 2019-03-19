eQTL CAUC-only sensitivity analysis
===================================

We carried out our eQTL analyses using a mixed population of subjects, most of which are either from AA or CAUC races, while adjusting for quantitative genomic ancestry components. This analysis design is in line with our previous work (Jaffe et al 2018). The statistical power gains by combined analysis, coupled with the overall lack of eQTL effect estimates by race, motivated these combined analyses here.

We however ran a sensitivity analysis by subsetting each of our eQTL models (DLPFC, HIPPO and interaction between the two brain regions) to just the CAUC samples (43.8 to 46.1% of the original sample size, see table below) at the gene expression feature level. We still adjusted for quantitative genomic ancestry components despite using individuals from the same race given their heterogeneous ancestry.

```{r}
sample_sizes
#     n  AA CAUC HISP AS CAUC_percent       model
# 1 792 417  356   10  9     44.94949 interaction
# 2 397 204  174   10  9     43.82872       DLPFC
# 3 395 213  182    0  0     46.07595       HIPPO
```

Each of the three models showed high directionality agreement (interaction: 94.5%, DLPFC: 94.1%, HIPPO: 94.5%) on the genotype regression coefficient for the originally identified eQTL SNP-gene pairs that we identified at FDR < 0.01. When considering both directionality and a nominal P-value < 0.05 in the CAUC re-analysis, we observed a decreased agreement (interaction 87.2%, DLPFC: 82.3%, HIPPO: 83.1%). This decrease is likely due to the decrease in power from using only CAUC individuals. 

The great majority of the SNP-gene pairs at FDR < 0.01 we identified in the CAUC-only re-analysis were observed at FDR < 0.01 in our original analysis (interaction 78.4%, DLPFC: 85.2%, HIPPO: 84.8%, see table below) with nearly identical directionality for those observed in both analyses (interaction 100%, DLPFC: 99.99%, HIPPO: 99.99%).

```{r}
cauc_sig_in_brainseq_tab
#             N_FALSE N_TRUE Percent_FALSE Percent_TRUE
# interaction    2393   8669      21.63262     78.36738
# dlpfc        124627 720099      14.75354     85.24646
# hippo         87323 485973      15.23175     84.76825
```
