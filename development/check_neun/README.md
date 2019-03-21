Check development after adjusting for RNA cell fractions
========================================================


We used the estimated RNA fraction (see [cellComp](../../cellComp) directory) for eight different cell types (fetal replicating, fetal quiescent, OPC, neurons, astrocytes, oligodendrocytes, microglia and endothelial cells). We then re-ran our DEG analysis for each expression feature adjusting for the first 7 cell types (to keep a full rank design model). Across features, we observed a high degree of concordance across our initial analysis and this re-analysis (see table below) at a threshold of Bonferroni<1%. This high concordance can be explained by the observation that the cell fractions are not significantly different (Bonferroni>5%) by age between DLPFC and HIPPO for each of the six age groups we considered, with the only exception being astrocytes for individuals 20 to 50 years old (details in `top_cell_by_agegroup` at [cellComp](../../cellComp), images at the end of this page). Thus while cell RNA fractions change over development, they are not a contributing factor in our development analysis comparing DLPFC and HIPPO.


The tables below show the number of features that are either null in both models, DE in the original model only, DE only in the model where we adjust for the cell type fractions, or DE in both models where DE is assed by a P-bonferroni <1%. In the second table we additionally require that the DE call replicate in BrainSpan by having a P-value <5%.

```R
make_table(tab_pbonf)
#   Null_both Original_only CellProp_only  Both feature        OR pval pval_bonf
# 1      3878          1382          3029 16363    gene  15.15875    0         0
# 2     88040          1049           877  2766      tx 264.70194    0         0
make_table(tab_pbonf_span)
#   Null_both Original_only CellProp_only  Both feature        OR pval pval_bonf
# 1     12501           748          1312 10091    gene  128.5415    0         0
# 2     90783           408           234  1307      tx 1242.8082    0         0
```


The following plot shows the F-statistics from the original model against the NeuN adjusted model with DE determined by P-bonferroni<1% and requiring replication in BrainSpan.

![F-statistics original vs CellProp-adjusted](f_original_vs_f_adjCellProp.png)


## Cell types RNA fraction over age

![Cell over age page 1](../../cellComp/proportions_over_age_by_cell_type_Page_1.png)
![Cell over age page 2](../../cellComp/proportions_over_age_by_cell_type_Page_2.png)
![Cell over age page 3](../../cellComp/proportions_over_age_by_cell_type_Page_3.png)
![Cell over age page 4](../../cellComp/proportions_over_age_by_cell_type_Page_4.png)
![Cell over age page 5](../../cellComp/proportions_over_age_by_cell_type_Page_5.png)
![Cell over age page 6](../../cellComp/proportions_over_age_by_cell_type_Page_6.png)
![Cell over age page 7](../../cellComp/proportions_over_age_by_cell_type_Page_7.png)
![Cell over age page 8](../../cellComp/proportions_over_age_by_cell_type_Page_8.png)
