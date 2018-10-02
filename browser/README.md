eQTL browser files
------------------

Location at JHPCE: `/dcl01/lieber/ajaffe/lab/brainseq_phase2/browser`.

This directory has all the files for creating the eQTL browser. As stated in the [log file](logs/extract_data.txt) the files and their number of lines are:

```
      396584 BrainSeqPhaseII_clean_expression_development_exon.txt
       24653 BrainSeqPhaseII_clean_expression_development_gene.txt
      297182 BrainSeqPhaseII_clean_expression_development_jxn.txt
       92733 BrainSeqPhaseII_clean_expression_development_tx.txt
      396584 BrainSeqPhaseII_clean_expression_eqtl_dlpfc_exon.txt
       24653 BrainSeqPhaseII_clean_expression_eqtl_dlpfc_gene.txt
      297182 BrainSeqPhaseII_clean_expression_eqtl_dlpfc_jxn.txt
       92733 BrainSeqPhaseII_clean_expression_eqtl_dlpfc_tx.txt
      396584 BrainSeqPhaseII_clean_expression_eqtl_hippo_exon.txt
       24653 BrainSeqPhaseII_clean_expression_eqtl_hippo_gene.txt
      297182 BrainSeqPhaseII_clean_expression_eqtl_hippo_jxn.txt
       92733 BrainSeqPhaseII_clean_expression_eqtl_hippo_tx.txt
      396584 BrainSeqPhaseII_clean_expression_eqtl_interaction_exon.txt
       24653 BrainSeqPhaseII_clean_expression_eqtl_interaction_gene.txt
      297182 BrainSeqPhaseII_clean_expression_eqtl_interaction_jxn.txt
       92733 BrainSeqPhaseII_clean_expression_eqtl_interaction_tx.txt
      396584 BrainSeqPhaseII_clean_expression_regionspecific_adult_exon.txt
       24653 BrainSeqPhaseII_clean_expression_regionspecific_adult_gene.txt
      297182 BrainSeqPhaseII_clean_expression_regionspecific_adult_jxn.txt
       92733 BrainSeqPhaseII_clean_expression_regionspecific_adult_tx.txt
      396584 BrainSeqPhaseII_clean_expression_regionspecific_prenatal_exon.txt
       24653 BrainSeqPhaseII_clean_expression_regionspecific_prenatal_gene.txt
      297182 BrainSeqPhaseII_clean_expression_regionspecific_prenatal_jxn.txt
       92733 BrainSeqPhaseII_clean_expression_regionspecific_prenatal_tx.txt
      396584 BrainSeqPhaseII_clean_expression_sczd_casecontrol_dlpfc_exon.txt
       24653 BrainSeqPhaseII_clean_expression_sczd_casecontrol_dlpfc_gene.txt
      297182 BrainSeqPhaseII_clean_expression_sczd_casecontrol_dlpfc_jxn.txt
       92733 BrainSeqPhaseII_clean_expression_sczd_casecontrol_dlpfc_tx.txt
      396584 BrainSeqPhaseII_clean_expression_sczd_casecontrol_hippo_exon.txt
       24653 BrainSeqPhaseII_clean_expression_sczd_casecontrol_hippo_gene.txt
      297182 BrainSeqPhaseII_clean_expression_sczd_casecontrol_hippo_jxn.txt
       92733 BrainSeqPhaseII_clean_expression_sczd_casecontrol_hippo_tx.txt
      396584 BrainSeqPhaseII_clean_expression_sczd_casecontrol_interaction_exon.txt
       24653 BrainSeqPhaseII_clean_expression_sczd_casecontrol_interaction_gene.txt
      297182 BrainSeqPhaseII_clean_expression_sczd_casecontrol_interaction_jxn.txt
       92733 BrainSeqPhaseII_clean_expression_sczd_casecontrol_interaction_tx.txt
    28863540 BrainSeqPhaseII_eQTL_dlpfc_full.txt
     3698269 BrainSeqPhaseII_eQTL_dlpfc_raggr.txt
   758634578 BrainSeqPhaseII_eQTL_dlpfc_replication_GTEx.txt
    22603132 BrainSeqPhaseII_eQTL_hippo_full.txt
     3698269 BrainSeqPhaseII_eQTL_hippo_raggr.txt
   758634578 BrainSeqPhaseII_eQTL_hippo_replication_GTEx.txt
    22603132 BrainSeqPhaseII_eQTL_interaction_full.txt
   398556091 BrainSeqPhaseII_eQTL_interaction_replication_GTEx.txt
      396584 BrainSeqPhaseII_feature_annotation_exon.txt
       24653 BrainSeqPhaseII_feature_annotation_gene.txt
      297182 BrainSeqPhaseII_feature_annotation_jxn.txt
       92733 BrainSeqPhaseII_feature_annotation_tx.txt
         902 BrainSeqPhaseII_sample_metadata.txt
     7023861 BrainSeqPhaseII_snp_annotation.txt
     7023861 BrainSeqPhaseII_snp_genotype.txt
```

### Sample information

```
         902 BrainSeqPhaseII_sample_metadata.txt
```

Includes many columns, but some of the more useful ones are:

* `SAMPLE_ID`: actual file ID
* `RNum`: internal LIBD RNA sequencing ID
* `BrNum`: internal LIBD brain (subject) ID
* `Region`: either HIPPO or DLPFC
* `Dx`: either Control or SCZD
* `Age`: numeric age in years
* `DLFPC_HIPPO_correlation_*` (`*` can be `geneRpkm`, `exonRpkm`, `jxnRp10m` or `txTpm`): HIPPO vs DLPFC correlation across the corresponding feature (only 265 subjects). Values are stored only for the `Region == DLPFC` entries.
* `analysis_regionspecific_adult`: logical, TRUE if the sample was used for the HIPPO vs DLPFC region analysis with adult controls.
* `analysis_regionspecific_prenatal`: similar to the previous one, but for prenatal age.
* `analysis_development`: logical, TRUE if the same was used for the analysis across development (age); they are only controls.
* `analysis_sczd_casecontrol_dlpfc`: logical, TRUE if used for the SCZD vs control DE analysis for the DLPFC region.
* `analysis_sczd_casecontrol_hippo`: logical, TRUE if used for the SCZD vs control DE analysis for the HIPPO region.
* `analysis_sczd_casecontrol_interaction`: logical, TRUE if used for the SCZD vs control and brain region interaction DE analysis across DLPFC and HIPPO.
* `analysis_eqtl_dlpfc`: logical, TRUE if used for the eQTL analysis with DLFPC samples.
* `analysis_eqtl_hippo`: logical, TRUE if used for the eQTL analysis with HIPPO samples.
* `analysis_eqtl_interaction`: logical, TRUE if used for the eQTL and brain region interaction analysis.

### SNP information

```
     7023861 BrainSeqPhaseII_snp_annotation.txt
```

Columns:

1. `snp`: snp ID
2. `chr_hg38`: chromosome in hg38 coordinates
3. `pos_hg38`: position in hg38 coordinates
4. `chr_hg19`: similar to `chr_hg38` but in hg19 coordinates
5. `pos_hg19`: hg19 position.
6. `cm`
7. `counted`: allele quantified
8. `alt`: alternative allele
9. `type`
10. `newref`: new reference allele (use this)
11. `newcount`: new counted allele (use this)
12. `name`: snp rs number (use this)
13. `rsnumguess`: _guessed_ snp name

### Genotype information

Columns are the subjects labeled by the `BrNum` from the sample metadata. Rows are the SNPs with labels corresponding to the `snp` column from the SNP information table. Entries are 0, 1, 2 or NA.

### Expression feature info

```
      396584 BrainSeqPhaseII_feature_annotation_exon.txt
       24653 BrainSeqPhaseII_feature_annotation_gene.txt
      297182 BrainSeqPhaseII_feature_annotation_jxn.txt
       92733 BrainSeqPhaseII_feature_annotation_tx.txt
```


Each files includes the annotation for each feature type. Since each feature is a bit different, each files has different columns. Some common ones are:

* `seqnames`: chromosome
* `start`: start position of the feature
* `end`: end position
* `strand`: strand of the feature
* `feature_id`: this matches the id used for the eQTL result tables.
* `Symbol` (`gene_name` for `tx`): symbol of the gene the feature corresponds to
* `Class` (except for `tx`): classification for the feature, this really only matters for `jxn` where the exon-exon junction could be un-annotated (novel), in Gencode, etc.
* `gencodeID` (`gene`, `exon`), `gencodeGeneID` (`jxn`), `gene_id` (`tx`): Gencode gene ID.


### eQTL result tables

```
 28863540 BrainSeqPhaseII_eQTL_dlpfc_full.txt
  3698269 BrainSeqPhaseII_eQTL_dlpfc_raggr.txt
 22603132 BrainSeqPhaseII_eQTL_hippo_full.txt
  3698269 BrainSeqPhaseII_eQTL_hippo_raggr.txt
 22603132 BrainSeqPhaseII_eQTL_interaction_full.txt
```

Each of these tables has the following columns:

1. `snp`: snp ID, matches the `snp` column from the SNP information file (`BrainSeqPhaseII_snp_annotation.txt`).
2. `feature_id`: feature ID, matches the `feature_id` column of the expression feature files.
3. `statistic`: eQTL statistic
4. `pvalue`: nominal p-value
5. `FDR`: FDR adjusted p-value
6. `beta`
7. `Type`: feature type in lower-case. This matches the expression feature files, for example `gene` for `BrainSeqPhaseII_feature_annotation_gene.txt`.

The 5 tables correspond to the following types of eQTLs:

* DLPFC using all the genome (`dlpfc_full`), so the samples with `TRUE` under the `analysis_eqtl_dlpfc` column from the sample info table. Shows only results with a nominal p-value <0.001.
* HIPPO using all the genome (`hippo_full`), so the samples with `TRUE` under the `analysis_eqtl_hippo` column from the sample info table. Shows only results with a nominal p-value <0.001.
* Brain region and eQTL interaction across all the genome (`interation_full`), so the samples with `TRUE` under the `analysis_eqtl_interaction` column from the sample info table. Shows only results with a nominal p-value <0.001.
* Sub-analysis using only the PGC2 SNPs and neighboring SNPs identified with rAggr (9,736 SNPs) for DLFPC (`dlpfc_raggr`). Given the smaller number of SNPs considered, the FDR is different compared to `dlpfc_full`. Shows all associations (so no nominal p-value filter).
* Sub-analysis using only the PGC2 SNPs and neighboring SNPs identified with rAggr (9,736 SNPs) for HIPPO (`hippo_raggr`). Given the smaller number of SNPs considered, the FDR is different compared to `hippo_full`. Shows all associations (so no nominal p-value filter).

We considered as significant associations those that had a FDR < 0.01.

#### GTEx eQTL replication tables

We also have 3 more tables that include all associations (no nominal p-value filter) between any of the features that had a significant eQTL in the 5 main analyses against any of the SNPs involved in those analyses. Note that some values are `NA` in these tables (likely to them not being observed in the GTEx samples).

```
758634578 BrainSeqPhaseII_eQTL_dlpfc_replication_GTEx.txt
758634578 BrainSeqPhaseII_eQTL_hippo_replication_GTEx.txt
398556091 BrainSeqPhaseII_eQTL_interaction_replication_GTEx.txt
```

The nominal p-value can be reported for replication purposes.

### Cleaned (_residualized_) expression tables

The tables contained the expression values for making boxplots or other types of graphs. The expression values are already normalized (RPKM: gene/exon, RP10M: jxn, TPM: tx) and log2 scaled (log2(x + 0.5)) with covariates for the corresponding analysis removed. Samples can be identified using the `analysis_*` columns from the sample table. The column names are the `RNum` from the sample metadata table and the rows correspond to the `feature_id` (as in the same order as the feature annotation tables).

#### HIPPO or DLPFC eQTLs

Can be visualized as 3 boxplots, one per allele type, for example CC, AC and AA.

[Sample visualization](../eQTL_full/hippo_top_eqtl_PGC_indexSNPs.pdf)


```
396584 BrainSeqPhaseII_clean_expression_eqtl_dlpfc_exon.txt
 24653 BrainSeqPhaseII_clean_expression_eqtl_dlpfc_gene.txt
297182 BrainSeqPhaseII_clean_expression_eqtl_dlpfc_jxn.txt
 92733 BrainSeqPhaseII_clean_expression_eqtl_dlpfc_tx.txt
396584 BrainSeqPhaseII_clean_expression_eqtl_hippo_exon.txt
 24653 BrainSeqPhaseII_clean_expression_eqtl_hippo_gene.txt
297182 BrainSeqPhaseII_clean_expression_eqtl_hippo_jxn.txt
 92733 BrainSeqPhaseII_clean_expression_eqtl_hippo_tx.txt
```

#### DLPFC and HIPPO interaction eQTLs

Can be visualized as 6 (3 x 2) boxplots, one per allele type, for example CC, AC and AA for each brain region. Brain region can be extracted from the sample metadata column.

[Sample visualization](../eQTL_full/interaction_top_eqtl_PGC_indexSNPs.pdf)


```
396584 BrainSeqPhaseII_clean_expression_eqtl_interaction_exon.txt
 24653 BrainSeqPhaseII_clean_expression_eqtl_interaction_gene.txt
297182 BrainSeqPhaseII_clean_expression_eqtl_interaction_jxn.txt
 92733 BrainSeqPhaseII_clean_expression_eqtl_interaction_tx.txt
```

#### DLPFC vs HIPPO differential expression (adult, prenatal)

Can be visualized as 2 boxplots, one per brain region. [R colors](http://bxhorn.com/r-color-tables/):

* DLPFC: `'dark orange'`
* HIPPO: `'skyblue3'`

[Sample visualization](../region_specific/pdf/tophits_adult_gene_cleaned.pdf)

```
396584 BrainSeqPhaseII_clean_expression_regionspecific_adult_exon.txt
 24653 BrainSeqPhaseII_clean_expression_regionspecific_adult_gene.txt
297182 BrainSeqPhaseII_clean_expression_regionspecific_adult_jxn.txt
 92733 BrainSeqPhaseII_clean_expression_regionspecific_adult_tx.txt
396584 BrainSeqPhaseII_clean_expression_regionspecific_prenatal_exon.txt
 24653 BrainSeqPhaseII_clean_expression_regionspecific_prenatal_gene.txt
297182 BrainSeqPhaseII_clean_expression_regionspecific_prenatal_jxn.txt
 92733 BrainSeqPhaseII_clean_expression_regionspecific_prenatal_tx.txt
```

Note that the analysis was done using limma-voom with log2(CPM + 0.5) instead of log2(RPKM + 1) for genes/exons or log2(RP10M + 1) for jxns.

#### Development differential expression

Can be visualized as a scatterplot of expression versus time split by the age linear spline cutoffs (age in years: 0, 1, 10, 20 and 50). Age below 0 (prenatal) can be transformed from years to the PCW scale using the function `round(range(age) * 52 + 40, 0)`. [R colors](http://bxhorn.com/r-color-tables/):

* DLPFC: `'dark orange'`
* HIPPO: `'skyblue3'`

[Sample visualization](../development/pdf/tophits_gene_cleaned.pdf)

```
396584 BrainSeqPhaseII_clean_expression_development_exon.txt
 24653 BrainSeqPhaseII_clean_expression_development_gene.txt
297182 BrainSeqPhaseII_clean_expression_development_jxn.txt
 92733 BrainSeqPhaseII_clean_expression_development_tx.txt
```

Note that the analysis was done using limma-voom with log2(CPM + 0.5) instead of log2(RPKM + 1) for genes/exons or log2(RP10M + 1) for jxns.

#### SCZD vs control differential expression (DLPFC, HIPPO)

Can be visualized as 2 boxplots, one for SCZD cases and one for controls. [R colors](http://bxhorn.com/r-color-tables/):

* SCZD case: `'aquamarine4'`
* control: `'orchid4'`

[Sample visualization](https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/pdf/top_50each_DLPFC_cleaned.pdf)

```
396584 BrainSeqPhaseII_clean_expression_sczd_casecontrol_dlpfc_exon.txt
 24653 BrainSeqPhaseII_clean_expression_sczd_casecontrol_dlpfc_gene.txt
297182 BrainSeqPhaseII_clean_expression_sczd_casecontrol_dlpfc_jxn.txt
 92733 BrainSeqPhaseII_clean_expression_sczd_casecontrol_dlpfc_tx.txt
396584 BrainSeqPhaseII_clean_expression_sczd_casecontrol_hippo_exon.txt
 24653 BrainSeqPhaseII_clean_expression_sczd_casecontrol_hippo_gene.txt
297182 BrainSeqPhaseII_clean_expression_sczd_casecontrol_hippo_jxn.txt
 92733 BrainSeqPhaseII_clean_expression_sczd_casecontrol_hippo_tx.txt
```

Note that the analysis was done using limma-voom with log2(CPM + 0.5) instead of log2(RPKM + 1) for genes/exons or log2(RP10M + 1) for jxns.

#### SCZD vs control and brain region interaction differential expression

Can be visualized as 4 boxplots, one for SCZD cases and one for controls for each brain region. [R colors](http://bxhorn.com/r-color-tables/):

* SCZD case: `'aquamarine4'` (box)
* control: `'orchid4'` (box)
* DLPFC: `'dark orange'` (dots)
* HIPPO: `'skyblue3'` (dots)

[Sample visualization](../casecontrolint/pdf/de_gene_cleaned.pdf)

```
396584 BrainSeqPhaseII_clean_expression_sczd_casecontrol_interaction_exon.txt
 24653 BrainSeqPhaseII_clean_expression_sczd_casecontrol_interaction_gene.txt
297182 BrainSeqPhaseII_clean_expression_sczd_casecontrol_interaction_jxn.txt
 92733 BrainSeqPhaseII_clean_expression_sczd_casecontrol_interaction_tx.txt
```

Note that the analysis was done using limma-voom with log2(CPM + 0.5) instead of log2(RPKM + 1) for genes/exons or log2(RP10M + 1) for jxns.

### Correlation

The `DLFPC_HIPPO_correlation_*` values from the sample metadata table (remember that they are only listed for the `DLPFC` samples, although the values for the `HIPPO` ones would be identical for the 265 subjects used in this analysis) can be visualized as 2 boxplots, one for SCZD cases and one for controls. [R colors](http://bxhorn.com/r-color-tables/):

* SCZD case: `'aquamarine4'`
* control: `'orchid4'`

[Sample visualization](../correlation/pdf/indv_box_sczd.pdf) (although this one doesn't have the right colors)
