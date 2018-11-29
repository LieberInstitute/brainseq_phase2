This directory contains the code used to make the supplementary tables of the BrainSeq Phase II manuscript. The tables are available at bioRxiv at https://www.biorxiv.org/content/early/2018/09/26/426213. More specifically at https://www.biorxiv.org/content/early/2018/09/26/426213.figures-only.


Below we describe the columns for some of these tables.

## Supplementary Table 2

### SCZD

#### gene

1. X (rowname): GENCODE gene id
2. Length: Coding length
3. gencodeID: GENCODE gene id
4. ensemblID: ENSEMBL gene id
5. gene_type: protein coding, etc
6. Symbol: Gene symbol
7. EntrezID: Gene ENTREZ id
8. Class: all are in the annotation (this is more useful for other features)
9. meanExprs: mean RPKM expression across all 900 BSP2 samples
10. NumTx: number of transcripts
11. gencodeTx: GENCODE transcript ids
12. passExprsCut: whether the gene passed the gene expression cutoff (true for all cases in this table)
13. logFC: log2 fold change between SCZD cases and non-psychiatric controls
14. AveExpr: average normalized expression calculated by `limma`.
15. t: t-statistic SCZD vs control.
16. P.Value: p-value of the t-statistic
17. adj.P.Val: FDR adjusted p-value
18. B: produced by `limma`
19. region: brain region, either DLPFC or SCZD
20. type: gene


#### exon

1. X (rowname): LIBD internal exon id.

#### jxn

1. X (rowname): chromsome, start and end positions of the junction and strand
2. inGencode: is it annotated in Gencode v25?
3. inGencodeStart: is the start of the junction annotated in Gencode v25?
4. inGencodeEnd: is the end  of the junction annotated in Gencode v25?
5. gencodeGeneID: Gencode v25 gene id
6. ensemblID: ENSEMBL gene id
7. Symbol: gene symbol
8. gencodeStrand: Annotated strand of the junction
9. gencodeTx: Gencode transcript ids
10. numTx: number of transcripts
11. Class: annotated, novel, etc
12. startExon: exon start
13. endExon: exon end
14. newGeneID: Gencode gene id assigned to the junction
15. newGeneSymbol: gene symbol(s) assigned to the junction (could be a fusion)
16. isFusion: is this a fusion of annotated genes?
17. meanExprs: mean RP10M expression of the junction in all 900 BSP2 samples
18. Length: 100 (set this way to calculate RP10M using the recount Bioconductor package)
19. passExprsCut: does the junction pass the expression cutoff?
20. logFC: (see gene table description)
21. AveExpr:
22. t:
23. P.Value:
24. adj.P.Val:
25. B:
26. region:
27. type: jxn

#### tx

1. X (rowname): Gencode transcript id.
2. logFC: (see gene table description)
3. AveExpr:
4. t:
5. P.Value:
6. adj.P.Val:
7. B: 
8. source: transcript source
9. type: tx
10. score:
11. phase:
12. gene_id:
13. gene_type:
14. gene_status:
15. gene_name:
16. level:
17. havana_gene:
18. transcript_id:
19. transcript_type:
20. transcript_status:
21. transcript_name:
22. transcript_support_level:
23. tag:
24. havana_transcript:
25. exon_number:
26. exon_id:
27. ont:
28. protein_id:
29. ccdsid:
30. meanExprs:
31. passExprsCut:
32. ensemblID:
33. region:

### Development

#### gene

1. X (rowname): gene."Gencode gene ID"
2. Age.RegionHIPPO: coefficient beta 13 from equation 3 (see the pre-print) which tests for an interaction between age (the first linear spline) and brain region (uses DLPFC as the reference and HIPPO as the contrast)
3. RegionHIPPO.fetal: coefficient beta 14 from equation 3
4. RegionHIPPO.birth: coefficient beta 15 from equation 3
5. RegionHIPPO.infant: coefficient beta 16 from equation 3
6. RegionHIPPO.child: coefficient beta 17 from equation 3
7. RegionHIPPO.teen: coefficient beta 18 from equation 3
8. RegionHIPPO.adult: coefficient beta 19 from equation 3
9. AveExpr: calculated by `limma` for the model in equation 3
10. F: F-statistic testing whether any of the coefficients from beta 13 to 19 in equation 3 are significantly different from 0
11. P.Value: p-value for the F-statistic
12. adj.P.Val: FDR adjusted p-value
13. type: gene
14. P.Bonf: bonferroni adjusted p-value (the one used in this analysis)
15. span_Age.RegionHIPPO: coefficient beta 13 from equation 3 applied to the BrainSpan data
16. span_RegionHIPPO.fetal: coefficient beta 14 from equation 3 applied to the BrainSpan data
17. span_RegionHIPPO.birth: coefficient beta 15 from equation 3 applied to the BrainSpan data
18. span_RegionHIPPO.infant: coefficient beta 16 from equation 3 applied to the BrainSpan data
19. span_RegionHIPPO.child: coefficient beta 17 from equation 3 applied to the BrainSpan data
20. span_RegionHIPPO.teen: coefficient beta 18 from equation 3 applied to the BrainSpan data
21. span_AveExpr: calculated by `limma` for the model in equation 3 when using the BrainSpan data
22. span_F: F-statistic when using the BrainSpan data
23. span_P.Value: p-value when using the BrainSpan data
24. span_adj.P.Val: FDR adjusted p-value when using the BrainSpan data
25. span_type: gene
26. span_P.Bonf: bonferroni adjusted p-value when using the BrainSpan data
27. seqnames: chromosome name
28. start: gene start
29. end: gene end
30. width: gene length (total, includes coding and non-coding)
31. strand: gene strand
32. Length: gene coding length
33. gencodeID: Gencode gene id
34. ensemblID: ENSEMBL gene id
35. gene_type: protein coding, etc
36. Symbol: gene symbol
37. EntrezID: gene Entrez id
38. Class: InGen since all are in Gencode (more useful for other features)
39. meanExprs: mean RPKM expression across all 900 BSP2 samples
40. NumTx: number of transcripts
41. gencodeTx: Gencode transcript ids
42. passExprsCut: whether it passes the expression cutoff (true for all cases in this table)
43. replicates_in_BrainSpan: whether the result replicates in BrainSpan

#### exon

Similar differences as in the gene and exon SCZD DE tables.

#### jxn

Similar differences as in the gene and jxn SCZD DE tables.

#### tx

Similar differences as in the gene and tx SCZD DE tables.


### Region-specific

#### gene

1. X (rowname): "age group (adult or prenatal)"_gene"Gencode gene id"
2. logFC: log2 fold change between HIPPO and DLPFC.
3. AveExpr: calculated by `limma` for the model in equation 2 (see pre-print)
4. t: t-statistic
5. P.Value: p-value
6. adj.P.Val: FDR adjusted p-value
7. B: calculated by `limma`
8. age: age group, either adult or prenatal
9. type: gene
10. P.Bonf: bonferroni adjusted p-value (the one used in the analysis)
11. span_logFC: similar as above, but for BrainSpan data
12. span_AveExpr:
13. span_t:
14. span_P.Value:
15. span_adj.P.Val:
16. span_B:
17. span_age:
18. span_type:
19. span_P.Bonf:
20. seqnames: chromosome
21. start: gene start
22. end: gene end
23. width: gene total length
24. strand: gene strand
25. Length: gene coding length
26. gencodeID: Gencode gene id
27. ensemblID: ENSEMBL gene id
28. gene_type: protein coding, etc
29. Symbol: gene symbol
30. EntrezID: gene ENTREZ id
31. Class: InGen since they are all in Gencode (more useful for other features)
32. meanExprs: mean RPKM expression across all 900 BSP2 samples
33. NumTx: number of transcripts
34. gencodeTx: Gencode transcript ids
35. passExprsCut: whether it passes the expression cutoff (true for all cases in this case)
36. replicates_in_BrainSpan: whether the result replicates in BrainSpan

#### exon

Similar differences as in the gene and exon SCZD DE tables.

#### jxn

Similar differences as in the gene and jxn SCZD DE tables.

#### tx

Similar differences as in the gene and tx SCZD DE tables.



